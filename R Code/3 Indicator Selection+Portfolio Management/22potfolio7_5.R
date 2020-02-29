library(xts) 
library(quantmod) 
library(PerformanceAnalytics)



########### Get Data ########################
setwd("C:/Users/m8sha/Desktop/DATA/Strategy")

X_lin <- read.csv('Strategy.Dailylinreturn.txt',header = TRUE,sep = '')
X_log <- read.csv('Strategy.Dailylogreturn.txt',header = TRUE,sep = '')
rownames(X_lin) <- substr(X_lin[,1],1,10)
rownames(X_log) <- substr(X_log[,1],1,10)
X_lin <- X_lin[,-1]
X_log <- X_log[,-1]
X_lin <- as.xts(X_lin)
X_log <- as.xts(X_log)


### Set portfolio Funtions ##############################
# Uniform, GMVP, Markowitz, MaxSR, DR, CVaR, MAX-DD,Ave-DD,CDaR portfolio 

# GMVP #
library(CVXR)
portolioGMVP <- function(Sigma) {
  w <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1)) 
  result <- solve(prob)
  return(as.vector(result$getValue(w))) }

# Markowitz #
portolioMarkowitz <- function(mu, Sigma, lmd = 0.5) {
  w <- Variable(nrow(Sigma))
  prob <- Problem(Maximize(t(mu) %*% w - lmd*quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1)) 
  result <- solve(prob)
  return(as.vector(result$getValue(w))) }
# MaxSR #
portolioMaxSharpeRatio <- function(mu, Sigma) { 
  w_ <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w_, Sigma)),
                  constraints = list(w_ >= 0, t(mu) %*% w_ == 1)) 
  result <- solve(prob)
  return(as.vector(result$getValue(w_)/sum(result$getValue(w_)))) }

# Downside Risk portfolio
portfolioDR <- function(X, lmd = 0.5, alpha) {
  T <- nrow(X)  
  N <- ncol(X)  
  X <- as.matrix(X)  
  mu <- colMeans(X)  
  w <- Variable(N)
  prob <- Problem(Maximize(t(w) %*% mu - (lmd/T) * sum(pos(t(mu) %*% w - X %*% w))^alpha), 
                  constraints = list(w >= 0, sum(w) == 1))  
  result <- solve(prob) 
  return(as.vector(result$getValue(w))) }

# CVaR portfolio 
portolioCVaR <- function(X, lmd = 0.5, alpha) {
  T <- nrow(X)  
  N <- ncol(X)  
  X <- as.matrix(X)  
  mu <- colMeans(X) # variables  
  w <- Variable(N)  
  z <- Variable(T)  
  zeta <- Variable(1) # problem  
  prob <- Problem(Maximize(t(w) %*% mu - lmd*zeta - (lmd/(T*(1-alpha))) * sum(z)), 
                  constraints = list(z >= 0, z >= -X %*% w - zeta,w >= 0, sum(w) == 1))  
  result <- solve(prob) 
  return(as.vector(result$getValue(w))) }


# Max-Drawdown portfolio 
portfolioMaxDD <- function(X, c) {
  T <- nrow(X)  
  N <- ncol(X)  
  X <- as.matrix(X)  
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)  
  mu <- colMeans(X) # variables  
  w <- Variable(N)  
  u <- Variable(T) # problem  
  prob <- Problem(Maximize(t(w) %*% mu), 
                  constraints = list(w >= 0, sum(w) == 1,u <= X_cum %*% w + c,
                                     u >= X_cum %*% w,u[-1] >= u[-T]))  
  result <- solve(prob) 
  return(as.vector(result$getValue(w))) } 


# Ave-DD portfolio 
portfolioAveDD <- function(X, c) {
  T <- nrow(X)  
  N <- ncol(X)  
  X <- as.matrix(X)  
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)  
  mu <- colMeans(X) # variables  
  w <- Variable(N)  
  u <- Variable(T) # problem  
  prob <- Problem(Maximize(t(w) %*% mu), 
                  constraints = list(w >= 0, sum(w) == 1, mean(u) <= mean(X_cum %*% w) + c,
                                     u >= X_cum %*% w,u[-1] >= u[-T]))  
  result <- solve(prob) 
  return(as.vector(result$getValue(w))) } 

# CDaR portfolio 
portfolioCDaR <- function(X, c, alpha) {
  T <- nrow(X)  
  N <- ncol(X)  
  X <- as.matrix(X)  
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)  
  mu <- colMeans(X) # variables  
  w <- Variable(N)  
  z <- Variable(T)  
  zeta <- Variable(1)  
  u <- Variable(T) # problem  
  prob <- Problem(Maximize(t(w) %*% mu), constraints = list(w >= 0, sum(w) == 1,zeta + (1/(T*(1-alpha))) * sum(z) <= c,
                                                            z >= 0, z >= u - X_cum %*% w - zeta,u >= X_cum %*% w,
                                                            u[-1] >= u[-T]))  
  result <- solve(prob) 
  return(as.vector(result$getValue(w))) } 

# inverse
portolioIV<-function(Sigma,mu=0){
  result<-(1/diag(Sigma))/sum( 1/diag(Sigma))
  return(as.vector(result))
}




### Set function for returns ###
CalculateWeights <- function(strategies_log_roll,N){
  
  # Mu & Sigma
  mu_roll <- colMeans(strategies_log_roll)
  Sigma_roll <- cov(strategies_log_roll)
  
  ### Uniform Portfolio ###
  w_unif <- rep(1/N, N)
  
  ### GMVP portfolio ###
  w_GMVP <- portolioGMVP(Sigma_roll)
  
  ### Markowitz Portfolio ###
  w_Markowitz <- portolioMarkowitz(mu_roll, Sigma_roll, lmd = 2)
  
  ### Maximum Sharpe Ratio Portfolio ###
  
  w_maxSR <- portolioMaxSharpeRatio(mu_roll, Sigma_roll)
  print(length(w_maxSR))
  ### Long or Long-Short quintile portfolios ###
  
  # find indices of sorted stocks
  i1 <- sort(mu_roll, decreasing = TRUE, index.return = TRUE)$ix
  i2 <- sort(mu_roll/diag(Sigma_roll), decreasing = TRUE, index.return = TRUE)$ix
  i3 <- sort(mu_roll/sqrt(diag(Sigma_roll)), decreasing = TRUE, index.return = TRUE)$ix
  # create portfolios
  w_Lquintile1 <- w_Lquintile2 <- w_Lquintile3 <- rep(0, N) 
  w_Lquintile1[i1[1:round(N/5)]] <- 1/round(N/5) 
  w_Lquintile2[i2[1:round(N/5)]] <- 1/round(N/5) 
  w_Lquintile3[i3[1:round(N/5)]] <- 1/round(N/5)
  # w_Lquintile <- cbind(w_Lquintile1, w_Lquintile2, w_Lquintile3)
  # rownames(w_Lquintile) <- colnames(strategies_log_roll)
  # colnames(w_Lquintile) <- c("Lquintile1", "Lquintile2", "Lquintile3")
  #w_Lquintile
  
  ## Downside Risk portfolio
  w_DR_alpha1 <- portfolioDR(strategies_log_roll, alpha = 1)
  w_DR_alpha2 <- portfolioDR(strategies_log_roll, alpha = 2)
  w_DR_alpha3 <- portfolioDR(strategies_log_roll, alpha = 3)
  print(length(w_DR_alpha3))
  # ### CVaR
  w_CVaR095 <- portolioCVaR(strategies_log_roll, alpha = 0.95)
  w_CVaR099 <- portolioCVaR(strategies_log_roll, alpha = 0.99)
  print(length(w_CVaR095))
  # ### MAX-DD
  w_MaxDD_c018 <- portfolioMaxDD(strategies_log_roll, c = 0.18)
  w_MaxDD_c021 <- portfolioMaxDD(strategies_log_roll, c = 0.21)
  w_MaxDD_c024 <- portfolioMaxDD(strategies_log_roll, c = 0.24)
  
  # ### Avg-DD
  w_AveDD_c004 <- portfolioAveDD(strategies_log_roll, c = 0.04)
  w_AveDD_c006 <- portfolioAveDD(strategies_log_roll, c = 0.06)
  w_AveDD_c008 <- portfolioAveDD(strategies_log_roll, c = 0.08)
  
  # ### CDaR
  w_CDaR095_c014 <- portfolioCDaR(strategies_log_roll, c = 0.14, alpha = 0.95)
  w_CDaR099_c016 <- portfolioCDaR(strategies_log_roll, c = 0.16, alpha = 0.99)
  w_CDaR095_c016 <- portfolioCDaR(strategies_log_roll, c = 0.16, alpha = 0.95)
  w_CDaR099_c018 <- portfolioCDaR(strategies_log_roll, c = 0.18, alpha = 0.99)
  #
  ### Return-risk tradeoff for all portfolios ###
  # put together all portfolios
  # w_all <- cbind(w_unif, w_GMVP, w_Markowitz, w_maxSR,
  #                w_Lquintile3,
  #                w_DR_alpha1,
  #                w_CVaR099,
  #                w_MaxDD_c018,
  #                w_AveDD_c004)
  w_all <- cbind(w_unif, w_GMVP, w_Markowitz, w_maxSR,
                 #w_Lquintile1,w_Lquintile2,w_Lquintile3,
                 w_DR_alpha1,w_DR_alpha2,w_DR_alpha3,
                 w_CVaR095,w_CVaR099,
                 w_MaxDD_c018,w_MaxDD_c021,w_MaxDD_c024,
                 w_AveDD_c004,w_AveDD_c006,w_AveDD_c008,
                 w_CDaR095_c014,w_CDaR099_c016,w_CDaR095_c016,w_CDaR099_c018)
  
  rownames(w_all) <- colnames(strategies_log_roll)
  colnames(w_all) <- c("unif", "GMVP", "Markowitz", "maxSR",
                       #"Lquintile1","Lquintile2","Lquintile3",
                       "DR_alpha1","DR_alpha2","DR_alpha3",
                       "CVaR095","CVaR099",
                       "MaxDD_c018","MaxDD_c021","MaxDD_c024",
                       "AveDD_c004","AveDD_c006","AveDD_c008",
                       "CDaR095_c014","CDaR099_c016","CDaR095_c016","CDaR099_c018")
  round(w_all, digits = 2)
  return(w_all)
}


########### »» 20 ##############################################
end_20month<-c()
T_trn = 60
t <- nrow(X_log)
X_= X_log[(T_trn+1):t,]
end_week<-endpoints(X_, on = "weeks")
end_month<-endpoints(X_, on = "months")
for (i in end_month){
  ind = max(0,(which(end_week==max(end_week[end_week<=i]))-2))
  end_20month<-c(end_20month,end_week[ind])
}


####################################################################

rolling <- function(X_log,X_lin,on,start){
  if (on == 'weeks'){
    T_trn <- 15}
  else{T_trn <- 60}
  t <- nrow(X_log)
  all_weights = list()
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):t, ]
  
  if (on == 'weeks'||on == 'months'){
    rebal_indices <- T_trn + endpoints(X_log_tst, on = on)}
  else{rebal_indices <- T_trn + c(0,end_20month)}
  
  result =matrix(ncol=22)[-1,]
  
  for (i in 1:(length(rebal_indices)-1)){
    #Training
    if(start==1){X_ <- X_log[(rebal_indices[1]-T_trn+1):rebal_indices[i], ]}
    else{X_ <- X_log[(rebal_indices[i]-T_trn+1):rebal_indices[i], ]}
    #get weights
    w_all <- CalculateWeights(X_,N=ncol(X_log))
    # get test return
    tst <- X_lin[(rebal_indices[i]+1):min(rebal_indices[i+1],nrow(X_log)),] 
    # calculate mean returns as result for each days
    result <- rbind(result,as.matrix(tst)%*%w_all)
    all_weights <- c(all_weights,list(w_all))
  }
  
  rownames(result) = rownames(as.matrix(X_log))[(T_trn+1):t]
  colnames(result) = colnames(w_all)
  
  return(list('result'=result,
              'weight'=all_weights))
}

rollingFri <- function(X_log,X_lin,week,start,T_trn){
  all_weights = list()
  
  t <- nrow(X_log)
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):t, ]
  # X_lin_trn <- X_lin[1:T_trn, ]
  # X_lin_tst <- X_lin[(T_trn+1):t, ]
  # Rank_rolling <- X_log
  # Rank_rolling[] <- NA
  end_Fri<-c()
  end_week<-endpoints(X_log_tst, on = "weeks")
  end_month<-endpoints(X_log_tst, on = "months")
  for (i in end_month){
    if (week=='1'){ind = max(0,(which(end_week==max(end_week[end_week<=i]))-3))}
    if (week=='2'){ind = max(0,(which(end_week==max(end_week[end_week<=i]))-2))}
    else {ind = which(end_week==max(end_week[end_week<=i]))-1} 
    end_Fri<-c(end_Fri,end_week[ind])
  }
  rebal_indices <- T_trn + c(0,end_Fri) 
  #rebal_indices2 <- T_trn + endpoints( X_logContract[243:357,], on = "weeks")
  # index(X_log)[rebal_indices]
  lookback <- T_trn
  result =matrix(ncol=22)[-1,]
  
  for (i in 1:(length(rebal_indices)-1)){
    #Training
    print(i)
    if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
    else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
    #get weights
    w_all <-  CalculateWeights(X_,N=ncol(X_log))
    # get test return
    if(i == (length(rebal_indices)-1)) {tst <- X_lin[(rebal_indices[i]+1):nrow(X_log),] }
    else{tst <- X_lin[(rebal_indices[i]+1):rebal_indices[i+1],] }
    # calculate mean returns as result for each days
    result <- rbind(result,as.matrix(tst)%*%w_all)
    all_weights <- c(all_weights,list(w_all))
  }
  
  rownames(result) = rownames(as.matrix(X_log))[(T_trn+1):t]
  colnames(result) = colnames(w_all)
  
  return(list('result'=result,
              'weight'=all_weights))
}


setwd("C:/Users/m8sha/Desktop/compare/Contracts")
on = "weeks"
week1819start<-rolling(X_log,X_lin,on,start = 1)
week1819block<-rolling(X_log,X_lin,on,start = 0)

on = "months"
month1819start<-rolling(X_log,X_lin,on,start = 1)
#write.csv(month1819start$weight,'month1819start.csv')
month1819block<-rolling(X_log,X_lin,on,start = 0)
#write.csv(month1819block$weight,'month1819block.csv')

on = "changes"
change1819week1start<-rollingFri(X_log,X_lin,week="1",start = 1,T_trn=60)
#write.csv(change1819week1start$weight,'change1819week1start.csv')
change1819week1block<-rollingFri(X_log,X_lin,week="1",start = 0,T_trn=60)
#write.csv(change1819week1block$weight,'change1819week1block.csv')

change1819week2start<-rollingFri(X_log,X_lin,week="2",start = 1,T_trn=60)
#write.csv(change1819week2start$weight,'change1819week2start.csv')
change1819week2block<-rollingFri(X_log,X_lin,week="2",start = 0,T_trn=60)
#write.csv(change1819week2block$weight,'change1819week2block.csv')

change1819week3start<-rollingFri(X_log,X_lin,week="3",start = 1,T_trn=60)
#write.csv(change1819week3start$weight,'change1819week3start.csv')
change1819week3block<-rollingFri(X_log,X_lin,week="3",start = 0,T_trn=60)
write.csv(change1819week3block$weight,'change1819week3block.csv')


layout(matrix(1,1,1))
#par(mfrow=c(3,2))
c_week1819start = chart.CumReturns(week1819start$result, main = "18-19 performance of 1-week fixed start",
                                   wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)
c_week1819block = chart.CumReturns(week1819block$result, main = "18-19 performance of 1-week fixed blocks",
                                   wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_month1819start = chart.CumReturns(month1819start$result, main = "18-19 performance of 1-month fixed start",
                                    wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_month1819block = chart.CumReturns(month1819block$result, main = "18-19 performance of 1-month fixed blocks",
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_change1819week1start = chart.CumReturns(change1819week1start$result, main = "18-19 performance of before week1 fixed start",
                                          wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_change1819week1block = chart.CumReturns(change1819week1block$result, main = "18-19 performance of before week1 fixed blocks",
                                          wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_change1819week2start = chart.CumReturns(change1819week2start$result, main = "18-19 performance of before week2 fixed start",
                                        wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_change1819week2block = chart.CumReturns(change1819week2block$result, main = "18-19 performance of before week2 fixed blocks",
                                          wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_change1819week3start = chart.CumReturns(change1819week3start$result, main = "18-19 performance of before week3 fixed start",
                                          wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_change1819week3block = chart.CumReturns(change1819week3block$result, main = "18-19 performance of before week3 fixed blocks",
                                          wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

layout(matrix(1,1,1))
par(mfrow=c(2,2))
c_week1819start
c_week1819block
c_month1819start
c_month1819block

par(mfrow=c(3,1))
c_change1819week1start
c_change1819week2start
c_change1819week3start

par(mfrow=c(3,1))
c_change1819week1block
c_change1819week2block
c_change1819week3block

######################################### 2019
X_log19tst_week = X_log[229:357,]
X_log19tst_month = X_log[184:357,]
X_lin19tst_week = X_lin[229:357,]
X_lin19tst_month = X_lin[184:357,]


T_trn = 60
t <- nrow(X_log19tst_month)
X_= X_log19tst_month[(T_trn+1):t,]
end_20month<-c()
end_week<-endpoints(X_, on = "weeks")
end_month<-endpoints(X_, on = "months")
for (i in end_month){
  ind = max(0,(which(end_week==max(end_week[end_week<=i]))-3))
  print(ind)
  end_20month<-c(end_20month,end_week[ind])
}


#rolling windows 19 tst
on = "weeks"
week19start<-rolling(X_log19tst_week,X_lin19tst_week,on,start = 1)
week19block<-rolling(X_log19tst_week,X_lin19tst_week,on,start = 0)

#rolling windows 19 tst month
on = "months"
month19start<-rolling(X_log19tst_month,X_lin19tst_month,on,start = 1)
#write.csv(month19start$weight,'month19start.csv')
month19block<-rolling(X_log19tst_month,X_lin19tst_month,on,start = 0)
#write.csv(month19block$weight,'month19block.csv')

on = "changes"
change19week1start<-rollingFri(X_log19tst_month,X_lin19tst_month,week="1",start = 1,T_trn=60)
write.csv(change19week1start$weight,'change19week1start.csv')
change19week1block<-rollingFri(X_log19tst_month,X_lin19tst_month,week="1",start = 0,T_trn=60)
write.csv(change19week1block$weight,'change19week1block.csv')

change19week2start<-rollingFri(X_log19tst_month,X_lin19tst_month,week="2",start = 1,T_trn=60)
write.csv(change19week2start$weight,'change19week2start.csv')
change19week2block<-rollingFri(X_log19tst_month,X_lin19tst_month,week="2",start = 0,T_trn=60)
write.csv(change19week2block$weight,'change19week2block.csv')

change19week3start<-rollingFri(X_log19tst_month,X_lin19tst_month,week="3",start = 1,T_trn=60)
write.csv(change19week3start$weight,'change19week3start.csv')
change19week3block<-rollingFri(X_log19tst_month,X_lin19tst_month,week="3",start = 0,T_trn=60)
write.csv(change19week3block$weight,'change19week3block.csv')


############# 
c_week19start = chart.CumReturns(week19start$result, main = "19 performance of 1-week fixed start",
                                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)
c_week19block = chart.CumReturns(week19block$result, main = "19 performance of 1-week fixed blocks",
                                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_month19start = chart.CumReturns(month19start$result, main = "19 performance of 1-month fixed start",
                                  wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_month19block = chart.CumReturns(month19block$result, main = "19 performance of 1-month fixed blocks",
                                  wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_change19week1start = chart.CumReturns(change19week1start$result, main = "19 performance of before week1 fixed start",
                                        wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_change19week1block = chart.CumReturns(change19week1block$result, main = "19 performance of before week1 fixed blocks",
                                        wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_change19week2start = chart.CumReturns(change19week2start$result, main = "19 performance of before week2 fixed start",
                                        wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_change19week2block = chart.CumReturns(change19week2block$result, main = "19 performance of before week2 fixed blocks",
                                        wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

c_change19week3start = chart.CumReturns(change19week3start$result, main = "19 performance of before week3 fixed start",
                                        wealth.index = T, legend.loc = "topleft", colorset = rainbow12equal)
c_change19week3block = chart.CumReturns(change19week3block$result, main = "19 performance of before week3 fixed blocks",
                                        wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow12equal)

layout(matrix(1,1,1))
par(mfrow=c(2,2))
c_week19start
c_week19block
c_month19start
c_month19block

par(mfrow=c(3,1))
c_change19week1start
c_change19week2start
c_change19week3start

par(mfrow=c(3,1))
c_change19week1block
c_change19week2block
c_change19week3block

####################################

indices <- function(return){
  res<-rbind(table.AnnualizedReturns(return),
             SharpeRatio(return, Rf = 0, p = 0.95, FUN = "StdDev"),
             maxDrawdown(return, weights = NULL, geometric = TRUE, invert = TRUE),
             SortinoRatio(return, MAR=0),
             CDD(return, weights = NULL, geometric = F, invert = TRUE, p = 0.95),
             CDD(return, weights = NULL, geometric = F, invert = TRUE, p = 0.5),
             VaR(return, p = 0.95,method = "historical", portfolio_method ="single",invert = TRUE),
             VaR(return, p = 0.99,method = "historical", portfolio_method ="single",invert = TRUE),
             Return.annualized(return,scale=252)/maxDrawdown(return, weights = NULL, geometric = TRUE, invert = TRUE),
             0.95*Return.annualized(return,scale=252)/CDD(return, weights = NULL, geometric = F, invert = TRUE, p = 0.95),
             0.5*Return.annualized(return,scale=252)/CDD(return, weights = NULL, geometric = F, invert = TRUE, p = 0.5),
             0.95*Return.annualized(return,scale=252)/VaR(return, p = 0.95,method = "historical", portfolio_method ="single",invert = TRUE),
             0.99*Return.annualized(return,scale=252)/VaR(return, p = 0.99,method = "historical", portfolio_method ="single",invert = TRUE))
  return(res)
}

setwd("C:/Users/m8sha/Desktop/compare/Contracts")

write.csv(rbind(t(indices(week1819start$result)),
                t(indices(week1819block$result)),
                t(indices(month1819start$result)),
                t(indices(month1819block$result)),
                t(indices(change1819week1start$result)),
                t(indices(change1819week1block$result)),
                t(indices(change1819week2start$result)),
                t(indices(change1819week2block$result)),
                t(indices(change1819week3start$result)),
                t(indices(change1819week3block$result)),
                t(indices(week19start$result)),
                t(indices(week19block$result)),
                t(indices(month19start$result)),
                t(indices(month19block$result)),
                t(indices(change19week1start$result)),
                t(indices(change19week1block$result)),
                t(indices(change19week2start$result)),
                t(indices(change19week2block$result)),
                t(indices(change19week3start$result)),
                t(indices(change19week3block$result))),"indices20.csv")






