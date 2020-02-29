library(xts) 
library(quantmod) 
library(PerformanceAnalytics)



########### Get Data ########################
setwd("C:/Users/m8sha/Desktop/DATA/Strategy")

strategies_linret <- read.csv('Strategy.Dailylinreturn.txt',header = TRUE,sep = '')
strategies_logret <- read.csv('Strategy.Dailylogreturn.txt',header = TRUE,sep = '')
rownames(strategies_linret) <- substr(strategies_linret[,1],1,10)
rownames(strategies_logret) <- substr(strategies_logret[,1],1,10)
strategies_linret <- strategies_linret[,-1]
strategies_logret <- strategies_logret[,-1]
strategies_linret <- as.xts(strategies_linret)
strategies_logret <- as.xts(strategies_logret)
str(strategies_linret)
str(strategies_logret)
# head(strategies_linret)
#head(strategies_logret)
#tail(strategies_linret)
#tail(strategies_logret)



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
CalculateTestReturns <- function(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll){
  
  strategies_log_rolltrn<- strategies_log_roll[1:T_trn_roll, ]
  strategies_log_rolltst<- strategies_log_roll[(T_trn_roll+1):T_roll, ]
  strategies_lin_rolltrn<- strategies_lin_roll[1:T_trn_roll, ]
  strategies_lin_rolltst<- strategies_lin_roll[(T_trn_roll+1):T_roll, ]
  
  # Mu & Sigma
  mu_roll <- colMeans(strategies_log_rolltrn)
  Sigma_roll <- cov(strategies_log_rolltrn)
  
  ### Uniform Portfolio ###
  w_unif <- rep(1/N, N)
  
  ### GMVP portfolio ###
  w_GMVP <- portolioGMVP(Sigma_roll)
  
  ### Markowitz Portfolio ###
  w_Markowitz <- portolioMarkowitz(mu_roll, Sigma_roll, lmd = 2)
  
  ### Maximum Sharpe Ratio Portfolio ###
  
  w_maxSR <- portolioMaxSharpeRatio(mu_roll, Sigma_roll)
  
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
  w_DR_alpha1 <- portfolioDR(strategies_log_rolltrn, alpha = 1)
  w_DR_alpha2 <- portfolioDR(strategies_log_rolltrn, alpha = 2)
  w_DR_alpha3 <- portfolioDR(strategies_log_rolltrn, alpha = 3)

  # ### CVaR
  w_CVaR095 <- portolioCVaR(strategies_log_rolltrn, alpha = 0.95)
  w_CVaR099 <- portolioCVaR(strategies_log_rolltrn, alpha = 0.99)

  # ### MAX-DD
  w_MaxDD_c018 <- portfolioMaxDD(strategies_log_rolltrn, c = 0.18)
  w_MaxDD_c021 <- portfolioMaxDD(strategies_log_rolltrn, c = 0.21)
  w_MaxDD_c024 <- portfolioMaxDD(strategies_log_rolltrn, c = 0.24)

  # ### Avg-DD
  w_AveDD_c004 <- portfolioAveDD(strategies_log_rolltrn, c = 0.04)
  w_AveDD_c006 <- portfolioAveDD(strategies_log_rolltrn, c = 0.06)
  w_AveDD_c008 <- portfolioAveDD(strategies_log_rolltrn, c = 0.08)

  # ### CDaR
  w_CDaR095_c014 <- portfolioCDaR(strategies_log_rolltrn, c = 0.14, alpha = 0.95)
  w_CDaR099_c016 <- portfolioCDaR(strategies_log_rolltrn, c = 0.16, alpha = 0.99)
  w_CDaR095_c016 <- portfolioCDaR(strategies_log_rolltrn, c = 0.16, alpha = 0.95)
  w_CDaR099_c018 <- portfolioCDaR(strategies_log_rolltrn, c = 0.18, alpha = 0.99)
  #
  ### Return-risk tradeoff for all portfolios ###
  # put together all portfolios
  w_all <- cbind(w_unif, w_GMVP, w_Markowitz, w_maxSR,
                 w_Lquintile3,
                 w_DR_alpha1,
                 w_CVaR099,
                 w_MaxDD_c018,
                 w_AveDD_c004)
  # w_all <- cbind(w_unif, w_GMVP, w_Markowitz, w_maxSR,
  #                w_Lquintile1,w_Lquintile2,w_Lquintile3,
  #                w_DR_alpha1,w_DR_alpha2,w_DR_alpha3,
  #                w_CVaR095,w_CVaR099,
  #                w_MaxDD_c018,w_MaxDD_c021,w_MaxDD_c024,
  #                w_AveDD_c004,w_AveDD_c006,w_AveDD_c008,
  #                w_CDaR095_c014,w_CDaR099_c016,w_CDaR095_c016,w_CDaR099_c018)

  rownames(w_all) <- colnames(strategies_log_roll)
  colnames(w_all) <- c("uniform", "GMVP", "Markowitz", "maxSR",
                       'Lquintile3',
                       'DR_alpha1',
                       'CVaR099',
                       'MaxDD_c018',
                       'AveDD_c004')
  round(w_all, digits = 2)
  
  ### Then compute the performance ###
  
  # compute returns of all portfolios
  ret_all <- xts(strategies_lin_roll %*% w_all, order.by=index(strategies_lin_roll))
  ret_all_trn <- ret_all[1:T_trn_roll, ]
  result <- ret_all[-c(1:T_trn_roll), ]
  return(result)
  # return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
  #             "MDD"=maxDrawdown(result, weights = NULL, geometric = TRUE, invert = TRUE),
  #             'Calmer'=CalmarRatio(result),
  #             'Annulized returns'=Return.annualized(result),
  #             # 'CVaR'=ETL(result, p = 0.95,method = c("modified")),
  #             "result" = result))
                     
 }

#-----------------------------------------------------------------------------------

#### Rolling ####
N <- ncol(X_log)
T_trn <- 15



### 1. For weeks ###

## 1.1 15-5 days test ##

# 1.1.1 fixed block #
N <- ncol(X_log)
T_roll <- 20
T_trn_roll <- 15

strategies_log_roll_1st <- strategies_logret[(1+(1-1)*5):(20+(1-1)*5),]
strategies_lin_roll_1st <- strategies_linret[(1+(1-1)*5):(20+(1-1)*5),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- (floor(nrow(strategies_logret)/5)-3)

for (j in 2:last_test){
  strategies_log_roll <- strategies_logret[(1+(j-1)*5):(20+(j-1)*5),]
  strategies_lin_roll <- strategies_linret[(1+(j-1)*5):(20+(j-1)*5),]
  ret_rolling <- rbind(ret_rolling,CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll))
}
r_5dfb <- ret_rolling

# performance--------
layout(matrix(1,1,1))
# all years #
c_5dfb_1819_1p <-chart.CumReturns(r_5dfb, main = "Performance of weeks fixed blocks rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
# 2019 first 6 months #
r_5dfb_19 <- r_5dfb[which(index(r_5dfb)>='2019-01-02 CST')]
c_5dfb_19_1p <- { chart.CumReturns(r_5dfb_19, main = "Performance of weeks fixed blocks rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
}


# 1.1.2 fixed start #

strategies_log_roll_1st <- strategies_logret[(1+(1-1)*5):(20+(1-1)*5),]
strategies_lin_roll_1st <- strategies_linret[(1+(1-1)*5):(20+(1-1)*5),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,20,15)
last_test <- (floor(nrow(strategies_logret)/5)-3)

for (j in 2:last_test){
  strategies_log_roll <- strategies_logret[1:(20+(j-1)*5),]
  strategies_lin_roll <- strategies_linret[1:(20+(j-1)*5),]
  T_roll <- (j-1)*5+20
  T_trn_roll <- (j-1)*5+15
  ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}
r_5dfs <- ret_rolling

# performance--------

# all years #
c_5dfs_1819_1p <-{ chart.CumReturns(r_5dfs, main = "Performance of weeks fixed start rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
}


# 2019 first 6 months #

r_5dfs_19 <- r_5dfs[which(index(r_5dfs)>='2019-01-02 CST')]
c_5dfs_19_1p <- { chart.CumReturns(r_5dfs_19, main = "Performance of weeks fixed start rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
}

### 2. for months ###

## 2.1 40-20 days test 

# 2.1.1 fixed block #
T_roll <- 60
T_trn_roll <- 40

# For first 20 d #
strategies_log_roll_1st <- strategies_logret[(1+(1-1)*20):(60+(1-1)*20),]
strategies_lin_roll_1st <- strategies_linret[(1+(1-1)*20):(60+(1-1)*20),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,60,40)
last_test <- (floor(nrow(strategies_logret)/20)-3)

for (j in 2:last_test){
  strategies_log_roll <- strategies_logret[(1+(j-1)*20):(60+(j-1)*20),]
  strategies_lin_roll <- strategies_linret[(1+(j-1)*20):(60+(j-1)*20),]
  ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

r_20dfb <- ret_rolling

# performance---------

# all years #
c_20dfb_1819_1p <-{ chart.CumReturns(r_20dfb, main = "Performance of months fixed blocks rolling", 
                                     wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
}


# 2019 first 6 months #

r_20dfb_19 <- r_20dfb[which(index(r_20dfb)>='2019-01-02 CST')]
c_20dfb_19_1p<-{ chart.CumReturns(r_20dfb_19, main = "Performance of months fixed blocks rolling", 
                                  wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
}


# 2.1.2  fixed start #

strategies_log_roll_1st <- strategies_logret[(1+(1-1)*20):(60+(1-1)*20),]
strategies_lin_roll_1st <- strategies_linret[(1+(1-1)*20):(60+(1-1)*20),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,60,40)
last_test <- (floor(nrow(strategies_logret)/20)-3)

for (j in 2:last_test){
  strategies_log_roll <- strategies_logret[1:(60+(j-1)*20),]
  strategies_lin_roll <- strategies_linret[1:(60+(j-1)*20),]
  T_roll <- (j-1)*20+60
  T_trn_roll <- (j-1)*20+40
  ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

r_20dfs <- ret_rolling

# performance----------

# all years #
c_20dfs_1819_1p<-{ chart.CumReturns(r_20dfs, main = "Performance of months fixed start rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
}


# 2019 first 6 months #

r_20dfs_19 <- r_20dfs[which(index(r_20dfs)>='2019-01-02 CST')]
c_20dfs_19_1p<-{ chart.CumReturns(r_20dfs_19, main = "Performance of months fixed start rolling", 
                                  wealth.index = TRUE, legend.loc = "topleft", color=c(1:9))
}



#list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
     #             "MDD"=maxDrawdown(result, weights = NULL, geometric = TRUE, invert = TRUE),
     #             'Calmer'=CalmarRatio(result),
     #             'Annulized returns'=Return.annualized(result),
     #             # 'CVaR'=ETL(result, p = 0.95,method = c("modified")),
     #             "result" = result)
     
     


### figure #####################
layout(matrix(1:4,2,2))

c_5dfb_19_1p
c_20dfb_19_1p
c_5dfs_19_1p
c_20dfs_19_1p

c_5dfb_1819_1p
c_20dfb_1819_1p
c_5dfs_1819_1p
c_20dfs_1819_1p

MAR = 0.05

# Table 5 days fix block

setwd("C:/Users/m8sha/Desktop/ratio sum")

### WEEK ########################################
week_fixblock_1819 = rbind(Return.annualized(r_5dfb),
                           SortinoRatio(r_5dfb, MAR),
                           SharpeRatio(r_5dfb, Rf = 0, p = 0.95, FUN = "StdDev"),
                           CDD(r_5dfb,p=0.95),
                           ETL(r_5dfb,p=0.95),
                           ETL(r_5dfb,p=0.99),
                           CalmarRatio(r_5dfb),
                           R_CDD=(0.95*Return.annualized(r_5dfb)/CDD(r_5dfb,p=0.95)))
write.csv(week_fixblock_1819,file = "week_fixblock_1819.csv")

week_fixblock_19 = rbind(Return.annualized(r_5dfb_19),
                           SortinoRatio(r_5dfb_19, MAR),
                           SharpeRatio(r_5dfb_19, Rf = 0, p = 0.95, FUN = "StdDev"),
                           CDD(r_5dfb_19,p=0.95),
                           ETL(r_5dfb_19,p=0.95),
                           ETL(r_5dfb_19,p=0.99),
                           CalmarRatio(r_5dfb_19),
                           R_CDD=(0.95*Return.annualized(r_5dfb_19)/CDD(r_5dfb_19,p=0.95)))
write.csv(week_fixblock_19,file = "week_fixblock_19.csv")

week_fixstart_1819 = rbind(Return.annualized(r_5dfs),
                           SortinoRatio(r_5dfs, MAR),
                           SharpeRatio(r_5dfs, Rf = 0, p = 0.95, FUN = "StdDev"),
                           CDD(r_5dfs,p=0.95),
                           ETL(r_5dfs,p=0.95),
                           ETL(r_5dfs,p=0.99),
                           CalmarRatio(r_5dfs),
                           R_CDD=(0.95*Return.annualized(r_5dfs)/CDD(r_5dfs,p=0.95)))
write.csv(week_fixstart_1819,file = "week_fixstart_1819.csv")

week_fixstart_19 = rbind(Return.annualized(r_5dfs_19),
                         SortinoRatio(r_5dfs_19, MAR),
                         SharpeRatio(r_5dfs_19, Rf = 0, p = 0.95, FUN = "StdDev"),
                         CDD(r_5dfs_19,p=0.95),
                         ETL(r_5dfs_19,p=0.95),
                         ETL(r_5dfs_19,p=0.99),
                         CalmarRatio(r_5dfs_19),
                         R_CDD=(0.95*Return.annualized(r_5dfs_19)/CDD(r_5dfs_19,p=0.95)))
write.csv(week_fixstart_19,file = "week_fixstart_19.csv")

### MONTH ########################################
month_fixblock_1819 = rbind(Return.annualized(r_20dfb),
                           SortinoRatio(r_20dfb, MAR),
                           SharpeRatio(r_20dfb, Rf = 0, p = 0.95, FUN = "StdDev"),
                           CDD(r_20dfb,p=0.95),
                           ETL(r_20dfb,p=0.95),
                           ETL(r_20dfb,p=0.99),
                           CalmarRatio(r_20dfb),
                           R_CDD=(0.95*Return.annualized(r_20dfb)/CDD(r_20dfb,p=0.95)))
write.csv(month_fixblock_1819,file = "month_fixblock_1819.csv")

month_fixblock_19 = rbind(Return.annualized(r_20dfb_19),
                         SortinoRatio(r_20dfb_19, MAR),
                         SharpeRatio(r_20dfb_19, Rf = 0, p = 0.95, FUN = "StdDev"),
                         CDD(r_20dfb_19,p=0.95),
                         ETL(r_20dfb_19,p=0.95),
                         ETL(r_20dfb_19,p=0.99),
                         CalmarRatio(r_20dfb_19),
                         R_CDD=(0.95*Return.annualized(r_20dfb_19)/CDD(r_20dfb_19,p=0.95)))
write.csv(month_fixblock_19,file = "month_fixblock_19.csv")

month_fixstart_1819 = rbind(Return.annualized(r_20dfs),
                           SortinoRatio(r_20dfs, MAR),
                           SharpeRatio(r_20dfs, Rf = 0, p = 0.95, FUN = "StdDev"),
                           CDD(r_20dfs,p=0.95),
                           ETL(r_20dfs,p=0.95),
                           ETL(r_20dfs,p=0.99),
                           CalmarRatio(r_20dfs),
                           R_CDD=(0.95*Return.annualized(r_20dfs)/CDD(r_20dfs,p=0.95)))
write.csv(month_fixstart_1819,file = "month_fixstart_1819.csv")

month_fixstart_19 = rbind(Return.annualized(r_20dfs_19),
                         SortinoRatio(r_20dfs_19, MAR),
                         SharpeRatio(r_20dfs_19, Rf = 0, p = 0.95, FUN = "StdDev"),
                         CDD(r_20dfs_19,p=0.95),
                         ETL(r_20dfs_19,p=0.95),
                         ETL(r_20dfs_19,p=0.99),
                         CalmarRatio(r_20dfs_19),
                         R_CDD=(0.95*Return.annualized(r_20dfs_19)/CDD(r_20dfs_19,p=0.95)))
write.csv(month_fixstart_19,file = "month_fixstart_19.csv")
              



names(sharpe_raw) = c('r_19weeks_fixblock',
                      'r_19weeks_fixstart',
                      'r_19months_fixblock',
                      'r_19months_fixblock',
                      'r_1819weeks_fixblock',
                      'r_1819weeks_norolling',
                      'r_1819months_fixblock',
                      'r_1819months_fixstart')

SM = rbind(sharpe_raw,MDD_raw)    
