library(xts) 
library(quantmod) 
library(PerformanceAnalytics)

### Set portfolio Funtions ###

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
w_Lquintile <- cbind(w_Lquintile1, w_Lquintile2, w_Lquintile3) 
rownames(w_Lquintile) <- colnames(strategies_log_roll)
colnames(w_Lquintile) <- c("Lquintile1", "Lquintile2", "Lquintile3") 
#w_Lquintile

### Return-risk tradeoff for all portfolios ###
# put together all portfolios
w_meanvar <- cbind(w_unif, w_GMVP, w_Markowitz, w_maxSR) 
rownames(w_meanvar) <- colnames(strategies_log_roll)
colnames(w_meanvar) <- c("uniform", "GMVP", "Markowitz", "maxSR") 
w_all <- cbind(w_meanvar, w_Lquintile)
round(w_all, digits = 2)

### Then compute the performance ###

# compute returns of all portfolios
ret_all <- xts(strategies_lin_roll %*% w_all, index(strategies_lin_roll))  
ret_all_trn <- ret_all[1:T_trn_roll, ]
ret_all_tst <- ret_all[-c(1:T_trn_roll), ]
return(ret_all_tst)
}

#-----------------------------------------------------------------------------------

#### I. strategies ####

### Get data ###
setwd('/Users/jingyuegao/Desktop')

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
head(strategies_linret)
#head(strategies_logret)
#tail(strategies_linret)
#tail(strategies_logret)
N <- ncol(strategies_logret)

#### Rolling ####

### 1. For weeks ###

## 1.1 15-5 days test ##

# 1.1.1 fixed block #
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

# performance----------

# all years #
c_5dfb_1819_1p <-{ chart.CumReturns(r_5dfb, main = "Performance of 5 day fixed blocks rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
  }

c_5dfb_1819_3p<- charts.PerformanceSummary(r_5dfb, main = "Performance of 5 day fixed blocks rolling",
                          wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_5dfb_19 <- r_5dfb[which(index(r_5dfb)>='2019-01-02 CST')]
c_5dfb_19_1p <- { chart.CumReturns(r_5dfb_19, main = "Performance of 5 day fixed blocks rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
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
c_5dfs_1819_1p <-{ chart.CumReturns(r_5dfs, main = "Performance of 5 day fixed start rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
  }

c_5dfs_1819_3p<- charts.PerformanceSummary(r_5dfs, main = "Performance of 5 day fixed start rolling",
                          wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_5dfs_19 <- r_5dfs[which(index(r_5dfs)>='2019-01-02 CST')]
c_5dfs_19_1p <- { chart.CumReturns(r_5dfs_19, main = "Performance of 5 day fixed start rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
  }

## 1.2 one week test ##

dayname <- substr(index(strategies_logret),1,10)
weekname <- dayname[1]
for (m in 1:(length(dayname)-1)){
  if (index(strategies_logret)[m+1]-index(strategies_logret)[m]!=1){
    weekname <- rbind(weekname, substr(dayname[m+1],1,10))
  }
}   # get the first day of each week 

# 1.2.1 fixed block #

# For first block #

T_roll <- length(dayname[which(dayname < weekname[5])])
T_trn_roll <- length(dayname[which(dayname < weekname[4])])

strategies_log_roll_1st <- strategies_logret[(1:T_roll),]
strategies_lin_roll_1st <- strategies_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(weekname)-4

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < weekname[j+4])])-length(dayname[which(dayname < weekname[j])])
  T_trn_roll <- length(dayname[which(dayname < weekname[j+3])])-length(dayname[which(dayname < weekname[j])])
  
  strategies_log_roll <- strategies_logret[(length(dayname[which(dayname < weekname[j])])+1):length(dayname[which(dayname < weekname[j+4])]),]
  strategies_lin_roll <- strategies_linret[(length(dayname[which(dayname < weekname[j])])+1):length(dayname[which(dayname < weekname[j+4])]),]
  ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)-length(dayname[which(dayname < weekname[length(weekname)-3])])
T_trn_roll <- length(dayname[which(dayname < weekname[length(weekname)])])-length(dayname[which(dayname < weekname[length(weekname)-3])])

strategies_log_roll <- strategies_logret[(length(dayname[which(dayname < weekname[length(weekname)-3])])+1):length(dayname),]
strategies_lin_roll <- strategies_linret[(length(dayname[which(dayname < weekname[length(weekname)-3])])+1):length(dayname),]
ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

r_1weekfb <- ret_rolling

# performance----------

# all years #
c_1weekfb_1819_1p<-{ chart.CumReturns(r_1weekfb, main = "Performance of 1 week fixed blocks rolling", 
                                     wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

c_1weekfb_1819_3p<-charts.PerformanceSummary(r_1weekfb, main = "Performance of 1 week fixed blocks rolling",
                                            wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_1weekfb_19 <- r_1weekfb[which(index(r_1weekfb)>='2019-01-02 CST')]
c_1weekfb_19_1p<-{ chart.CumReturns(r_1weekfb_19, main = "Performance of 1 week fixed blocks rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

# 1.2.2 fixed start #

# For first block #

T_roll <- length(dayname[which(dayname < weekname[5])])
T_trn_roll <- length(dayname[which(dayname < weekname[4])])

strategies_log_roll_1st <- strategies_logret[(1:T_roll),]
strategies_lin_roll_1st <- strategies_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(weekname)-4

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < weekname[j+4])])
  T_trn_roll <- length(dayname[which(dayname < weekname[j+3])])
  
  strategies_log_roll <- strategies_logret[1:length(dayname[which(dayname < weekname[j+4])]),]
  strategies_lin_roll <- strategies_linret[1:length(dayname[which(dayname < weekname[j+4])]),]
  ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)
T_trn_roll <- length(dayname[which(dayname < weekname[length(weekname)])])

strategies_log_roll <- strategies_logret[1:length(dayname),]
strategies_lin_roll <- strategies_linret[1:length(dayname),]
ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

r_1weekfs <- ret_rolling

# performance----------

# all years #
c_1weekfs_1819_1p<-{ chart.CumReturns(r_1weekfs, main = "Performance of 1 week fixed start rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

c_1weekfs_1819_3p<-charts.PerformanceSummary(r_1weekfs, main = "Performance of 20 day fixed start rolling",
                                           wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_1weekfs_19 <- r_1weekfs[which(index(r_1weekfs)>='2019-01-02 CST')]
c_1weekfs_19_1p<-{ chart.CumReturns(r_1weekfs_19, main = "Performance of 1 week fixed start rolling", 
                                  wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
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
c_20dfb_1819_1p <-{ chart.CumReturns(r_20dfb, main = "Performance of 20 day fixed blocks rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
 }

c_20dfb_1819_3p<-charts.PerformanceSummary(r_20dfb, main = "Performance of different portfolios",
                          wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_20dfb_19 <- r_20dfb[which(index(r_20dfb)>='2019-01-02 CST')]
c_20dfb_19_1p<-{ chart.CumReturns(r_20dfb_19, main = "Performance of 20 day fixed blocks rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
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
c_20dfs_1819_1p<-{ chart.CumReturns(r_20dfs, main = "Performance of 20 day fixed start rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

c_20dfs_1819_3p<-charts.PerformanceSummary(r_20dfs, main = "Performance of 20 day fixed start rolling",
                          wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_20dfs_19 <- r_20dfs[which(index(r_20dfs)>='2019-01-02 CST')]
c_20dfs_19_1p<-{ chart.CumReturns(r_20dfs_19, main = "Performance of 20 day fixed start rolling", 
                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

## 2.2 one month test ##

dayname <- substr(index(strategies_logret),1,10)
monthname <- dayname[1]
for (m in 1:(length(dayname)-1)){
  if (substr(dayname[m],6,7)!=substr(dayname[m+1],6,7)){
    monthname <- rbind(monthname, substr(dayname[m+1],1,10))
  }
}   # get the first day of each month 

# 2.2.1 fixed block #

# For first block #

T_roll <- length(dayname[which(dayname < monthname[4])])
T_trn_roll <- length(dayname[which(dayname < monthname[3])])

strategies_log_roll_1st <- strategies_logret[(1:T_roll),]
strategies_lin_roll_1st <- strategies_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(monthname)-3

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < monthname[j+3])])-length(dayname[which(dayname < monthname[j])])
  T_trn_roll <- length(dayname[which(dayname < monthname[j+2])])-length(dayname[which(dayname < monthname[j])])
  
  strategies_log_roll <- strategies_logret[(length(dayname[which(dayname < monthname[j])])+1):length(dayname[which(dayname < monthname[j+3])]),]
  strategies_lin_roll <- strategies_linret[(length(dayname[which(dayname < monthname[j])])+1):length(dayname[which(dayname < monthname[j+3])]),]
  ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)-length(dayname[which(dayname < monthname[length(monthname)-2])])
T_trn_roll <- length(dayname[which(dayname < monthname[length(monthname)])])-length(dayname[which(dayname < monthname[length(monthname)-2])])

strategies_log_roll <- strategies_logret[(length(dayname[which(dayname < monthname[length(monthname)-2])])+1):length(dayname),]
strategies_lin_roll <- strategies_linret[(length(dayname[which(dayname < monthname[length(monthname)-2])])+1):length(dayname),]
ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

r_1monfb <- ret_rolling

# performance----------

# all years #
c_1monfb_1819_1p<-{ chart.CumReturns(r_1monfb, main = "Performance of 1 month fixed blocks rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

c_1monfb_1819_3p<-charts.PerformanceSummary(r_1monfb, main = "Performance of 1 month fixed blocks rolling",
                                           wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_1monfb_19 <- r_1monfb[which(index(r_1monfb)>='2019-01-02 CST')]
c_1monfb_19_1p<-{ chart.CumReturns(r_1monfb_19, main = "Performance of 1 month fixed blocks rolling", 
                                  wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}


# 2.2.2 fixed start #

# For first block #

T_roll <- length(dayname[which(dayname < monthname[4])])
T_trn_roll <- length(dayname[which(dayname < monthname[3])])

strategies_log_roll_1st <- strategies_logret[(1:T_roll),]
strategies_lin_roll_1st <- strategies_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(strategies_log_roll_1st,strategies_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(monthname)-3

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < monthname[j+3])])
  T_trn_roll <- length(dayname[which(dayname < monthname[j+2])])
  
  strategies_log_roll <- strategies_logret[1:length(dayname[which(dayname < monthname[j+3])]),]
  strategies_lin_roll <- strategies_linret[1:length(dayname[which(dayname < monthname[j+3])]),]
  ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)
T_trn_roll <- length(dayname[which(dayname < monthname[length(monthname)])])

strategies_log_roll <- strategies_logret[1:length(dayname),]
strategies_lin_roll <- strategies_linret[1:length(dayname),]
ret_add <- CalculateTestReturns(strategies_log_roll,strategies_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

r_1monfs <- ret_rolling

# performance----------

# all years #
c_1monfs_1819_1p<-{ chart.CumReturns(r_1monfs, main = "Performance of 1 month fixed start rolling", 
                                     wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

c_1monfs_1819_3p<-charts.PerformanceSummary(r_1monfs, main = "Performance of 1 month fixed start rolling",
                                            wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

r_1monfs_19 <- r_1monfs[which(index(r_1monfs)>='2019-01-02 CST')]
c_1monfs_19_1p<-{ chart.CumReturns(r_1monfs_19, main = "Performance of 1 month fixed start rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}


#-----------------------------------------------------------------------------

#### II. categories ####


### Get Data ###
categories_linret <- read.csv('Contract.Dailylinreturn.txt',header = TRUE,sep = '')
categories_logret <- read.csv('Contract.Dailylogreturn.txt',header = TRUE,sep = '')
rownames(categories_linret) <- substr(categories_linret[,1],1,10)
rownames(categories_logret) <- substr(categories_logret[,1],1,10)

categories_linret <- categories_linret[,-1]
categories_logret <- categories_logret[,-1]
categories_linret <- as.xts(categories_linret)
categories_logret <- as.xts(categories_logret)
str(categories_linret)
str(categories_logret)
head(categories_linret)
#head(categories_logret)
#tail(categories_linret)
#tail(categories_logret)
N <- ncol(categories_logret)






#### Rolling ####

### 1. For weeks ###

## 1.1 15-5 days test ##

# 1.1.1 fixed block #
T_roll <- 20
T_trn_roll <- 15

categories_log_roll_1st <- categories_logret[(1+(1-1)*5):(20+(1-1)*5),]
categories_lin_roll_1st <- categories_linret[(1+(1-1)*5):(20+(1-1)*5),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,T_roll,T_trn_roll) 
last_test <- (floor(nrow(categories_logret)/5)-3)

categories_log_roll <- categories_logret[(1+(4-1)*5):(20+(4-1)*5),]
categories_lin_roll <- categories_linret[(1+(4-1)*5):(20+(4-1)*5),]
ret_rolling <- rbind(ret_rolling,CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll))



for (j in 2:last_test){
  categories_log_roll <- categories_logret[(1+(j-1)*5):(20+(j-1)*5),]
  categories_lin_roll <- categories_linret[(1+(j-1)*5):(20+(j-1)*5),]
  ret_rolling <- rbind(ret_rolling,CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll))
}
cr_5dfb <- ret_rolling

# performance----------

# all years #
cc_5dfb_1819_1p <-{ chart.CumReturns(cr_5dfb, main = "Performance of 5 day fixed blocks rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_5dfb_1819_3p<- charts.PerformanceSummary(cr_5dfb, main = "Performance of 5 day fixed blocks rolling",
                                           wealth.index = TRUE, colorset = rich8equal)


# 2019 first 6 months #

cr_5dfb_19 <- cr_5dfb[which(index(cr_5dfb)>='2019-01-02 CST')]
cc_5dfb_19_1p <- { chart.CumReturns(cr_5dfb_19, main = "Performance of 5 day fixed blocks rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

# 1.1.2 fixed start #

categories_log_roll_1st <- categories_logret[(1+(1-1)*5):(20+(1-1)*5),]
categories_lin_roll_1st <- categories_linret[(1+(1-1)*5):(20+(1-1)*5),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,20,15)
last_test <- (floor(nrow(categories_logret)/5)-3)

for (j in 2:last_test){
  categories_log_roll <- categories_logret[1:(20+(j-1)*5),]
  categories_lin_roll <- categories_linret[1:(20+(j-1)*5),]
  T_roll <- (j-1)*5+20
  T_trn_roll <- (j-1)*5+15
  ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}
cr_5dfs <- ret_rolling

# performance--------

# all years #
cc_5dfs_1819_1p <-{ chart.CumReturns(cr_5dfs, main = "Performance of 5 day fixed start rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_5dfs_1819_3p<- charts.PerformanceSummary(cr_5dfs, main = "Performance of 5 day fixed start rolling",
                                           wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

cr_5dfs_19 <- cr_5dfs[which(index(cr_5dfs)>='2019-01-02 CST')]
cc_5dfs_19_1p <- { chart.CumReturns(cr_5dfs_19, main = "Performance of 5 day fixed start rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

## 1.2 one week test ##

dayname <- substr(index(categories_logret),1,10)
weekname <- dayname[1]
for (m in 1:(length(dayname)-1)){
  if (index(categories_logret)[m+1]-index(categories_logret)[m]!=1){
    weekname <- rbind(weekname, substr(dayname[m+1],1,10))
  }
}   # get the first day of each week 

# 1.2.1 fixed block #

# For first block #

T_roll <- length(dayname[which(dayname < weekname[5])])
T_trn_roll <- length(dayname[which(dayname < weekname[4])])

categories_log_roll_1st <- categories_logret[(1:T_roll),]
categories_lin_roll_1st <- categories_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(weekname)-4

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < weekname[j+4])])-length(dayname[which(dayname < weekname[j])])
  T_trn_roll <- length(dayname[which(dayname < weekname[j+3])])-length(dayname[which(dayname < weekname[j])])
  
  categories_log_roll <- categories_logret[(length(dayname[which(dayname < weekname[j])])+1):length(dayname[which(dayname < weekname[j+4])]),]
  categories_lin_roll <- categories_linret[(length(dayname[which(dayname < weekname[j])])+1):length(dayname[which(dayname < weekname[j+4])]),]
  ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)-length(dayname[which(dayname < weekname[length(weekname)-3])])
T_trn_roll <- length(dayname[which(dayname < weekname[length(weekname)])])-length(dayname[which(dayname < weekname[length(weekname)-3])])

categories_log_roll <- categories_logret[(length(dayname[which(dayname < weekname[length(weekname)-3])])+1):length(dayname),]
categories_lin_roll <- categories_linret[(length(dayname[which(dayname < weekname[length(weekname)-3])])+1):length(dayname),]
ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

cr_1weekfb <- ret_rolling

# performance----------

# all years #
cc_1weekfb_1819_1p<-{ chart.CumReturns(cr_1weekfb, main = "Performance of 1 week fixed blocks rolling", 
                                      wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_1weekfb_1819_3p<-charts.PerformanceSummary(cr_1weekfb, main = "Performance of 1 week fixed blocks rolling",
                                             wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

cr_1weekfb_19 <- cr_1weekfb[which(index(cr_1weekfb)>='2019-01-02 CST')]
cc_1weekfb_19_1p<-{ chart.CumReturns(cr_1weekfb_19, main = "Performance of 1 week fixed blocks rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

# 1.2.2 fixed start #

# For first block #

T_roll <- length(dayname[which(dayname < weekname[5])])
T_trn_roll <- length(dayname[which(dayname < weekname[4])])

categories_log_roll_1st <- categories_logret[(1:T_roll),]
categories_lin_roll_1st <- categories_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(weekname)-4

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < weekname[j+4])])
  T_trn_roll <- length(dayname[which(dayname < weekname[j+3])])
  
  categories_log_roll <- categories_logret[1:length(dayname[which(dayname < weekname[j+4])]),]
  categories_lin_roll <- categories_linret[1:length(dayname[which(dayname < weekname[j+4])]),]
  ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)
T_trn_roll <- length(dayname[which(dayname < weekname[length(weekname)])])

categories_log_roll <- categories_logret[1:length(dayname),]
categories_lin_roll <- categories_linret[1:length(dayname),]
ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

cr_1weekfs <- ret_rolling

# performance----------

# all years #
cc_1weekfs_1819_1p<-{ chart.CumReturns(cr_1weekfs, main = "Performance of 1 week fixed start rolling", 
                                      wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_1weekfs_1819_3p<-charts.PerformanceSummary(cr_1weekfs, main = "Performance of 20 day fixed start rolling",
                                             wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

cr_1weekfs_19 <- cr_1weekfs[which(index(cr_1weekfs)>='2019-01-02 CST')]
cc_1weekfs_19_1p<-{ chart.CumReturns(cr_1weekfs_19, main = "Performance of 1 week fixed start rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

### 2. for months ###

## 2.1 40-20 days test 

# 2.1.1 fixed block #
T_roll <- 60
T_trn_roll <- 40

# For first 20 d #
categories_log_roll_1st <- categories_logret[(1+(1-1)*20):(60+(1-1)*20),]
categories_lin_roll_1st <- categories_linret[(1+(1-1)*20):(60+(1-1)*20),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,60,40)
last_test <- (floor(nrow(categories_logret)/20)-3)

for (j in 2:last_test){
  categories_log_roll <- categories_logret[(1+(j-1)*20):(60+(j-1)*20),]
  categories_lin_roll <- categories_linret[(1+(j-1)*20):(60+(j-1)*20),]
  ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

cr_20dfb <- ret_rolling

# performance---------

# all years #
cc_20dfb_1819_1p <-{ chart.CumReturns(cr_20dfb, main = "Performance of 20 day fixed blocks rolling", 
                                     wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_20dfb_1819_3p<-charts.PerformanceSummary(cr_20dfb, main = "Performance of different portfolios",
                                           wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

cr_20dfb_19 <- cr_20dfb[which(index(cr_20dfb)>='2019-01-02 CST')]
cc_20dfb_19_1p<-{ chart.CumReturns(cr_20dfb_19, main = "Performance of 20 day fixed blocks rolling", 
                                  wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}


# 2.1.2  fixed start #

categories_log_roll_1st <- categories_logret[(1+(1-1)*20):(60+(1-1)*20),]
categories_lin_roll_1st <- categories_linret[(1+(1-1)*20):(60+(1-1)*20),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,60,40)
last_test <- (floor(nrow(categories_logret)/20)-3)

for (j in 2:last_test){
  categories_log_roll <- categories_logret[1:(60+(j-1)*20),]
  categories_lin_roll <- categories_linret[1:(60+(j-1)*20),]
  T_roll <- (j-1)*20+60
  T_trn_roll <- (j-1)*20+40
  ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

cr_20dfs <- ret_rolling

# performance----------

# all years #
cc_20dfs_1819_1p<-{ chart.CumReturns(cr_20dfs, main = "Performance of 20 day fixed start rolling", 
                                    wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_20dfs_1819_3p<-charts.PerformanceSummary(cr_20dfs, main = "Performance of 20 day fixed start rolling",
                                           wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

cr_20dfs_19 <- cr_20dfs[which(index(cr_20dfs)>='2019-01-02 CST')]
cc_20dfs_19_1p<-{ chart.CumReturns(cr_20dfs_19, main = "Performance of 20 day fixed start rolling", 
                                  wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

## 2.2 one month test ##

dayname <- substr(index(categories_logret),1,10)
monthname <- dayname[1]
for (m in 1:(length(dayname)-1)){
  if (substr(dayname[m],6,7)!=substr(dayname[m+1],6,7)){
    monthname <- rbind(monthname, substr(dayname[m+1],1,10))
  }
}   # get the first day of each month 

# 2.2.1 fixed block #

# For first block #

T_roll <- length(dayname[which(dayname < monthname[4])])
T_trn_roll <- length(dayname[which(dayname < monthname[3])])

categories_log_roll_1st <- categories_logret[(1:T_roll),]
categories_lin_roll_1st <- categories_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(monthname)-3

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < monthname[j+3])])-length(dayname[which(dayname < monthname[j])])
  T_trn_roll <- length(dayname[which(dayname < monthname[j+2])])-length(dayname[which(dayname < monthname[j])])
  
  categories_log_roll <- categories_logret[(length(dayname[which(dayname < monthname[j])])+1):length(dayname[which(dayname < monthname[j+3])]),]
  categories_lin_roll <- categories_linret[(length(dayname[which(dayname < monthname[j])])+1):length(dayname[which(dayname < monthname[j+3])]),]
  ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)-length(dayname[which(dayname < monthname[length(monthname)-2])])
T_trn_roll <- length(dayname[which(dayname < monthname[length(monthname)])])-length(dayname[which(dayname < monthname[length(monthname)-2])])

categories_log_roll <- categories_logret[(length(dayname[which(dayname < monthname[length(monthname)-2])])+1):length(dayname),]
categories_lin_roll <- categories_linret[(length(dayname[which(dayname < monthname[length(monthname)-2])])+1):length(dayname),]
ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

cr_1monfb <- ret_rolling

# performance----------

# all years #
cc_1monfb_1819_1p<-{ chart.CumReturns(cr_1monfb, main = "Performance of 1 month fixed blocks rolling", 
                                     wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_1monfb_1819_3p<-charts.PerformanceSummary(cr_1monfb, main = "Performance of 1 month fixed blocks rolling",
                                            wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

cr_1monfb_19 <- cr_1monfb[which(index(cr_1monfb)>='2019-01-02 CST')]
cc_1monfb_19_1p<-{ chart.CumReturns(cr_1monfb_19, main = "Performance of 1 month fixed blocks rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}


# 2.2.2 fixed start #

# For first block #

T_roll <- length(dayname[which(dayname < monthname[4])])
T_trn_roll <- length(dayname[which(dayname < monthname[3])])

categories_log_roll_1st <- categories_logret[(1:T_roll),]
categories_lin_roll_1st <- categories_linret[(1:T_roll),]
ret_rolling <- CalculateTestReturns(categories_log_roll_1st,categories_lin_roll_1st,N,T_roll,T_trn_roll)
last_test <- nrow(monthname)-3

for (j in 2: last_test){
  T_roll <- length(dayname[which(dayname < monthname[j+3])])
  T_trn_roll <- length(dayname[which(dayname < monthname[j+2])])
  
  categories_log_roll <- categories_logret[1:length(dayname[which(dayname < monthname[j+3])]),]
  categories_lin_roll <- categories_linret[1:length(dayname[which(dayname < monthname[j+3])]),]
  ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
  ret_rolling <- rbind(ret_rolling,ret_add)
}

# For last block #

T_roll <- length(dayname)
T_trn_roll <- length(dayname[which(dayname < monthname[length(monthname)])])

categories_log_roll <- categories_logret[1:length(dayname),]
categories_lin_roll <- categories_linret[1:length(dayname),]
ret_add <- CalculateTestReturns(categories_log_roll,categories_lin_roll,N,T_roll,T_trn_roll)
ret_rolling <- rbind(ret_rolling,ret_add)

cr_1monfs <- ret_rolling

# performance----------

# all years #
cc_1monfs_1819_1p<-{ chart.CumReturns(cr_1monfs, main = "Performance of 1 month fixed start rolling", 
                                     wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

cc_1monfs_1819_3p<-charts.PerformanceSummary(cr_1monfs, main = "Performance of 1 month fixed start rolling",
                                            wealth.index = TRUE, colorset = rich8equal)

# 2019 first 6 months #

cr_1monfs_19 <- cr_1monfs[which(index(cr_1monfs)>='2019-01-02 CST')]
cc_1monfs_19_1p<-{ chart.CumReturns(cr_1monfs_19, main = "Performance of 1 month fixed start rolling", 
                                   wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)
}

### show graphs ###

layout(matrix(1:8,4,2))

c_5dfb_19_1p
c_1weekfb_19_1p
c_20dfb_19_1p
c_1monfb_19_1p
c_5dfs_19_1p
c_1weekfs_19_1p
c_20dfs_19_1p
c_1monfs_19_1p

cc_5dfb_19_1p
cc_1weekfb_19_1p
cc_20dfb_19_1p
cc_1monfb_19_1p
cc_5dfs_19_1p
cc_1weekfs_19_1p
cc_20dfs_19_1p
cc_1monfs_19_1p

layout(matrix(1,1,1))

