library('clusterSim')
library('ade4')
library('PerformanceAnalytics')
library(CVXR)



####### 1 JSD ####################################
## dist JSD alculates the Jensen Shannon Divergence
# input: inMatrix: a numeric matrix
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- dim(inMatrix)[2]
  matrixRowSize <- dim(inMatrix)[1]
  #matrixColSize <- length(colnames(inMatrix))
  #matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}


####### 2.1 performance ############################
# performance returns a matrix containing 31 ratio based on input signal and return,
# input:
# 1. allSig: a matrix of signal per minute for all strategies
# 2. rt: a xts of log return perminute for all strategies
# 3. rt.day: a xts of log return per day for all strategies

performance<-function(allSig,rt,rt.day){
  perf.m= matrix(nrow = 0,ncol = ncol(rt.day)) 
  colnames(perf.m) = colnames(rt.day)
  
  # performance
  # function for total number& gain.loss ratio & expectation
  CalculateNumbers<-function(allSig,rt,rt.day){
    ####连续亏损 (day)
    tmp = matrix(0,nrow(rt.day),ncol(rt.day))
    for(j in 1:ncol(rt.day)){
      for(i in 1:nrow(rt.day)){
        if(rt.day[i,j]<0){
          tmp[i,j] = -1
        }
      }
    }
    
    for(j in 1:ncol(rt.day)){
      for(i in 2:nrow(rt.day)){
        if(tmp[i,j]<0 && tmp[i-1,j]<0){
          tmp[i,j] = tmp[i-1,j]-1
        }
      }
    }
    
    contloss = abs(apply(tmp,2,min))
    #
    #   ## 交易次数 胜率 期望
    totalnumber = numeric(length = ncol(allSig))
    rtByTrade = matrix(0,500,ncol(allSig))
    
    for(j in 1:ncol(allSig)){
      for (i in 2:nrow(allSig)) {
        print(i)
        if((allSig[i,j]==1 || allSig[i,j] == -1) && allSig[i,j] != allSig[i-1,j]){
          totalnumber[j] = totalnumber[j]+1
          rtByTrade[totalnumber[j],j] = i
        }
      }
    }
    
    
    for(j in 1:ncol(rtByTrade)){
      for(i in 1:nrow(rtByTrade)){
        if( rtByTrade[i,j] != 0){
          if( rtByTrade[i+1,j] == 0){ rtByTrade[i,j] = sum(rt[rtByTrade[i,j]:nrow(rt),j]) }
          else{
            rtByTrade[i,j] = sum(rt[rtByTrade[i,j]:(rtByTrade[i+1,j]-1),j])
          }
        }
      }
    }
    
    gain.ratio = numeric(length = ncol(allSig))
    expe = numeric(length = ncol(allSig))
    gain.avg = numeric(length = ncol(allSig))
    loss.avg = numeric(length = ncol(allSig))
    for(j in 1:ncol(rtByTrade)){
      gainNum = sum(rtByTrade[,j]>0)
      lossNum = sum(rtByTrade[,j]<0)
      totalGain = sum(rtByTrade[rtByTrade[,j]>0,j])
      totalLoss = sum(rtByTrade[rtByTrade[,j]<0,j])
      gain.ratio[j] = gainNum/totalnumber[j]
      expe[j] = (gainNum*totalGain + lossNum*totalLoss)/totalnumber[j]
      gain.avg[j] = totalGain/gainNum
      loss.avg[j] = totalLoss/lossNum
    }
    return(list('contlossdays'=contloss,
                'totalnumber'=totalnumber,
                'gainratio'=gain.ratio,
                'expectation'=expe,
                'gainavg'=gain.avg,
                'lossavg'=loss.avg))
  }
  numrelated = CalculateNumbers(allSig,rt,rt.day)
  
  # information ratio
  benchmark = matrix(0,nrow = dim(rt.day)[1], ncol = dim(rt.day)[2])
  colnames(benchmark) =colnames(rt.day)
  rownames(benchmark) = rownames(as.matrix(rt.day))
  
  # for (i in 1:ncol(rt.day)) {
  #   benchmark[,i] = rep(0.035,nrow(rt.day))
  # }
  
  benchmark<-as.xts(benchmark)
  
  inf_ratio <- InformationRatio(rt.day, benchmark, scale = 252)
  inf_ratio = c(inf_ratio[2,1],inf_ratio[1,2:ncol(rt.day)])
  
  MAR=0
  
  perf.m=rbind(perf.m,
               colMeans(rt.day),
               diag(cov(rt.day)),
               table.AnnualizedReturns(rt.day),
               BernardoLedoitRatio(rt.day),
               BurkeRatio(rt.day, Rf = 0, modified = FALSE),
               CalmarRatio(rt.day, scale = NA),
               #CDD(rt.day, weights = NULL, geometric = TRUE, invert = TRUE, p = 0.95),
               DownsideDeviation(rt.day, MAR = 0, method = c("full", "subset"),potential = FALSE),
               DRatio(rt.day),
               DrawdownDeviation(rt.day),
               inf_ratio,
               KellyRatio(rt.day, Rf = 0, method = "half"),
               MartinRatio(rt.day, Rf = 0),
               maxDrawdown(rt.day, weights = NULL, geometric = TRUE, invert = TRUE),
               MeanAbsoluteDeviation(rt.day),
               Omega(rt.day, L = 0, method = "simple", output = "point", Rf = 0),
               #PainRatio(rt.day, Rf = 0),
               ProspectRatio(rt.day,MAR),
               SharpeRatio(rt.day, Rf = 0, p = 0.95, FUN = "StdDev"),
               SkewnessKurtosisRatio(rt.day),
               #Skewness(rt.day),
               #Kurtosis(rt.day),
               SmoothingIndex(rt.day, neg.thetas = FALSE, MAorder = 2, verbose = FALSE),
               SortinoRatio(rt.day, MAR),
               UlcerIndex(rt.day),
               UpsideFrequency(rt.day, MAR = 0),
               VaR(rt.day, p = 0.95,method = "historical", portfolio_method ="single",invert = TRUE),
               numrelated$contloss,
               numrelated$totalnumber,numrelated$gainratio,numrelated$expectation,
               numrelated$gainavg,numrelated$lossavg)
  
  rownames(perf.m)[1:2] = c('mu','cov')
  rownames(perf.m)[12] = c('Information ratio')
  rownames(perf.m)[26:31] = c('continuous loss days',"total numbers", "gain ratio",
                              "expectation",'avg gain','avg loss')
  
  return(perf.m)
}


####### 2.2 PCA of performance #####################
# PCAperf determines the most optimal number of factors and
# plots the eigen values for a principal components and the factor solution 
## input
# 1. inperf is a matrix from function performance
# 2. fact = c('ratio','strategy')
# 3. rotate = c("none", "varimax", "quartimax", "promax", "oblimin", "simplimax","cluster" )

PCAperf<-function(inPerf,fact,rotate){
  if (fact == 'ratio'){
    sigma<-cor(as.matrix(t(inPerf)))
    n = nrow(inPerf)}
  else{
    sigma<-cor(as.matrix(inPerf))
    n = ncol(inPerf)}
  sigma<-cor(as.matrix(t(perf)))
  #e<-eigen(sigma)
  # eigenvalue
  # e2$values
  # eigenvector
  # e2$vectors
  
  # 确定主成分个数
  #画碎石图
  a = fa.parallel(as.matrix(sigma),n.obs=n,fa='pc')
  abline(h=1)
  
  # 提取主成分
  fit<-principal(sigma,
                 nfactors = (a$ncomp+1),
                 rotate='varimax',  # max variance
                 scores=T)
  
  # 绘制主成分的载荷矩阵，查看各个主成分的综合构成变量
  diag = fa.diagram(fit,digits=2)
  
  return(list('fit'=fit,
              'diag'=diag))
}




####### 3.1 top10 => top5 & rolling ################################
# rolling window is either weekly or monthly
## input:
# 1&2. X_log&X_lin are matricies of log or linear return
# 3. on = c('months','weeks')
# 4. start = c(0,1),0 stands for block, 1 stands for fix start
rollingTop <- function(X_log,X_lin,on,start){
  if (on == 'week'){T_trn <- 15}
  else{T_trn <-60}
  t <- nrow(X_log)
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):t, ]
  X_lin_trn <- X_lin[1:T_trn, ]
  X_lin_tst <- X_lin[(T_trn+1):t, ]
  # Rank_rolling <- X_log
  # Rank_rolling[] <- NA
  rebal_indices <- T_trn + endpoints(X_log_tst, on = on)
  #rebal_indices2 <- T_trn + endpoints(X_log_tst, on = "months")
  # index(X_log)[rebal_indices]
  lookback <- 3*5  
  result =c()
  for (i in 1:(length(rebal_indices)-1)){
    #Training
    if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
    else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
    #get perfornamce: return and MDD
    perf_rt_mdd = rbind(table.AnnualizedReturns(X_)[1,], maxDrawdown(X_, weights = NULL, geometric = TRUE, invert = TRUE))
    # get top 10 return strategies
    rt_top10 = sort(as.numeric(perf_rt_mdd[1,]), decreasing= TRUE, index.return = TRUE)
    rt_top10 = perf_rt_mdd[,c(rt_top10$ix[1:10])]
    # get top 5 MDD
    mdd_top5 = colnames(rt_top10)[sort(as.numeric(rt_top10[2,]), decreasing= FALSE, index.return = TRUE)$ix[1:5]]
    # get test strategies with responding days
    tst <- X_lin[(rebal_indices[i]+1):min(rebal_indices[i+1],nrow(X_log)),mdd_top5] 
    # calculate mean returns as result for each days
    result <- c(result,apply(matrix(tst,ncol=5), 1, mean))
  }
  
  names(result) = rownames(X_log)[(T_trn+1):t]
  
  return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
              "MDD"=maxDrawdown(result, weights = NULL, geometric = TRUE, invert = TRUE),
              "result" = result))
  
}


####### 3.2 Portfolio optimization #################################
######## 3.2.1.1 Set portfolio Funtions ##############################
# Uniform, GMVP, Markowitz, MaxSR, DR, CVaR, MAX-DD,Ave-DD,CDaR portfolio 

# GMVP #
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




######## 3.2.1.2 calculate weights for all kinds of portfolio ##########
# uniform portfolio, GMVP portfolio, Markowitz portfolio, Long-Short quintile portfolio
# Alternative Risk Measures based portfolio

CalculateWeight1 <- function(strategies_log_roll){
  
  N = ncol(strategies_log_roll)
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

## Risk-Parity portfolio
CalculateWeight2 <-function(X_log_trn){
  mu <- colMeans(X_log_trn)
  Sigma <- cov(as.matrix(X_log_trn))
  
  w_Markowitz <- portolioMarkowitz(mu, Sigma)
  w_GMVP <- portolioGMVP(Sigma)
  w_all <- cbind(w_GMVP, w_Markowitz)
  rownames(w_all) <- colnames(X_log_trn)
  colnames(w_all) <- c("GMVP", "Markowitz")
  
  ## risk-parity-diag
  w_diag <- 1/sqrt(diag(Sigma))
  w_diag <- w_diag/sum(w_diag)
  w_all <- cbind(w_all, "risk-parity-diag" = w_diag)
  
  ## risk-parity-convex
  N <- nrow(Sigma)
  x0 <- rep(1/N, N) # initial point
  fn_convex <- function(x, Sigma) {
    N <- nrow(Sigma)
    0.5 * t(x) %*% Sigma %*% x - (1/N)*sum(log(x))
  }
  # general solver
  result <- optim(par = x0, fn = fn_convex, Sigma = Sigma,
                  method = "BFGS")
  x_convex <- result$par
  w_convex <- x_convex/sum(x_convex)
  
  w_all <- cbind(w_all, "risk-parity-convex" = w_convex)
  
  ## risk-parity-gen-solver
  w0 <- rep(1/N, N) # initial point
  fn <- function(w, Sigma) {
    N <- length(w)
    risks <- w * (Sigma %*% w)
    g <- rep(risks, times = N) - rep(risks, each = N)
    return(sum(g^2))
  }
  # general solver
  result <- optim(par = w0, fn = fn, Sigma = Sigma,
                  method = "BFGS")
  w_gen_solver <- result$par
  w_gen_solver <- w_gen_solver/sum(w_gen_solver)
  w_all <- cbind(w_all, "risk-parity-gen-solver" = w_gen_solver)
  
  ## risk-parity-package
  # devtools::install_github("dppalomar/riskParityPortfolio")
  #library(riskParityPortfolio)
  rpp <- riskParityPortfolio(Sigma)
  w_all <- cbind(w_all, "risk-parity-package" = rpp$w)
  
  ### SCA
  max_iter <- 40
  tau <- 1e-6
  zeta <- 0.1
  gamma <- 0.99
  # initial point
  obj_value <- NULL
  w_SCA <- rep(1/N, N)
  for (k in 1:max_iter) {
    # compute parameters
    gA <- compute_gA(w_SCA, Sigma)
    g <- gA$g
    A <- gA$A
    Q <- 2 * t(A) %*% A + tau*diag(N) # crossprod(A) = t(A) %*% A
    q <- 2 * t(A) %*% g - Q %*% w_SCA
    obj_value <- c(obj_value, sum(g^2))
    # solve problem with CVXR
    w_ <- Variable(N)
    prob <- Problem(Minimize(0.5*quad_form(w_, Q) + t(q) %*% w_),
                    constraints = list(sum(w_) == 1))
    result <- solve(prob)
    w_ <- as.vector(result$getValue(w_))
    # solve the problem with solve.QP()
    w__ <- solve.QP(Q, -q, matrix(1, N, 1), 1, meq = 1)$solution
    # solve problem in closed form
    C <- matrix(1, 1, N)
    c <- 1
    CinvQ <- C %*% solve(Q)
    lmd <- solve(CinvQ %*% t(C), -(CinvQ %*% q + c))
    w___ <- solve(Q, -(q + t(C) %*% lmd))
    #sanity checks for different solvers
    if ((err <- norm(w_ - w__, "2")/norm(w_, "2")) > 1e-2)
      cat("CVXR and solve.QP do not match:", err, "\n")
    if ((err <- norm(w_ - w___, "2")/norm(w_, "2")) > 1e-2)
      cat("Closed-form solution and CVXR do not match:", err, "\n")
    # next w
    gamma <- gamma*(1 - zeta*gamma)
    w_SCA_prev <- w_SCA
    w_SCA <- w_SCA + gamma*(w_ - w_SCA)
    # stopping criterion
    #if ((norm(w-w_prev, "2")/norm(w_prev, "2")) < 1e-6)
    # break
    if (k>1 && abs((obj_value[k]-obj_value[k-1])/obj_value[k-1]) < 1e-1)
      break
  }
  # cat("Number of iterations:", k)
  # 
  # plot(obj_value, type = "b",
  #      main = "Convergence of SCA", xlab = "iteration", ylab = "objective value")
  # 
  w_all <- cbind(w_all, "risk-parity-SCA" = w_SCA)
  
  # mu <- colMeans(X_)
  # Sigma <- cov(X_)
  rpp_mu <- riskParityPortfolio(Sigma, mu = mu, lmd_mu = 1e-6, formulation = "rc-double-index")
  w_all <- cbind(w_all,"risk-parity-mu" = rpp_mu$w)
  # 
  return(w_all)
}

###### 3.2.2 weekly,monthly,& 1\2\3weekly + rolling ####################
######## 3.2.2.1 weekly and monthly rolling #######################
# input:
# 1&2. X_log&X_lin are xts of log or linear return
# 3. on = c('weeks','months')
# 4. start = c(0,1),0 stands for block, 1 stands for fix start
rollingWM <- function(X_log,X_lin,on,start){
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
######## 3.2.2.2 endpoints is 1/2/3 week's Friday rolling #######################
# input:
# 1&2. X_log&X_lin are xts of log or linear return
# 3. on = c('1','2','3')
# 4. start = c(0,1),0 stands for block, 1 stands for fix start
rollingFri <- function(X_log,X_lin,week,start){
  T_trn = 60
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

####### 3.4 calculate rank #################################
CalculateRank <- function(perf){
  ratio = matrix(ncol=ncol(perf),nrow=nrow(perf))
  colnames(ratio)=colnames(perf)
  rownames(ratio)=rownames(perf)
  for (i in 1: nrow(perf)){
    print(i)
    if(rownames(perf)[i]=="cov"||
       rownames(perf)[i]=="Annualized Std Dev"||
       rownames(perf)[i]=="Downside Deviation (MAR = 0%)"||
       rownames(perf)[i]=="d ratio"||
       rownames(perf)[i]=="Drawdown Deviation"||
       rownames(perf)[i]=="Worst Downside"||
       rownames(perf)[i]=="Mean absolute deviation"||
       rownames(perf)[i]=="Ulcer Index"||
       rownames(perf)[i]=="continuous loss days"){
      s =sort(as.numeric(perf[i,]),decreasing= F, index.return = TRUE)
    }
    else{s =sort(as.numeric(perf[i,]), decreasing= T, index.return = TRUE)}
    tmp=perf[,s$ix]
    for (j in 1:ncol(tmp)){
      ratio[i,which(colnames(ratio)==colnames(tmp[j]))]=j}
  }
  return(ratio)
}
######## 3.4.1 rolling1 #################################
WeightUnif <- function(strategies_log_roll){
  ### Uniform Portfolio ###
  N = ncol(strategies_log_roll)
  w_unif <- as.matrix(rep(1/N, N))
  
  rownames(w_unif) <- colnames(strategies_log_roll)
  colnames(w_unif) <- c("unif")
  round(w_unif, digits = 2)
  return(w_unif)
}
rollingUnif <- function(X_log,X_lin,allSig,rt,start){
  
  t <- nrow(X_log)
  strategy = list()
  perf.list = list()
  rank.list = list()
  T_trn = which(index(X_log)=='2018-12-28 CST')
  
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):t, ]
  
  rebal_indice1 <- (c(0,endpoints(X_log_trn, on = 'month')[2:7])+1)
  rebal_indice2 <- T_trn + endpoints(X_log_tst, on = 'month')
  
  
  result =matrix(ncol=1)[-1,]
  
  for (i in 1:(length(rebal_indice2)-1)){
    #Training
    if(start==1 ||i ==1){
      print(i)
      print('start')
      X_log.day <- X_log[rebal_indice1[1]:rebal_indice2[i],]
      name1 = 1
      name2 = paste(rownames(as.matrix(X_log.day))[nrow(X_log.day)],'15:15:00',sep = " ")
      X <- rt[1:which(rownames(as.matrix(rt))==name2),]
      Sig <- allSig[1:which(rownames(as.matrix(allSig))==name2),]}
    else{
      print('block')
      X_log.day <- X_log[rebal_indice1[i]:rebal_indice2[i],]
      if(i ==5){name1 = "2018-05-02 09:01:00"}
      else{name1 = paste(rownames(as.matrix(X_log[(rebal_indice1[i]-1),])),'21:01:00',sep = " ")}
      name2 = paste(rownames(as.matrix(X_log.day))[nrow(X_log.day)],'15:15:00',sep = " ")
      X <- rt[which(rownames(as.matrix(rt))==name1):which(rownames(as.matrix(rt))==name2),]
      Sig <- allSig[which(rownames(as.matrix(allSig))==name1):which(rownames(as.matrix(allSig))==name2),]}
    
    # get performance and rank
    perf = performance(Sig,X,X_log.day)
    rank.perf = CalculateRank(perf)
    
    perf.list=c(perf.list,perf)
    rank.list=c(rank.list,rank.perf)
    
    # select top
    indx <- c("Annualized Std Dev","Annualized Sharpe (Rf=0%)","Drawdown Deviation",
              "Kelly Ratio","Worst Drawdown","gain ratio")
    
    top5<- apply(rank.perf[indx,1:5], 1, function(x){colnames(rank.perf)[x]})
    stra <- unique(as.vector(top5))
    
    #get weights
    w <- WeightUnif(X_log.day[,stra])
    
    # get test return
    tst <- X_lin[(rebal_indice2[i]+1):min(rebal_indice2[i+1],nrow(X_lin)),stra]
    #tst <- X_lin[(rebal_indice2[i]+1):nrow(X_lin),stra]
    # calculate mean returns as result for each days
    result <- c(result,as.matrix(tst)%*%w)
    strategy <- c(strategy,list(stra))
  }
  
  result = as.matrix(result)
  rownames(result) = rownames(as.matrix(X_log))[(T_trn+1):t]
  colnames(result) = c("Unif")
  
  
  return(list('result'=result,
              'strategy'=strategy,
              'performances' = perf.list,
              'rank.perf' = rank.list))
}


####### Indices #################################################
# Function indices returns a matrix contaning 13 ratios, which indicates the 
# performance of input return
## input: a xts of returns
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


####### 手数 #################################
####### Based on Strategy #################################
shoushuStra <- function(x,monthPrice){
  x <- round(x,digits = 4)
  numPf = 8
  num <- ncol(x)/numPf
  
  shoushu<-c()
  all_cap <- c()
  
  for (k in 1:num) {
    ## 交易单位 * 价格
    monthCoe <- matrix(0,nrow(x),numPf,dimnames = list(rownames(x)))
    multi <- coe[colnames(monthPrice),1]*monthPrice[k,]
    for (i in 1:length(multi)) {
      monthCoe[substr(rownames(monthCoe),1,2) == names(multi)[i],] = multi[i]
    }
    
    #计算手数
    x0 <- x[,((k-1)*numPf+1):(k*numPf)]
    Capital <- matrix(0,nrow(x0),ncol(x0),dimnames = list(rownames(x0),colnames(x0)))
    T_ <- x0[substr(rownames(x0),1,2) == "T_",]
    for (j in 1:numPf) {
      stra <- rownames(T_)[which(T_[,j] == min(T_[T_[,j] > 0, j]))][1]
      initialCapital <- monthPrice[k,"T_"]*10000*1/x0[stra,j]
      Capital[,j] <- initialCapital*x0[,j]
    }
    all_cap <- cbind(all_cap,round(Capital/monthCoe)*monthCoe)
    shoushu <- cbind(shoushu,Capital/monthCoe)
  }
  return(list("capital" = all_cap , "shoushu" = shoushu))
}
####### Based on contract #################################
shoushuStraC <- function(x,monthPrice,numT){
  x<-x[,-which(substr(colnames(x),1,4) == "Mark")]
  x <- round(x,digits = 4)
  for (i in 1:nrow(x)) {
    if(nchar(rownames(x)[i])==1){
      rownames(x)[i] = paste(rownames(x)[i],"_",sep="")
    }
  }
  numPf = 7
  num <- ncol(x)/numPf
  
  shoushu<-c()
  all_cap <- c()
  
  # contract <- c("a_","al","hc","i_","IF","j_","PP","rb","SR","T_","ZC","zn","jm")
  # accum <- matrix(0, nrow = length(contract), ncol = numPf)
  # rownames(accum) <- contract
  
  for (k in 1:num) {
    ## 交易单位 * 价格
    monthCoe <- matrix(0,nrow(x),numPf,dimnames = list(rownames(x)))
    multi <- coe[colnames(monthPrice),1]*monthPrice[k,]
    for (i in 1:length(multi)) {
      monthCoe[rownames(monthCoe) == names(multi)[i],] = multi[i]
    }
    
    #计算手数
    ## 取 一次 rolling weights
    x0 <- x[,((k-1)*numPf+1):(k*numPf)]
    ## 计算品种 weights
    # for(i in 1:length(contract)){
    #   accum[i,] <- apply(x0[substr(rownames(x0),1,2) == contract[i],],2,sum)
    # }
    # for(j in 1:ncol(accum)){
    #   accum[,j] <- accum[,j]/accum["T_",j]
    # }
    # ## 计算策略在品种中weights比例
    # rate <- matrix(0,nrow=nrow(x0),ncol = ncol(x0))
    # rownames(rate) = rownames(x0)
    # for (i in 1:length(contract)) {
    #   for (j in 1:ncol(rate)) {
    #     rate[substr(rownames(rate),1,2) == contract[i],j] <- x0[substr(rownames(x0),1,2) == contract[i],j]/accum[i,j]
    #   }
    # }
    # 
    ## 计算 T 市值
    Capital <- matrix(0,nrow(x0),ncol(x0),dimnames = list(rownames(x0),colnames(x0)))
    T_ <- x0["T_",]
    for (j in 1:numPf) {
      # stra <- rownames(T_)[which(T_[,j] == min(T_[T_[,j] > 0, j]))][1]
      initialCapital <- monthPrice[k,"T_"]*10000*numT/T_[,j]
      #T_Cap <- initialCapital*T_[,j]
      Capital[,j] <- initialCapital*x0[,j]
    }
    all_cap <- cbind(all_cap,round(Capital/monthCoe)*monthCoe)
    shoushu <- cbind(shoushu,Capital/monthCoe)
  }
  return(list("capital" = all_cap , "shoushu" = shoushu))
}
