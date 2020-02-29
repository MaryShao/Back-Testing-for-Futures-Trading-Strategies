library(riskParityPortfolio)
library(quadprog) 
library(xts) 
library(quantmod) 
library(PerformanceAnalytics)
library(CVXR)

## Markowitz  & GMVP
portolioMarkowitz <- function(mu, Sigma, lmd = 0.5) {
  w <- Variable(nrow(Sigma))
  prob <- Problem(Maximize(t(mu) %*% w - lmd*quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}
portolioGMVP <- function(Sigma) {
  w <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

##### sca
compute_gA <- function(w, Sigma) {
  N <- length(w)
  g <- rep(NA, N^2)
  A <- matrix(NA, N^2, N)
  for (i in 1:N) {
    Mi <- matrix(0, N, N)
    Mi[i, ] <- Sigma[i, ]
    for (j in 1:N) {
      Mj <- matrix(0, N, N)
      Mj[j, ] <- Sigma[j, ]
      #g[i + (j-1)*N] <- t(w) %*% (Mi - Mj) %*% w
      g[i + (j-1)*N] <- w[i]*(Sigma[i, ] %*% w) - w[j]*(Sigma[j, ] %*% w)
      A[i + (j-1)*N, ] <- (Mi + t(Mi) - Mj - t(Mj)) %*% w
      #A[i + (j-1)*N, ] <- (Sigma[i, ] %*% w
    }
  }
  # # this is much faster:
  # wSw <- w * (Sigma %*% w)
  # g <- rep(wSw, times = N) - rep(wSw, each = N) # N^2 different g_{i,j}
  return(list(g = g, A = A))
}

##### weights (w_all)
wght<-function(X_log_trn){
  mu <- colMeans(X_log_trn)
  Sigma <- cov(as.matrix(X_log_trn))
  
  w_Markowitz <- portolioMarkowitz(mu, Sigma)
  w_GMVP <- portolioGMVP(Sigma)
  
  w_all <- cbind(w_GMVP, w_Markowitz)
  rownames(w_all) <- colnames(X_log_trn)
  colnames(w_all) <- c("GMVP", "Markowitz")
  
  
  # ## risk-parity-diag
  # w_diag <- 1/sqrt(diag(Sigma))
  # w_diag <- w_diag/sum(w_diag)
  # # plot
  # w_all <- cbind(w_all, "risk-parity-diag" = w_diag)
  # 
  
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
  # plot
  w_all <- cbind(w_all, "risk-parity-gen-solver" = w_gen_solver)
  
  # ## risk-parity-package
  # #devtools::install_github("dppalomar/riskParityPortfolio")
  # library(riskParityPortfolio)
  # rpp <- riskParityPortfolio(Sigma)
  # w_all <- cbind(w_all, "risk-parity-package" = rpp$w)
  # 
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
  
  mu <- colMeans(X_)
  Sigma <- cov(X_)
  rpp_mu <- riskParityPortfolio(Sigma, mu = mu, lambda = 1e-7, formulation = "rc-double-index")
  w_all <- cbind(w_all,"risk-parity-mu" = rpp_mu$w)
  
  return(w_all)
}

#####  rolling windows
rolling <- function(X_log,on,start){
  T_trn <- 15
  t <- nrow(X_log)
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):t, ]
  # X_lin_trn <- X_lin[1:T_trn, ]
  # X_lin_tst <- X_lin[(T_trn+1):t, ]
  # Rank_rolling <- X_log
  # Rank_rolling[] <- NA
  rebal_indices <- T_trn + endpoints(X_log_tst, on = on)
  #rebal_indices2 <- T_trn + endpoints(X_log_tst, on = "months")
  # index(X_log)[rebal_indices]
  lookback <- 15
  result =matrix(ncol=6)[-1,]
  
  for (i in 1:(length(rebal_indices)-1)){
    #Training
    if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
    else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
    #get weights
    w_all <- wght(X_)
    # get test return
    tst <- X_lin[(rebal_indices[i]+1):min(rebal_indices[i+1],nrow(X_log)),] 
    # calculate mean returns as result for each days
    result <- rbind(result,as.matrix(tst)%*%w_all)
  }
  
  rownames(result) = rownames(X_log)[(T_trn+1):t]
  colnames(result) = colnames(w_all)
  
  return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
              "MDD"=maxDrawdown(result, weights = NULL, geometric = FALSE, invert = TRUE),
               "result" = result))
}

on = "weeks"
week1819start<-rolling(X_log,on,start = 1)
week1819block<-rolling(X_log,on,start = 0)

on = "months"
month1819start<-rolling(X_log,on,start = 1)
month1819block<-rolling(X_log,on,start = 0)


par(mfrow=c(2,2))
chart.CumReturns(week1819start$result, main = "18-19 performance of 1-week fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(week1819block$result, main = "18-19 performance of 1-week fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month1819start$result, main = "18-19 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month1819block$result, main = "18-19 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)


######################################### 2019
X_log19tst = X_log[229:357,]
#rolling windows 19 tst
on = "weeks"
week19start<-rolling(X_log19tst,on,start = 1)
week19block<-rolling(X_log19tst,on,start = 0)


week19startMU<-rollingMU(X_log19tst,on,start = 1)
#rolling windows 19 tst month
on = "months"
month19start<-rolling(X_log19tst,on,start = 1)
month19block<-rolling(X_log19tst,on,start = 0)

month19startMU<-rollingMU(X_log19tst,on,start = 1)
on = 'quarters'
quater19block<-rollingMU(X_log19tst,on,start = 0)

chart.CumReturns(week19start$result, main = "2019 performance of 1-week fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(week19block$result, main = "2019 performance of 1-week fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month19start$result, main = "2019 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month19block$result, main = "2019 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)


colnames(quater19block$result) = c("1","2","3","4","5","6","7","8","9","10")

par(mfrow=c(3,1))
chart.CumReturns(week19startMU$result, main = "2019 performance of weeks fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow10equal)
chart.CumReturns(month19startMU$result, main = "2019 performance of months fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow10equal)
chart.CumReturns(quater19block$result, main = "2019 performance of quarters fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow10equal)


####################################
w_all <- wght(X_log[1:250,])
mu <- colMeans(X_log[224:243,])
Sigma <- cov(as.matrix(X_log[224:243,]))
## risk-parity-diag
w_diag <- 1/sqrt(diag(Sigma))
w_diag <- w_diag/sum(w_diag)
# plot
w_all <- cbind(w_all, "risk-parity-diag" = w_diag)
## risk-parity-package
#devtools::install_github("dppalomar/riskParityPortfolio")
library(riskParityPortfolio)
rpp <- riskParityPortfolio(Sigma)
w_all <- cbind(w_all, "risk-parity-package" = rpp$w)




return <-as.matrix( X_lin)%*%w_all
chart.CumReturns(return, main = "1819 noRolling",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)

indices <- function(return){
  res<-rbind(table.AnnualizedReturns(return),
             SharpeRatio(return, Rf = 0, p = 0.95, FUN = "StdDev"),
             maxDrawdown(return, weights = NULL, geometric = TRUE, invert = TRUE),
             CDD(return, weights = NULL, geometric = TRUE, invert = TRUE, p = 0.95),
             CDD(return, weights = NULL, geometric = TRUE, invert = TRUE, p = 0.5),
             VaR(return, p = 0.95,method = "historical", portfolio_method ="single",invert = TRUE),
             VaR(return, p = 0.99,method = "historical", portfolio_method ="single",invert = TRUE))
  return(res)
}

write.csv(rbind(t(indices(week1819start$result)),
                t(indices(week1819block$result)),
                t(indices(month1819start$result)),
                t(indices(month1819block$result)),
                t(indices(week19start$result)),
                t(indices(week19block$result)),
                t(indices(month19start$result)),
                t(indices(month19block$result))),"indices.csv")

write.csv(rbind(
                t(indices(week19startMU$result)),
                t(indices(month19startMU$result)),
                t(indices(quater19block$result))),"indicesMU.csv")




write.csv(t(res),"noRolling indices.csv")



rollingMU <- function(X_log,on,start){
  T_trn <- 20
  t <- nrow(X_log)
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):t, ]
  # X_lin_trn <- X_lin[1:T_trn, ]
  # X_lin_tst <- X_lin[(T_trn+1):t, ]
  # Rank_rolling <- X_log
  # Rank_rolling[] <- NA
  rebal_indices <- T_trn + endpoints(X_log_tst, on = on)
  #rebal_indices2 <- T_trn + endpoints(X_log_tst, on = "months")
  # index(X_log)[rebal_indices]
  lookback <- 20  
  result =matrix(ncol=10)[-1,]
  
  for (i in 1:(length(rebal_indices)-1)){
    #Training
    if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
    else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
    #get weights
    mu <- colMeans(X_)
    Sigma <- cov(X_)
    w_all <- cbind(riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-1,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-2,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-3,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-4,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-5,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-6,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-7,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-8,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-9,  formulation = "rc-double-index")$w,
                   riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-10,  formulation = "rc-double-index")$w)
    # get test return
    tst <- X_lin[(rebal_indices[i]+1):min(rebal_indices[i+1],nrow(X_log)),] 
    # calculate mean returns as result for each days
    result <- rbind(result,as.matrix(tst)%*%w_all)
  }
  
  rownames(result) = rownames(X_log)[(T_trn+1):t]
  colnames(result) = c("1","2","3","4","5","6","7","8","9","10")
  
  return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
              "MDD"=maxDrawdown(result, weights = NULL, geometric = FALSE, invert = TRUE),
              "result" = result))
}




################################################## contract ##############################################
X_logContract <- Contract.Dailylogreturn
X_linContract <- Contract.Dailylinreturn

mu <- colMeans(X_logContract[1:20,])
Sigma <- cov(X_logContract[1:20,])
w_all <- riskParityPortfolio(Sigma, mu = mu, lmd_mu =1e-7,  formulation = "rc-double-index")$w
               
# get test return
tst <- X_lin[(rebal_indices[i]+1):min(rebal_indices[i+1],nrow(X_log)),] 
# calculate mean returns as result for each days
result <- rbind(result,as.matrix(tst)%*%w_all)
## 18-19
on = "weeks"
week1819startC<-rolling(X_logContract,on,start = 1)
week1819blockC<-rolling(X_logContract,on,start = 0)

on = "months"
month1819startC<-rolling(X_logContract,on,start = 1)
month1819blockC<-rolling(X_logContract,on,start = 0)


par(mfrow=c(2,2))
chart.CumReturns(week1819startC$result, main = "Contract 18-19 performance of 1-week fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(week1819blockC$result, main = "Contract 18-19 performance of 1-week fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month1819startC$result, main = "Contract 18-19 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month1819blockC$result, main = "Contract 18-19 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)

######################################### 2019
X_log19tst = X_log[229:357,]


#rolling windows 19 tst
on = "weeks"
week19start<-rolling(X_log19tst,on,start = 1)
week19block<-rolling(X_log19tst,on,start = 0)


week19startMU<-rollingMU(X_log19tst,on,start = 1)
#rolling windows 19 tst month
on = "months"
month19start<-rolling(X_log19tst,on,start = 1)
month19block<-rolling(X_log19tst,on,start = 0)

month19startMU<-rollingMU(X_log19tst,on,start = 1)
on = 'quarters'
quater19block<-rollingMU(X_log19tst,on,start = 0)

chart.CumReturns(week19start$result, main = "2019 performance of 1-week fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(week19block$result, main = "2019 performance of 1-week fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month19start$result, main = "2019 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)
chart.CumReturns(month19block$result, main = "2019 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow6equal)


colnames(quater19block$result) = c("1","2","3","4","5","6","7","8","9","10")

par(mfrow=c(3,1))
chart.CumReturns(week19startMU$result, main = "2019 performance of weeks fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow10equal)
chart.CumReturns(month19startMU$result, main = "2019 performance of months fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow10equal)
chart.CumReturns(quater19block$result, main = "2019 performance of quarters fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow10equal)




