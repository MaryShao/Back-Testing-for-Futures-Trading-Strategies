
library(riskParityPortfolio)
library(quadprog) 
library(xts) 
library(quantmod) 
library(PerformanceAnalytics)
library(CVXR)

indices <- function(return){
    res<-rbind(table.AnnualizedReturns(return),
               SharpeRatio(return, Rf = 0, p = 0.95, FUN = "StdDev"),
               SortinoRatio(return,MAR =0),
               maxDrawdown(return, weights = NULL, geometric = FALSE, invert = TRUE),
               CDD(return, weights = NULL, geometric = FALSE, invert = TRUE, p = 0.95),
               CDD(return, weights = NULL, geometric = FALSE, invert = TRUE, p = 0.5),
               VaR(return, p = 0.95,method = "historical", portfolio_method ="single",invert = TRUE),
               VaR(return, p = 0.99,method = "historical", portfolio_method ="single",invert = TRUE))
    return(res)
}
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


######################### read data #############################
#X_log
X_log = Strategy.Dailylogreturn1
rownames(X_log) = substr(rownames(X_log),1,10)
#X_lin
X_lin = Strategy.Dailylinreturn1
rownames(X_lin) = substr(rownames(X_lin),1,10)



####################################rolling month ############################################

rolling <- function(X_log,X_lin,start,T_trn,on){
    t <- nrow(X_log)
    X_log_trn <- X_log[1:T_trn, ]
    X_log_tst <- X_log[(T_trn+1):t, ]
    rebal_indices <- T_trn + endpoints(X_log_tst, on = on)
    lookback <- T_trn
    result =matrix(ncol=8)[-1,]
    all_weights=list()
    for (i in 1:(length(rebal_indices)-1)){
        #Training
        if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
        else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
        #get weights
        w_all <- wght(X_)
        
        # get test return
        if(i == (length(rebal_indices)-1)) {tst <- X_lin[(rebal_indices[i]+1):nrow(X_log),] }
        else{tst <- X_lin[(rebal_indices[i]+1):rebal_indices[i+1],] }
        print(rownames(tst))
        result <- rbind(result,as.matrix(tst)%*%w_all)
        
        all_weights <- c(all_weights,list(w_all))
    }
    
    rownames(result) = rownames(X_log)[(T_trn+1):t]
    colnames(result) = colnames(w_all)
    
    return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
                "MDD"=maxDrawdown(result, weights = NULL, geometric = FALSE, invert = TRUE),
                "result" = result,
                "weights" = all_weights))
}


T_trn = 15
on = "weeks"
week1819start<-rolling(X_log,X_lin,start = 1,T_trn,on)
week1819block<-rolling(X_log,X_lin,start = 0,T_trn,on)
T_trn = 60
on = "months"
month1819start<-rolling(X_log,X_lin,start = 1,T_trn,on)
month1819block<-rolling(X_log,X_lin,start = 0,T_trn,on)

par(mfrow=c(2,2))
chart.CumReturns(week1819start$result, main = "18-19 performance of 1-week fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal,geometric = FALSE)
chart.CumReturns(week1819block$result, main = "18-19 performance of 1-week fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month1819start$result, main = "18-19 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month1819block$result, main = "18-19 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


######################################### 2019

X_log19tstWeek = X_log[229:357,]
X_lin19tstWeek = X_lin[229:357,]
T_trn = 15
on = "weeks"
week19start<-rolling(X_log19tstWeek,X_lin19tstWeek,start = 1,T_trn,on)
week19block<-rolling(X_log19tstWeek,X_lin19tstWeek,start = 0,T_trn,on)


X_log19tstMonth = X_log[163:357,]
X_lin19tstMonth = X_lin[163:357,]
T_trn = 242
on = "months"
month19start<-rolling(X_log,X_lin,start = 1,T_trn,on)
month19block<-rolling(X_log,X_lin,start = 0,T_trn,on)


chart.CumReturns(week19start$result, main = "2019 performance of 1-week fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(week19block$result, main = "2019 performance of 1-week fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month19start$result, main = "2019 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month19block$result, main = "2019 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(week1819start$result)),
                t(indices(week1819block$result)),
                t(indices(month1819start$result)),
                t(indices(month1819block$result)),
                t(indices(week19start$result)),
                t(indices(week19block$result)),
                t(indices(month19start$result)),
                t(indices(month19block$result))),"indicesWeekMonth.csv")


write.csv(week1819start$weights,"week1819start.csv")
write.csv(week1819block$weights,"week1819block.csv")
write.csv(month1819start$weights,"month1819start.csv")
write.csv(month1819block$weights,"month1819block.csv")
write.csv(week19start$weights,"week19start.csv")
write.csv(week19block$weights,"week19block.csv")
write.csv(month19start$weights,"month19start.csv")
write.csv(month19block$weights,"month19block.csv")

############################################# rolling 3rd Fri #########################################################
rolling3Fri <- function(X_log,X_lin,start,T_trn){
    t <- nrow(X_log)
    X_log_trn <- X_log[1:T_trn, ]
    X_log_tst <- X_log[(T_trn+1):t, ]
    
    end_20month<-c()
    end_week<-endpoints(X_log_tst, on = "weeks")
    end_month<-endpoints(X_log_tst, on = "months")
    for (i in end_month){
        ind = which(end_week==max(end_week[end_week<=i]))-1
        if(ind>0)
        {end_20month<-c(end_20month,end_week[ind])}
    }
    rebal_indices <- T_trn + c(0,end_20month) 
    if(rebal_indices[length(rebal_indices)]<nrow(X_log)){
        rebal_indices = c(rebal_indices, nrow(X_log))
    }
    lookback <- T_trn
    result =matrix(ncol=8)[-1,]
    all_weights=list()
    
    for (i in 1:(length(rebal_indices)-1)){
        #Training
        if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
        else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
        #get weights
        w_all <- wght(X_)
        # get test return
        if(i == (length(rebal_indices)-1)) {tst <- X_lin[(rebal_indices[i]+1):nrow(X_log),] }
        else{tst <- X_lin[(rebal_indices[i]+1):rebal_indices[i+1],] }
        result <- rbind(result,as.matrix(tst)%*%w_all)
        print(rownames(tst))
        all_weights <- c(all_weights,list(w_all))
    }
    
    rownames(result) = rownames(X_log)[(T_trn+1):t]
    colnames(result) = colnames(w_all)
    
    return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
                "MDD"=maxDrawdown(result, weights = NULL, geometric = FALSE, invert = TRUE),
                "result" = result,
                "weights" = all_weights))
}

T_trn = 60
thirdFRI1819start<-rolling3Fri(X_log,X_lin,start = 1,T_trn)
thirdFRI1819block<-rolling3Fri(X_log,X_lin,start = 0,T_trn)


par(mfrow=c(2,2))
chart.CumReturns(thirdFRI1819start$result, main = "3rd Fri 18-19 fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal,geometric = FALSE)
chart.CumReturns(thirdFRI1819block$result, main = "3rd Fri 18-19 fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


########### 2019
X_log19tst = X_log[163:357,]
X_lin19tst = X_lin[163:357,]

T_trn =242
thirdFRI19start<-rolling3Fri(X_log,X_lin,start = 1,T_trn)
thirdFRI19block<-rolling3Fri(X_log,X_lin,start = 0,T_trn)

chart.CumReturns(thirdFRI19start$result, main = "3rd Fri 2019 fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(thirdFRI19block$result, main = "3rd Fri 2019 fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(thirdFRI1819start$result)),
                t(indices(thirdFRI1819block$result)),
                t(indices(thirdFRI19start$result)),
                t(indices(thirdFRI19block$result))),"indices3Fri.csv")

write.csv(thirdFRI1819start$weights,"thirdFRI1819start.csv")
write.csv(thirdFRI1819block$weights,"thirdFRI1819block.csv")
write.csv(thirdFRI19start$weights,"thirdFRI19start.csv")
write.csv(thirdFRI19block$weights,"thirdFRI19block.csv")

########################################## rolling 2nd Fri###########################################################
rolling2Fri <- function(X_log,X_lin,start,T_trn){
    t <- nrow(X_log)
    X_log_trn <- X_log[1:T_trn, ]
    X_log_tst <- X_log[(T_trn+1):t, ]

    end_20month<-c()
    end_week<-endpoints(X_log_tst, on = "weeks")
    end_month<-endpoints(X_log_tst, on = "months")
    for (i in end_month){
        ind = which(end_week==max(end_week[end_week<=i]))-2
        if(ind>0)
        {end_20month<-c(end_20month,end_week[ind])}
    }
    rebal_indices <- T_trn + c(0,end_20month) 
    if(rebal_indices[length(rebal_indices)]<nrow(X_log)){
        rebal_indices = c(rebal_indices, nrow(X_log))
    }

    lookback <- T_trn
    result =matrix(ncol=8)[-1,]
    all_weights=list()
    
    for (i in 1:(length(rebal_indices)-1)){
        #Training
        if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
        else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
        #get weights
        w_all <- wght(X_)
        # get test return
        if(i == (length(rebal_indices)-1)) {tst <- X_lin[(rebal_indices[i]+1):nrow(X_log),] }
        else{tst <- X_lin[(rebal_indices[i]+1):rebal_indices[i+1],] }
        print(rownames(tst))
        
        result <- rbind(result,as.matrix(tst)%*%w_all)
        
        all_weights <- c(all_weights,list(w_all))
    }
    
    rownames(result) = rownames(X_log)[(T_trn+1):t]
    colnames(result) = colnames(w_all)
    
    return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
                "MDD"=maxDrawdown(result, weights = NULL, geometric = FALSE, invert = TRUE),
                "result" = result,
                "weights" = all_weights))
}

T_trn = 60
secFRI1819start<-rolling2Fri(X_log,X_lin,start = 1,T_trn)
secFRI1819block<-rolling2Fri(X_log,X_lin,start = 0,T_trn)

par(mfrow=c(2,2))
chart.CumReturns(secFRI1819start$result, main = "2nd Fri 18-19 fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(secFRI1819block$result, main = "2rd Fri 18-19 fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


###### 2019
X_log19tst = X_log[163:357,]
X_lin19tst = X_lin[163:357,]

T_trn = 242
secFRI19start<-rolling2Fri(X_log,X_lin,start = 1,T_trn)
secFRI19block<-rolling2Fri(X_log,X_lin,start = 0,T_trn)

chart.CumReturns(secFRI19start$result, main = "2nd Fri 2019 fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(secFRI19block$result, main = "2nd Fri 2019 fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(secFRI1819start$result)),
                t(indices(secFRI1819block$result)),
                t(indices(secFRI19start$result)),
                t(indices(secFRI19block$result))),"indices2Fri.csv")

write.csv(secFRI1819start$weights,"secFRI1819start.csv")
write.csv(secFRI1819block$weights,"secFRI1819block.csv")
write.csv(secFRI19start$weights,"secFRI19start.csv")
write.csv(secFRI19block$weights,"secFRI19block.csv")

############################################### rolling 1st Fri ###################################################
rolling1Fri <- function(X_log,X_lin,start,T_trn){
    t <- nrow(X_log)
    X_log_trn <- X_log[1:T_trn, ]
    X_log_tst <- X_log[(T_trn+1):t, ]
    # X_lin_trn <- X_lin[1:T_trn, ]
    # X_lin_tst <- X_lin[(T_trn+1):t, ]
    # Rank_rolling <- X_log
    # Rank_rolling[] <- NA
    end_20month<-c()
    end_week<-endpoints(X_log_tst, on = "weeks")
    end_month<-endpoints(X_log_tst, on = "months")
    for (i in end_month){
        ind = which(end_week==max(end_week[end_week<=i]))-3
        if(ind>0)
        {end_20month<-c(end_20month,end_week[ind])}
    }
    rebal_indices <- T_trn + c(0,end_20month) 
    if(rebal_indices[length(rebal_indices)]<nrow(X_log)){
        rebal_indices = c(rebal_indices, nrow(X_log))
    }
    #rebal_indices2 <- T_trn + endpoints( X_logContract[243:357,], on = "weeks")
    # index(X_log)[rebal_indices]
    lookback <- T_trn
    result =matrix(ncol=8)[-1,]
    all_weights=list()
    
    for (i in 1:(length(rebal_indices)-1)){
        #Training
        #print(i)
        if(start==1){X_ <- X_log[(rebal_indices[1]-lookback+1):rebal_indices[i], ]}
        else{X_ <- X_log[(rebal_indices[i]-lookback+1):rebal_indices[i], ]}
        #get weights
        w_all <- wght(X_)
        # get test return
        if(i == (length(rebal_indices)-1)) {tst <- X_lin[(rebal_indices[i]+1):nrow(X_log),] }
        else{tst <- X_lin[(rebal_indices[i]+1):rebal_indices[i+1],] }
        print(rownames(tst))
        # calculate mean returns as result for each days
        result <- rbind(result,as.matrix(tst)%*%w_all)
        all_weights <- c(all_weights,list(w_all))
    }
    
    rownames(result) = rownames(X_log)[(T_trn+1):t]
    colnames(result) = colnames(w_all)
    
    return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
                "MDD"=maxDrawdown(result, weights = NULL, geometric = FALSE, invert = TRUE),
                "result" = result,
                "weights" = all_weights))
}

T_trn = 60
firstFRI1819start<-rolling1Fri(X_log,X_lin,start = 1,T_trn)
firstFRI1819block<-rolling1Fri(X_log,X_lin,start = 0,T_trn)

par(mfrow=c(2,2))
chart.CumReturns(firstFRI1819start$result, main = "1st Fri 18-19 fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(firstFRI1819block$result, main = "1st Fri 18-19 fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


###### 2019
X_log19tst = X_log[163:357,]
X_lin19tst = X_lin[163:357,]

T_trn = 242
firstFRI19start<-rolling1Fri(X_log,X_lin,start = 1,T_trn)
firstFRI19block<-rolling1Fri(X_log,X_lin,start = 0,T_trn)

chart.CumReturns(firstFRI19start$result, main = "1st Fri 2019 fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(firstFRI19block$result, main = "1st Fri 2019 fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(firstFRI1819start$result)),
                t(indices(firstFRI1819block$result)),
                t(indices(firstFRI19start$result)),
                t(indices(firstFRI19block$result))),"indices1Fri.csv")

write.csv(firstFRI1819start$weights,"firstFRI1819start.csv")
write.csv(firstFRI1819block$weights,"firstFRI1819block.csv")
write.csv(firstFRI19start$weights,"firstFRI19start.csv")
write.csv(firstFRI19block$weights,"firstFRI19block.csv")

###############################################################  Contract  #################################################################

############## data ##############
#X_logC
X_logC = Contract.Dailylogreturn
rownames(X_logC) = substr(rownames(X_logC),1,10)
#X_linC
X_linC = Contract.Dailylinreturn
rownames(X_linC) = substr(rownames(X_linC),1,10)


################################### rolling months/weeks ##############################################
## 18-19
on = "weeks"
T_trn = 15
week1819startC<-rolling(X_logC,X_linC,start = 1,T_trn,on)
week1819blockC<-rolling(X_logC,X_linC,start = 0,T_trn,on)

on = "months"
T_trn = 60
month1819startC<-rolling(X_logC,X_linC,start = 1,T_trn,on)
month1819blockC<-rolling(X_logC,X_linC,start = 0,T_trn,on)


par(mfrow=c(2,2))
chart.CumReturns(week1819startC$result, main = "Contract 18-19 performance of 1-week fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(week1819blockC$result, main = "Contract 18-19 performance of 1-week fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month1819startC$result, main = "Contract 18-19 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month1819blockC$result, main = "Contract 18-19 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


#### 2019
X_log19tstC = X_logC[229:357,]
X_lin19tstC = X_linC[229:357,]
on = "weeks"
T_trn = 15
week19startC<-rolling(X_log19tstC,X_lin19tstC,start = 1,T_trn,on)
week19blockC<-rolling(X_log19tstC,X_lin19tstC,start = 0,T_trn,on)


X_log19tstC = X_logC[184:357,]
X_lin19tstC = X_linC[184:357,]
T_trn = 242
on = "months"
month19startC<-rolling(X_logC,X_linC,start = 1,T_trn,on)
month19blockC<-rolling(X_logC,X_linC,start = 0,T_trn,on)


par(mfrow=c(2,2))
chart.CumReturns(week19startC$result, main = "Contract 2019 performance of 1-week fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(week19blockC$result, main = "Contract 2019 performance of 1-week fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month19startC$result, main = "Contract 2019 performance of 1-month fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(month19blockC$result, main = "Contract 2019 performance of 1-month fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(week1819startC$result)),
                t(indices(week1819blockC$result)),
                t(indices(month1819startC$result)),
                t(indices(month1819blockC$result)),
                t(indices(week19startC$result)),
                t(indices(week19blockC$result)),
                t(indices(month19startC$result)),
                t(indices(month19blockC$result))),"indicesC.csv")

write.csv(week1819startC$weights,"week1819startC.csv")
write.csv(week1819blockC$weights,"week1819blockC.csv")
write.csv(month1819startC$weights,"month1819startC.csv")
write.csv(month1819blockC$weights,"month1819blockC.csv")
write.csv(week19startC$weights,"week19startC.csv")
write.csv(week19blockC$weights,"week19blockC.csv")
write.csv(month19startC$weights,"month19startC.csv")
write.csv(month19blockC$weights,"month19blockC.csv")


###########################################  rolling 2nd Fri ##########################################################

T_trn = 60
secFRI1819startC<-rolling2Fri(X_logC,X_linC,start = 1,T_trn)
secFRI1819blockC<-rolling2Fri(X_logC,X_linC,start = 0,T_trn)

par(mfrow=c(2,2))
chart.CumReturns(secFRI1819startC$result, main = "Contract 2rd Fri 18-19 fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(secFRI1819blockC$result, main = "Contract 2rd Fri 18-19 fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


###### 2019
X_log19tstC = X_logC[163:357,]
X_lin19tstC = X_linC[163:357,]

T_trn = 242
secFRI19startC<-rolling2Fri(X_logC,X_linC,start = 1,T_trn)
secFRI19blockC<-rolling2Fri(X_logC,X_linC,start = 0,T_trn)

chart.CumReturns(secFRI19startC$result, main = "Contract 2rd Fri 2019 fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(secFRI19blockC$result, main = "Contract 2rd Fri 2019 fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(secFRI1819startC$result)),
                t(indices(secFRI1819blockC$result)),
                t(indices(secFRI19startC$result)),
                t(indices(secFRI19blockC$result))),"indices2FriC.csv")

write.csv(secFRI1819startC$weights,"secFRI1819startC.csv")
write.csv(secFRI1819blockC$weights,"secFRI1819blockC.csv")
write.csv(secFRI19startC$weights,"secFRI19startC.csv")
write.csv(secFRI19blockC$weights,"secFRI19blockC.csv")

###############################################  rolling 1st Fri ###################################################

T_trn = 60
firstFRI1819startC<-rolling1Fri(X_logC,X_linC,start = 1,T_trn)
firstFRI1819blockC<-rolling1Fri(X_logC,X_linC,start = 0,T_trn)

par(mfrow=c(2,2))
chart.CumReturns(firstFRI1819startC$result, main = "Contract 1st Fri 18-19 fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(firstFRI1819blockC$result, main = "Contract 1st Fri 18-19 fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


###### 2019
X_log19tstC = X_logC[163:357,]
X_lin19tstC = X_linC[163:357,]

T_trn = 242
firstFRI19startC<-rolling1Fri(X_logC,X_linC,start = 1,T_trn)
firstFRI19blockC<-rolling1Fri(X_logC,X_linC,start = 0,T_trn)

chart.CumReturns(firstFRI19startC$result, main = "Contract 1st Fri 2019 fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(firstFRI19blockC$result, main = "Contract 1st Fri 2019 fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(firstFRI1819startC$result)),
                t(indices(firstFRI1819blockC$result)),
                t(indices(firstFRI19startC$result)),
                t(indices(firstFRI19blockC$result))),"indices1FriC.csv")

write.csv(firstFRI1819startC$weights,"firstFRI1819startC.csv")
write.csv(firstFRI1819blockC$weights,"firstFRI1819blockC.csv")
write.csv(firstFRI19startC$weights,"firstFRI19startC.csv")
write.csv(firstFRI19blockC$weights,"firstFRI19blockC.csv")

#########################################   rolling 3rd Fri  #######################################################
T_trn = 60
thirdFRI1819startC<-rolling3Fri(X_logC,X_linC,start = 1,T_trn)
thirdFRI1819blockC<-rolling3Fri(X_logC,X_linC,start = 0,T_trn)


par(mfrow=c(2,2))
chart.CumReturns(thirdFRI1819startC$result, main = "Contract 3rd Fri 18-19 fixed start",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(thirdFRI1819blockC$result, main = "Contract 3rd Fri 18-19 fixed blocks",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)


########### 2019
X_log19tstC = X_logC[163:357,]
X_lin19tstC = X_linC[163:357,]

T_trn = 242
thirdFRI19startC<-rolling3Fri(X_logC,X_linC,start = 1,T_trn)
thirdFRI19blockC<-rolling3Fri(X_logC,X_linC,start = 0,T_trn)

chart.CumReturns(thirdFRI19startC$result, main = "Contract 3rd Fri 2019 fixed start",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rainbow8equal)
chart.CumReturns(thirdFRI19blockC$result, main = "Contract 3rd Fri 2019 fixed blocks",
                 wealth.index = TRUE,  legend.loc = "topleft", colorset = rainbow8equal)



write.csv(rbind(t(indices(thirdFRI1819startC$result)),
                t(indices(thirdFRI1819blockC$result)),
                t(indices(thirdFRI19startC$result)),
                t(indices(thirdFRI19blockC$result))),"indices3FriC.csv")

write.csv(thirdFRI1819startC$weights,"thirdFRI1819startC.csv")
write.csv(thirdFRI1819blockC$weights,"thirdFRI1819blockC.csv")
write.csv(thirdFRI19startC$weights,"thirdFRI19startC.csv")
write.csv(thirdFRI19blockC$weights,"thirdFRI19blockC.csv")

