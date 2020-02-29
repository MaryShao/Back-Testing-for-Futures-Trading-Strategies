library(xts)
library(quantmod)
library(PerformanceAnalytics)
#install.packages('psych')
library('psych')
#install.packages('clusterSim')
library('clusterSim')
#install.packages('ade4')
library('ade4')


########### 分ratio
setwd("C:/Users/m8sha/Desktop/DATA")
perf <- read.csv("performance2018.txt", sep="", stringsAsFactors=FALSE)
#perf <- read.csv("performance.txt", sep="", stringsAsFactors=FALSE)

perf = perf[,!colnames(perf)%in%c("MA_ind30.02.LN_30", "v_LN_indv15_03_15","MA_ind120.07.LN_120","MA_SYL.MA.5min.index.1_5","MA_WS01.GT.15m60m_15","MA_WS07.GT.15m30m60m_60",
                                "ni_ind10.04.LN_10" ,"ni_LJ.Avanti02PP.5min_60","ni_LJ.MB02PP.5min_30" ,"ni_LJ.ninexiuPP.5min_60","ni_LJ.TW02PP.1H_15",
                                "OI_01.4H_240","OI_02.2H_120","OI_LJ.Kelther.Tsi01T.30min_240",
                                "OI_LJ.multsig03PP.15min_120","OI_LJ.multsig.A.4H_120","OI_LJ.multsig.A.4H_240","OI_LJ.threeswordPP.5min_30")]


perf.scale = perf
for (i in 1:nrow(perf)) {
  minimum = min(perf[i,])
  maximum = max(perf[i,])
  perf.scale[i,] = (perf[i,] - minimum)/(maximum - minimum)
}

perf.dist=dist(perf.scale,method ="euclidean",upper =T,diag=T)

sigma<-cor(as.matrix(perf.dist))
e<-eigen(sigma)
#特征值
#e$values
#特征向量
#e$vectors

# 确定主成分个数
#画碎石图
fa.parallel(as.matrix(sigma),n.obs=ncol(perf),fa='pc')
abline(h=1)

# 提取主成分
fit<-principal(sigma,
               nfactors=4,
               rotate='varimax',  # max variance
               scores=T)
fit

fa.diagram(fit,digits=2)

fit.cl<-principal(sigma,
                  nfactors=6,
                  rotate='cluster',  # max variance
                  scores=T)
fit.cl

fa.diagram(fit.cl,digits=2)


## RANK1234 PCA  ######################
# METHOD 1: 4 VARIMAX
ratio.sum = matrix(ncol=ncol(perf),nrow=4)
colnames(ratio.sum)=colnames(perf)

# factor 1: Burke Ratio
s1 =sort(as.numeric(perf[which(rownames(perf)=='Burke ratio (Risk free = 0)'),]), 
          decreasing= T, index.return = TRUE)
top20_burkeratio = perf[,s1$ix[1:20]]
c1=colnames(top20_burkeratio)
tmp=perf[,s1$ix]
for (i in 1:ncol(perf)){
  ratio.sum[1,which(colnames(ratio.sum)==colnames(tmp[i]))]=i}
  
# factor 2: Drawdown Deviation
s2 = sort(as.numeric(perf[which(rownames(perf)=='Drawdown Deviation'),]), 
          decreasing= F, index.return = TRUE)
top20_DDD = perf[,s2$ix[1:20]]
c2=colnames(top20_DDD)
tmp=perf[,s2$ix]
for (i in 1:ncol(perf)){
  ratio.sum[2,which(colnames(ratio.sum)==colnames(tmp[i]))]=i}

# factor 3: expectation
s3 = sort(as.numeric(perf[which(rownames(perf)=='expectation'),]), 
          decreasing= T, index.return = TRUE)
top20_exp = perf[,s3$ix[1:20]]
c3=colnames(top20_exp)
tmp=perf[,s3$ix]
for (i in 1:ncol(perf)){
  ratio.sum[3,which(colnames(ratio.sum)==colnames(tmp[i]))]=i}

# factor 4: cont loss days
s4=sort(as.numeric(perf[which(rownames(perf)=='continuous loss days'),]), 
        decreasing= F, index.return = TRUE)
top20_contlossday = perf[,s4$ix[1:20]]
c4=colnames(top20_contlossday)
tmp=perf[,s4$ix]
for (i in 1:ncol(perf)){
  ratio.sum[4,which(colnames(ratio.sum)==colnames(tmp[i]))]=i}
  
rownames(ratio.sum)=c('Burke ratio (Risk free = 0)','Drawdown Deviation','expectation','continuous loss days')


# factor 5: 25% factors
rank25 = matrix(nrow=1,ncol=ncol(perf))
for (i in 1:ncol(perf)){
  rank25[1,i]=sum(ratio.sum[,i])
}
colnames(rank25)=colnames(perf)

top20_rank25 = perf[,sort(as.numeric(rank25), 
                              decreasing= F, index.return = TRUE)$ix[1:20]]
c5=colnames(top20_rank25)  
  
  

# method 2:
ratio = matrix(ncol=ncol(perf),nrow=4)
colnames(ratio)=colnames(perf)

# factor 1: return
s5 =sort(as.numeric(perf[which(rownames(perf)=='Annualized Return'),]), 
         decreasing= T, index.return = TRUE)
tmp=perf[,s5$ix]
for (i in 1:ncol(perf)){
  ratio[1,which(colnames(ratio)==colnames(tmp[i]))]=i}

# factor 2: sharpe
s6 = sort(as.numeric(perf[which(rownames(perf)=='StdDev Sharpe (Rf=0%, p=95%):'),]), 
          decreasing= T, index.return = TRUE)
tmp=perf[,s6$ix]
for (i in 1:ncol(perf)){
  ratio[2,which(colnames(ratio)==colnames(tmp[i]))]=i}

# factor 3: MDD
s7 = sort(as.numeric(perf[which(rownames(perf)=='Worst Drawdown'),]), 
          decreasing= F, index.return = TRUE)
tmp=perf[,s7$ix]
for (i in 1:ncol(perf)){
  ratio[3,which(colnames(ratio)==colnames(tmp[i]))]=i} 

# factor 4: loss ratio
s8 = sort(as.numeric(perf[which(rownames(perf)=='gain ratio'),]), 
          decreasing= F, index.return = TRUE)
tmp=perf[,s8$ix]
for (i in 1:ncol(perf)){
  ratio[4,which(colnames(ratio)==colnames(tmp[i]))]=i}  

############# return50%+sharpe40%+MDD10% ##################################
rank541 = matrix(nrow=1,ncol=ncol(perf))
for (i in 1:ncol(perf)){
  rank541[1,i]=0.5*ratio541[1,i]+0.4*ratio541[2,i]+0.1*ratio541[2,i]
}
colnames(rank541)=colnames(perf)

top20_rank541 = perf[,sort(as.numeric(rank541), 
                          decreasing= F, index.return = TRUE)$ix[1:20]]
c6=colnames(top20_rank541)  

############# loss ratio 75%+sharpe15%+MDD10% ##################################
rank751 = matrix(nrow=1,ncol=ncol(perf))
for (i in 1:ncol(perf)){
  rank751[1,i]=0.1*ratio[3,i]+0.15*ratio541[2,i]+0.75*ratio[4,i]
}
colnames(rank751)=colnames(perf)

top20_rank751 = perf[,sort(as.numeric(rank751), 
                           decreasing= F, index.return = TRUE)$ix[1:20]]
c7=colnames(top20_rank751)  


top20_ratio96 = rbind(c1,c2,c3,c4,c5,c6)
colnames(top20_ratio96)=c(1:20)
rownames(top20_ratio96)=c('Burke ratio (Risk free = 0)','Drawdown Deviation','expectation',
                          'continuous loss days','25% of 4 factors',
                          '50%Return+40%Sharpe+10%MDD')

setwd("C:/Users/m8sha/Desktop")
write.csv(top20_ratio96,'top20_ratio78.csv')




#######################################################################
setwd("C:/Users/m8sha/Desktop/DATA/Strategy")
X_log <- read.csv("Strategy.Dailylogreturn.txt", sep="", stringsAsFactors=FALSE)
X_lin <- read.csv("Strategy.Dailylinreturn.txt", sep="", stringsAsFactors=FALSE)

X_log = X_log[,!colnames(X_log)%in%c("MA_ind30.02.LN_30", "v_LN_indv15_03_15","MA_ind120.07.LN_120","MA_SYL.MA.5min.index.1_5","MA_WS01.GT.15m60m_15","MA_WS07.GT.15m30m60m_60",
                                  "ni_ind10.04.LN_10" ,"ni_LJ.Avanti02PP.5min_60","ni_LJ.MB02PP.5min_30" ,"ni_LJ.ninexiuPP.5min_60","ni_LJ.TW02PP.1H_15","OI_01.4H_240","OI_02.2H_120","OI_LJ.Kelther.Tsi01T.30min_240",
                                  "OI_LJ.multsig03PP.15min_120","OI_LJ.multsig.A.4H_120","OI_LJ.multsig.A.4H_240","OI_LJ.threeswordPP.5min_30")]
X_lin = X_lin[,!colnames(X_lin)%in%c("MA_ind30.02.LN_30", "v_LN_indv15_03_15","MA_ind120.07.LN_120","MA_SYL.MA.5min.index.1_5","MA_WS01.GT.15m60m_15","MA_WS07.GT.15m30m60m_60",
                                  "ni_ind10.04.LN_10" ,"ni_LJ.Avanti02PP.5min_60","ni_LJ.MB02PP.5min_30" ,"ni_LJ.ninexiuPP.5min_60","ni_LJ.TW02PP.1H_15","OI_01.4H_240","OI_02.2H_120","OI_LJ.Kelther.Tsi01T.30min_240",
                                  "OI_LJ.multsig03PP.15min_120","OI_LJ.multsig.A.4H_120","OI_LJ.multsig.A.4H_240","OI_LJ.threeswordPP.5min_30")]
rownames(X_log)=rownames(X_lin)=X_log[,1]
X_log18=X_log[1:which(rownames(X_log)=='2018-12-28 15:15:00'),2:ncol(X_log)]
X_lin19=X_log[which(rownames(X_lin)=='2019-01-02 15:15:00'):nrow(X_lin),2:ncol(X_lin)]






ReturnCalculate <- function(X_log,X_lin){
  w_all <- CalculateWeights(X_log,N=ncol(X_log))
  result <- as.matrix(X_lin)%*%w_all
  
  rownames(result) = rownames(as.matrix(X_lin))
  colnames(result) = colnames(w_all)
  
  return(list('result'=result,
              'weight'=w_all))
}


# METHOD 1: 4 VARIMAX
# factor 1: Burke Ratio
X_log1=X_log18[,s1$ix[1:20]]
X_lin1=X_lin19[,s1$ix[1:20]]
rt1 = ReturnCalculate(X_log1,X_lin1)
c_rt1 = chart.CumReturns(rt1$result, main = "Using Factor1 Burke Ratio top20",
                         wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)

# factor 2: Drawdown Deviation
X_log2=X_log18[,s2$ix[1:20]]
X_lin2=X_lin19[,s2$ix[1:20]]
rt2 = ReturnCalculate(X_log2,X_lin2)
c_rt2 = chart.CumReturns(rt2$result, main = "Using Factor2 Drawdown Deviation top20",
                         wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)

# factor 3: expectation
X_log3=X_log18[,s3$ix[1:20]]
X_lin3=X_lin19[,s3$ix[1:20]]
rt3 = ReturnCalculate(X_log3,X_lin3)
c_rt3 = chart.CumReturns(rt3$result, main = "Using Factor3 expectation top20",
                         wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)


# factor 4: cont loss days
X_log4=X_log18[,s4$ix[1:20]]
X_lin4=X_lin19[,s4$ix[1:20]]
rt4 = ReturnCalculate(X_log4,X_lin4)
c_rt4 = chart.CumReturns(rt4$result, main = "Using Factor4 cont loss days top20",
                         wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)


# 25% factors
# X_log5=X_log18[,sort(as.numeric(rank25), 
#                      decreasing= F, index.return = TRUE)$ix[1:20]]
# X_lin5=X_lin19[,sort(as.numeric(rank25), 
#                      decreasing= F, index.return = TRUE)$ix[1:20]]
# rt5 = ReturnCalculate(X_log5,X_lin5)
rt_25 = (rt1$result+rt2$result+rt3$result+rt4$result)/4
c_rt25 = chart.CumReturns(rt_25, main = "Using 25% factors top20",
                         wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)


# factor 5: return
X_log5=X_log18[,s5$ix[1:20]]
X_lin5=X_lin19[,s5$ix[1:20]]
rt5 = ReturnCalculate(X_log5,X_lin5)

# factor 6: sharpe
X_log6=X_log18[,s6$ix[1:20]]
X_lin6=X_lin19[,s6$ix[1:20]]
rt6 = ReturnCalculate(X_log6,X_lin6)

# factor 7: MDD
X_log7=X_log18[,s7$ix[1:20]]
X_lin7=X_lin19[,s7$ix[1:20]]
rt7 = ReturnCalculate(X_log7,X_lin7)

# factor 8: loss ratio
X_log8=X_log18[,s8$ix[1:20]]
X_lin8=X_lin19[,s8$ix[1:20]]
rt8 = ReturnCalculate(X_log8,X_lin8) 

##return50%+sharpe40%+MDD10%
rt_541 = 0.5*rt5$result+0.4*rt6$result+0.1*rt7$result
c_rt541 = chart.CumReturns(rt_541, main = "Using return50%+sharpe40%+MDD10% top20",
                          wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)

##loss ratio 75%+sharpe15%+MDD10%
rt_751 = 0.75*rt8$result+0.15*rt6$result+0.1*rt7$result
c_rt751 = chart.CumReturns(rt_751, main = "Using loss ratio 75%+sharpe15%+MDD10% top20",
                           wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)









# # method 2:
# 
# #return50%+sharpe40%+MDD10%
# 
# X_log6=X_log18[,sort(as.numeric(rank541),
#                      decreasing= F, index.return = TRUE)$ix[1:20]]
# X_lin6=X_lin19[,sort(as.numeric(rank541),
#                      decreasing= F, index.return = TRUE)$ix[1:20]]
# rt6 = ReturnCalculate(X_log6,X_lin6)
# c_rt6 = chart.CumReturns(rt6$result, main = "Using return50%+sharpe40%+MDD10% top20",
#                          wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)
# 
# # loss ratio 75%+sharpe15%+MDD10%
# X_log7=X_log18[,sort(as.numeric(rank751),
#                      decreasing= F, index.return = TRUE)$ix[1:20]]
# X_lin7=X_lin19[,sort(as.numeric(rank751),
#                      decreasing= F, index.return = TRUE)$ix[1:20]]
# rt7 = ReturnCalculate(X_log7,X_lin7)
# c_rt7 = chart.CumReturns(rt7$result, main = "Using loss ratio 75%+sharpe15%+MDD10% top20",
#                          wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)


## ALL
rt_all = ReturnCalculate(X_log18,X_lin19)

# par(mfrow=c(3,2))
# c_rt1
# c_rt2
# c_rt3
# c_rt4
# c_rt5
# c_rt6

par(mfrow=c(1,3))
c_rt25
c_rt541
c_rt751

layout(matrix(1,1,1))


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

ind1 = indices(rt1$result)
ind2= indices(rt2$result)
ind3= indices(rt3$result)
ind4 = indices(rt4$result)
ind25 = indices(rt_25)
ind541 = indices(rt_541)
ind751 = indices(rt_751)
indall= indices(rt_all$result)

setwd("C:/Users/m8sha/Desktop")
write.csv(t(ind25),'25%(BurkeRatio+DDD+exp+contloss).csv')
write.csv(t(ind541),'return50%+sharpe40%+MDD10%.csv')
write.csv(t(ind751),"lossratio75%+sharpe15%+MDD10%.csv")
write.csv(t(indall),"all.csv")
write.csv(t(ind1),"Burke Ratio.csv")
write.csv(t(ind2),"Drawdown Deviation.csv")
write.csv(t(ind3),"expectation.csv")
write.csv(t(ind4),"contlossday.csv")




rolling <- function(X_log,X_lin,on,start,T_trn){
  T_trn <- which(rownames(X_log)=='2018-12-28 15:15:00')
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
    tst <- X_log[(rebal_indices[i]+1):min(rebal_indices[i+1],nrow(X_log)),mdd_top5] 
    # calculate mean returns as result for each days
    result <- c(result,apply(matrix(tst,ncol=5), 1, mean))
  }
  
  names(result) = rownames(X_log)[(T_trn+1):t]
  
  return(list("sharpe"=SharpeRatio(result, Rf = 0, p = 0.95, FUN = "StdDev"),
              "MDD"=maxDrawdown(result, weights = NULL, geometric = TRUE, invert = TRUE),
              "result" = result))
  
}












