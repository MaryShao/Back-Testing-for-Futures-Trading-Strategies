library(xts)
library(quantmod)
library(PerformanceAnalytics)
#install.packages('psych')
library('psych')
#install.packages('clusterSim')
library('clusterSim')
#install.packages('ade4')
library('ade4')

######## preparing function ########################

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

################ data import #######################
setwd("C:/Users/m8sha/Desktop/updated data")
########## Daily
# rt.day
X_lin <- read.csv('Strategy.Dailylinreturn.txt',header = TRUE,sep = '')
X_log <- read.csv('Strategy.Dailylogreturn.txt',header = TRUE,sep = '')
rownames(X_lin) <- substr(X_lin[,1],1,10)
rownames(X_log) <- substr(X_log[,1],1,10)
X_lin <- X_lin[,-1]
X_log <- X_log[,-1]
X_lin <- as.xts(X_lin)
X_log <- as.xts(X_log)

############## By minute
# allSig
Strategy.signals1 <- read.csv("Strategy.signals1.txt", row.names=NULL, sep="")
name1 = Strategy.signals1[,1]
names2 = Strategy.signals1[,2]
rownames(Strategy.signals1)=paste(name1,names2,sep=' ')
allSig = Strategy.signals1[2:nrow(Strategy.signals1),3:ncol(Strategy.signals1)]

# rt
Strategy.logcr <- read.csv("Strategy.logcr.txt", sep="", stringsAsFactors=FALSE)
rt = Strategy.logcr[,1:ncol(allSig)]
rownames(rt)=rownames(allSig)
rt = as.xts(rt)


################ rolling ###########################

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


### block
blockunif = rollingUnif(X_log,X_lin,allSig,rt,start=0)
setwd("C:/Users/m8sha/Desktop")
# write.csv(blockunif$performances,'performances.csv')
# write.csv(blockunif$rank.perf,'rank.perf.csv')
write.csv(blockunif$strategy[[1]],'strategyblock1.csv')
write.csv(blockunif$strategy[[2]],'strategyblock2.csv')
write.csv(blockunif$strategy[[3]],'strategyblock3.csv')
write.csv(blockunif$strategy[[4]],'strategyblock4.csv')
write.csv(blockunif$strategy[[5]],'strategyblock5.csv')
write.csv(blockunif$strategy[[6]],'strategyblock6.csv')
chart.CumReturns(blockunif$result, main = "18block unif",wealth.index = T, colorset = rainbow12equal)

ind_block = indices(blockunif$result)
write.csv(ind_block,'ind_block.csv')

### fix
fixunif = rollingUnif(X_log,X_lin,allSig,rt,start=1)
setwd("C:/Users/m8sha/Desktop")
# write.csv(fixunif$performances,'performances.csv')
# write.csv(fixunif$rank.perf,'rank.perf.csv')
write.csv(fixunif$strategy[[1]],'strategyfix1.csv')
write.csv(fixunif$strategy[[2]],'strategyfix2.csv')
write.csv(fixunif$strategy[[3]],'strategyfix3.csv')
write.csv(fixunif$strategy[[4]],'strategyfix4.csv')
write.csv(fixunif$strategy[[5]],'strategyfix5.csv')
write.csv(fixunif$strategy[[6]],'strategyfix6.csv')
chart.CumReturns(fixunif$result, main = "18fix unif",wealth.index = T, colorset = rainbow12equal)

ind_fix = indices(fixunif$result)
write.csv(ind_fix,'ind_fix.csv')

