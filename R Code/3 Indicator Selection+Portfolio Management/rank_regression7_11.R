
library(xts)
library(quantmod)
library(PerformanceAnalytics)
#install.packages('psych')
library('psych')
#install.packages('clusterSim')
library('clusterSim')
#install.packages('ade4')
library('ade4')

####################### DATA #############################
# allSig
Strategy.signals1 <- read.csv("C:/Users/m8sha/Desktop/updated data/Strategy.signals1.txt", row.names=NULL, sep="")
Signals.zn <- read.csv("C:/Users/m8sha/Desktop/updated data/Signals.zn.txt", sep="", stringsAsFactors=FALSE)
name1 = Strategy.signals1[,1]
names2 = Strategy.signals1[,2]
rownames(Strategy.signals1)=paste(name1,names2,sep=' ')
allSig = Strategy.signals1[2:nrow(Strategy.signals1),3:ncol(Strategy.signals1)]
#sig.zn=Signals.zn[3:ncol(Signals.zn)]
#allSig = as.xts(cbind(Strategy.signals1[,3:ncol(Strategy.signals1)],sig.zn))

# rt
Strategy.logcr <- read.csv("C:/Users/m8sha/Desktop/updated data/Strategy.logcr.txt", sep="", stringsAsFactors=FALSE)
rt = Strategy.logcr[,1:ncol(allSig)]
rownames(rt)=rownames(allSig)
#rt.day
Strategy.Dailylogreturn <- read.csv("C:/Users/m8sha/Desktop/updated data/Strategy.Dailylogreturn.txt", sep="", stringsAsFactors=FALSE)
rt.day =Strategy.Dailylogreturn[,2:ncol(Strategy.Dailylogreturn)]
rownames(rt.day)=substr(Strategy.Dailylogreturn[,1],1,10)

#allSig = as.xts(allSig)
rt = as.xts(rt)
rt.day = as.xts(rt.day)
allSig18 = allSig[1:which(rownames(allSig)=='2018-12-28 15:15:00'),]
rt18 = rt[1:which(index(rt)=='2018-12-28 15:15:00 CST'),]
rt.day18 = rt.day[1:which(index(rt.day)=='2018-12-28 CST'),]
########################BENCHMARK###############################

# calculate performance############################################
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
    
    ## 交易次数 胜率 期望
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

perf18=performance(allSig18,rt18,rt.day18)
# setwd("C:/Users/m8sha/Desktop")
# write.table(perf18,'perf18.txt')

###############Adjusted perf################

scale<-function(R){
  R.scale = R
  for (i in 1:nrow(R)) {
    minimum = min(R[i,])
    maximum = max(R[i,])
    R.scale[i,] = (R[i,] - minimum)/(maximum - minimum)
  }
  return(as.matrix(R))
}


check.perf = cor(t(perf18))
perf.adj = perf18[-c(1,2,6,7,8,11,12,14,19,20,9,16,22),]

########## Rank #######################################################
CalculateRank <- function(perf){
  ratio = matrix(ncol=ncol(perf),nrow=nrow(perf))
  colnames(ratio)=colnames(perf)
  rownames(ratio)=rownames(perf)
  for (i in 1: nrow(perf)){
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
    for (j in 1:ncol(perf)){
      ratio[i,which(colnames(ratio)==colnames(tmp[j]))]=j}
  }
  return(ratio)
}

rank.perf=CalculateRank(perf18)
# setwd("C:/Users/m8sha/Desktop")
# write.table(rank.perf,'rankperf18.txt')
rank.perfadj = CalculateRank(perf.adj)


############ Regression#########################
#install.packages('CRAN')
# install.packages(Rfit)
library(Rfit)
library(tseries)
library(car)

## linear regression #############################

### after cor selection ####################

calmer =  t(perf18[8,])
dat = as.data.frame(t(rank.perfadj))


#  (lm1: rank_annul.rt~.,data=rank)
lm1 = lm(rank.perfadj[1,]~.,data = as.data.frame(t(rank.perfadj[2:nrow(rank.perfadj),])))
summary(lm1)
vif(lm1, digits = 3)

lm1.step <- step(lm1, direction = "backward")
summary(lm1.step)

#  (lm2: scale.mu~rank)
scaled.rt =  t(scale(perf18[3,]))
scaled.dat = as.data.frame(cbind(1000*scaled.rt,t(rank.perfadj)))
lm2 = lm(scaled.dat[,1]~.,data = as.data.frame(scaled.dat[,2:ncol(scaled.dat)]))
vif(lm2, digits = 3)

lm2.step <- step(lm2, direction = "backward")
summary(lm2.step)

#  (lm3: calmer~rank)

lm3 = lm(calmer~.,data=dat)
summary(lm3)
vif(lm3, digits = 3)

lm3.step <- step(lm3, direction = "backward")
summary(lm3.step)


## rank-based regression #
#(rm1: rank_annul.rt~.,data=rank)
rm1 = rfit(rank.perfadj[1,]~.,data = as.data.frame(t(rank.perfadj[2:nrow(rank.perfadj),])))
summary(rm1)

#(rm2: calmer~.,data=rank)
rm2 = rfit(calmer~.,data=dat)
summary(rm2)

################################ before cor selection ###############

#  (lm4: calmer~rankall)
dat = as.data.frame(t(rank.perf[-c(8),]))
lm4 = lm(calmer~.,data=dat)
summary(lm4)
vif(lm4, digits = 3)

lm4.step <- step(lm4, direction = "backward")
summary(lm4.step)


## rank-based regression *(rm1: rank_annul.rt~.,data=rank)


# rm3 = rfit(calmer~.,data=dat)
# summary(rm3)

summary(lm4.step)
summary(rm2)
summary(lm3.step)


# top5_sharpe = perf[,sort(as.numeric(perf[which(rownames(perf)=='Annualized Sharpe (Rf=0%)'),]),
#                                          decreasing= T, index.return = TRUE)$ix[1:5]]
# top5_MDD = perf[,sort(as.numeric(perf[which(rownames(perf)=='Worst Drawdown'),]),
#                          decreasing= F, index.return = TRUE)$ix[1:5]]
# top5_prospect = perf[,sort(as.numeric(perf[which(rownames(perf)=='Prospect ratio (MAR = 0%)'),]),
#                          decreasing= T, index.return = TRUE)$ix[1:5]]
# top5_upsidefq = perf[,sort(as.numeric(perf[which(rownames(perf)=='Upside Frequency (MAR = 0%)'),]),
#                          decreasing= T, index.return = TRUE)$ix[1:5]]

rank25 = t(as.matrix(apply(rank.perf[c(5,15,18,24),],2,sum)))
colnames(rank25)=colnames(perf)

top20_rank25 = perf[,sort(as.numeric(rank25), 
                          decreasing= F, index.return = TRUE)$ix[1:20]]
top20_rt= rt.day[,sort(as.numeric(rank25), 
                     decreasing= F, index.return = TRUE)$ix[1:20]]

chart.CumReturns(top20_rt[which(index(top20_rt)=='2018-12-28 CST'):nrow(top20_rt),],
                 main = "top20 during 2019",wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)

top10_rt= rt.day[,sort(as.numeric(rank25), 
                       decreasing= F, index.return = TRUE)$ix[1:10]]

chart.CumReturns(top10_rt[which(index(top20_rt)=='2018-12-28 CST'):nrow(top10_rt),],
                 main = "top10 during 2019",wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow12equal)
