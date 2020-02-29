########################################## PP 31 INDEX && DAILY/MONTHLY DATA ##########################################
# daily rt
pp[,1] = as.character(pp[,1])
day = unique(substr(pp[,1],1,10))


pp.day = matrix(0,346,11)
for(k in 1:length(day)){
    pp.day[k,] = apply(pp[substr(pp[,1],1,10)==day[k],2:12],2,sum)
}

rownames(pp.day) = day
colnames(pp.day) = colnames(pp)[-1]
pp.day=pp.day[-346,]

#monthly rt
pp.month = matrix(0,17,11)
m = unique(substr(pp[,1],1,7))[-18]
rownames(pp.month) = m
colnames(pp.month) = colnames(pp)[-1]

for(k in 1:length(m)){
  pp.month[k,] = apply(pp[substr(pp[,1],1,7)==m[k],2:12],2,sum)
}


####连续亏损
c = numeric(length = nrow(pp.month))
tmp = matrix(0,17,11)
for(j in 1:11){
  for(i in 1:nrow(pp.month)){
    if(pp.month[i,j]<0){
      tmp[i,j] = -1
    }
  }
}

for(j in 1:11){
  for(i in 2:nrow(pp.month)){
    if(tmp[i,j]<0 && tmp[i-1,j]<0){
      tmp[i,j] = tmp[i-1,j]-1
    }
  }
}

c = abs(apply(tmp,2,min))

###### 交易次数 胜率 期望
totalnumber = numeric(length = ncol(allSig))
rtByTrade = matrix(0,500,ncol(allSig))

for(j in 1:ncol(allSig)){
  for (i in 1:nrow(allSig)) {
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
for(j in 1:ncol(rtByTrade)){
  gainNum = sum(rtByTrade[,j]>0)
  lossNum = sum(rtByTrade[,j]<0)
  totalGain = sum(rtByTrade[rtByTrade[,j]>0,j])
  totalLoss = sum(rtByTrade[rtByTrade[,j]>0,j])
  gain.ratio[j] = gainNum/totalnumber[j]
  expe[j] = (gainNum*totalGain + lossNum*totalLoss)/totalnumber[j]
}


########################BENCHMARK###############################
benchmark = matrix(0,nrow = dim(pp.day)[1], ncol = dim(pp.day)[2])
colnames(benchmark) =colnames(pp)[-1]
rownames(benchmark) = rownames(pp.day)

for (i in 1:11) {
  benchmark[,i] = rep(0.035,345)
}
########################BENCHMARK###############################

benchmark<-as.xts(benchmark)

perf= matrix(nrow = 0,ncol = 11) 
colnames(perf) = colnames(pp)[-1]

# calculate performance###################################################################

inf_ratio <- InformationRatio(pp.day, benchmark, scale = 252)
inf_ratio = c(inf_ratio[2,1],inf_ratio[1,2:11])

MAR=0.05

perf= rbind(perf,
            colMeans(pp.day),
            diag(cov(pp.day)),
            table.AnnualizedReturns(pp.day),
            BernardoLedoitRatio(pp.day),
            BurkeRatio(pp.day, Rf = 0, modified = FALSE),
            CalmarRatio(pp.day, scale = NA),
            #CDD(pp.day, weights = NULL, geometric = TRUE, invert = TRUE, p = 0.95),
            DownsideDeviation(pp.day, MAR = 0, method = c("full", "subset"),potential = FALSE),
            DRatio(pp.day),
            DrawdownDeviation(pp.day),
            inf_ratio,
            KellyRatio(pp.day, Rf = 0, method = "half"),
            MartinRatio(pp.day, Rf = 0),
            maxDrawdown(pp.day, weights = NULL, geometric = TRUE, invert = TRUE),
            MeanAbsoluteDeviation(pp.day),
            Omega(pp.day, L = 0, method = "simple", output = "point", Rf = 0),
            #PainRatio(pp.day, Rf = 0),
            ProspectRatio(pp.day,MAR),
            SharpeRatio(pp.day, Rf = 0, p = 0.95, FUN = "StdDev"),
            SkewnessKurtosisRatio(pp.day),
            SmoothingIndex(pp.day, neg.thetas = FALSE, MAorder = 2, verbose = FALSE),
            SortinoRatio(pp.day, MAR),
            UlcerIndex(pp.day),
            UpsideFrequency(pp.day, MAR = 0),
            VaR(pp.day, p = 0.95,method = "gaussian", portfolio_method ="single",invert = TRUE),
            c)

rownames(perf)[1:2] = c('mu','cov')
rownames(perf)[12] = c('Information ratio')
rownames(perf)[26] = c('continuous loss months')

setwd("~/Tencent Files/1060573664/FileRecv")
append.perf = as.matrix(read_excel("PP.xlsx")[,2:12])
colnames(append.perf) = colnames(perf)
perf = rbind(perf,append.perf)
rownames(perf)[27:29] = c("total numbers", "gain ratio", "expectation")

################################################ ALL STRATEGIES ######################################################

all = Strategy.signals[,3:84]


# daily rt
all[,1] = as.character(all[,1])
date_day = unique(substr(all[,1],1,10))


all.day = matrix(0,length(date_day),n_strategy)
for(k in 1:length(date_day)){
  all.day[k,] = apply(all[substr(all[,1],1,10)==date_day[k],2:ncol(all)],2,sum)
}

rownames(all.day) = date_day
colnames(all.day) = colnames(pp)[-1]
pp.day=pp.day[-346,]

#monthly rt
pp.month = matrix(0,17,11)
m = unique(substr(pp[,1],1,7))[-18]
rownames(pp.month) = m
colnames(pp.month) = colnames(pp)[-1]

for(k in 1:length(m)){
  pp.month[k,] = apply(pp[substr(pp[,1],1,7)==m[k],2:12],2,sum)
}


############ all strategy signal cluster #############################
allSig = Strategy.signals[,3:84]
rownames(allSig) = paste(Strategy.signals[,1],Strategy.signals[,2])

allSig.Euc = dist(t(allSig))
plot(hclust(allSig.Euc,"ward"))

# hour
allSig.hour = allSig[substr(rownames(allSig), start = 15, stop = 19) == "00:00",]

# dist
allSigHour.dtw = dist(t(allSig.hour),method = "dtw")
plot(hclust(allSigHour.dtw,"complete"))


allSigHour.Euc = dist(t(allSig.hour))
plot(hclust(allSigHour.Euc,"ward"))


# JSD
all.hourProb = matrix(0,nrow=3,ncol = ncol(allSig.hour))
colnames(all.hourProb) = colnames(allSig.hour)
rownames(all.hourProb) = c(0,-1,1)
for(i in 1:ncol(allSig.hour)){
  for (k in 1:3) {
    all.hourProb[k,i] = sum(allSig.hour[,i] == rownames(all.hourProb)[k])/nrow(allSig.hour)
  }
}

allSigHour.JSD = dist.JSD(all.hourProb)
plot(hclust(allSigHour.JSD,"ward"))

rt.pam <- pam(allSigHour.Euc, k=12, diss = TRUE)
plot(rt.pam)

# day
allSig.day = allSig[substr(rownames(allSig), start = 12, stop = 19) == "15:15:00",]
allSigDay.Euc = dist(t(allSig.day))
plot(hclust(allSigDay.Euc,"ward"))



