########################################
ord18<-c()


for (i in 1:nrow(perf2018)){
  if (row.names(perf2018)[i]=="cov"||
      row.names(perf2018)[i]=="Annualized Std Dev"||
      row.names(perf2018)[i]=="Downside Deviation (MAR = 0%)"||
      row.names(perf2018)[i]=="d ratio"||
      row.names(perf2018)[i]=="Drawdown Deviation"||
      row.names(perf2018)[i]=="Worst Downside"||
      row.names(perf2018)[i]=="Mean absolute deviation"||
      row.names(perf2018)[i]=="Ulcer Index"||
      row.names(perf2018)[i]=="continuous loss days"){
    ord18 = rbind(ord18,sort(as.numeric(perf2018[i,]), decreasing= F, index.return = TRUE)$ix)
  }
  else{
    ord18 = rbind(ord18,sort(as.numeric(perf2018[i,]), decreasing = TRUE, index.return = TRUE)$ix)
  }
}

rownames(ord18)=rownames(perf2018)
colnames(ord18)=colnames(perf2018)


# w = matrix(0,nrow(perf2018),ncol(perf2018))
# rownames(w)=rownames(perf2018)
# colnames(w)=colnames(perf2018)
# w = t(w)
# 
# for(j in 1:ncol(w)){
#   w[ord[j,1],j] <- 1/round(N/4)
# }
# 
# rt_all = xts(2, rownames(rt.day))


top20<- apply(ord18[c("Sortino Ratio (MAR = 0%)","Drawdown Deviation","Upside Frequency (MAR = 0%)","gain ratio","Annualized Return","StdDev Sharpe (Rf=0%, p=95%):","Worst Drawdown"),1:20], 1, function(x){colnames(ord18)[x]})


# #fourComp <- apply(ord18[c("Sortino Ratio (MAR = 0%)","Drawdown Deviation","Upside Frequency (MAR = 0%)","gain ratio"),], 1, function(x){colnames(ord18)[x]})
# score18<-matrix(0,7,78)
# rownames(score18) <-c("Sortino Ratio (MAR = 0%)","Drawdown Deviation","Upside Frequency (MAR = 0%)","gain ratio","Annualized Return","StdDev Sharpe (Rf=0%, p=95%):","Worst Drawdown")
# colnames(score18) <- colnames(ord18)
# 
# for (i in 1:nrow(score18)) {
#   score18[i,ord18[rownames(score18)[i],]] = 78:1
# }
# 
# top20=cbind(top20,names(sort(apply(score18[1:4,],2,mean), decreasing= TRUE, index.return = TRUE)$x)[1:20],
#       names(sort(apply(score18[5:7,],2,function(x){0.5*x[1]+0.4*x[2]+0.1*x[3]}), decreasing= TRUE, index.return = TRUE)$x)[1:20])


write.csv(t(top20),"top20.csv")



#18 train 19 test
i=4
X_log18 <- X_log[1:242,]
X_lin19 <- X_lin[-(1:242),]
w_all <- wght(X_log18)
# get test return
result <- as.matrix(X_lin19)%*%w_all
chart.CumReturns(rank2, main = "18-19 performance",
                 wealth.index = TRUE,legend.loc = "topleft", colorset = rainbow8equal,geometric = FALSE)


avg <- (result1+result2+result3+result4)/4
rank2 <- 0.5*result5 + 0.4*result6 + 0.1*result7
rank3 <- 0.15*result6 + 0.1*result7 + 0.75*result4

write.csv(rbind(t(indices(result)),
                t(indices(result1)),
                t(indices(result2)),
                t(indices(result3)),
                t(indices(result4)),
                t(indices(avg)),
                t(indices(rank2)),
                t(indices(rank3))),"indicesTop20.csv")

unique(as.vector(top20[,c(1:4)]))
unique(as.vector(top20[,c(5:7)]))
unique(as.vector(top20[,c(4,6,7)]))



rolling20 <- function(X_log,X_lin,start,T_trn,on){
  t <- nrow(X_log)
  X_log_trn <- X_log[1:T_trn, ]
  X_log_tst <- X_log[(T_trn+1):t, ]
  rebal_indices <- T_trn + endpoints(X_log_tst, on = on)
  lookback <- T_trn
  result =matrix(ncol=7)[-1,]
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
###################################################################################################################

singleStr = c()

for(j in 1:nrow(ord18)){
  singleStr =c(singleStr,colnames(perf2018)[ord18[j,1]])
}
singleStr_indx = paste(singleStr,rownames(perf2018))


split = 243
rt_all =  rt.day[,singleStr]
rt_all_trn = rt_all[1:split,]
rt_all_tst = rt_all[-(1:split),]
colnames(rt_all_tst) = singleStr_indx

# eva = matrix(0,nrow(rt_all)-1,ncol(rt_all))
# 
# for(i in 2:nrow(rt_all)){
#   eva[i-1,] = CalmarRatio(rt_all[1:i,], scale = NA)
# }
top8 = names(sort(apply(rt_all_tst, 2, sum),decreasing = TRUE, index.return = TRUE)$x[1:10])

rt19pos = rt_all_tst[,top8]

delDrawdown = rt19pos[,maxDrawdown(rt19pos, weights = NULL, geometric = TRUE, invert = TRUE)<0.05]



chart.CumReturns(rt19pos, main = "Performance of different portfolios",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)



############## combine top 5 of one index

fiveStr = matrix(0,nrow(perf2018),5)
rownames(fiveStr) = rownames(perf2018)

for(j in 1:nrow(ord18)){
  fiveStr[j,] =colnames(perf2018)[ord18[j,1:5]]
  #fiveStr_indx = paste(fiveStr,rownames(perf2018))
}

i=1


split = 243
rt_fiveStr = matrix(0, nrow = nrow(rt.day),ncol = nrow(perf))
rownames(rt_fiveStr) = rownames(rt.day)
colnames(rt_fiveStr) = rownames(perf)

for(j in 1:ncol(rt_fiveStr)){
  rt_fiveStr[,j] =  apply(0.2*rt.day[,fiveStr[j,]],1,sum)
}

rt_five_trn = rt_fiveStr[1:split,]
rt_five_tst = rt_fiveStr[-(1:split),]
#colnames(rt_five_tst) = fiveStr[i,]

# eva = matrix(0,nrow(rt_all)-1,ncol(rt_all))
# 
# for(i in 2:nrow(rt_all)){
#   eva[i-1,] = CalmarRatio(rt_all[1:i,], scale = NA)
# }
top10 = names(sort(apply(rt_five_tst, 2, sum),decreasing = TRUE, index.return = TRUE)$x[1:10])

rt_five_19pos = rt_five_tst[,top10]

delDrawdown_five = rt_five_19pos[,maxDrawdown(rt_five_19pos, weights = NULL, geometric = TRUE, invert = TRUE)<0.05]



chart.CumReturns(delDrawdown_five, main = "Performance of different portfolios",
                 wealth.index = TRUE, legend.loc = "topleft", colorset = rich8equal)

