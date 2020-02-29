
# Data import

# 每手
contract<- read.csv("C:/Users/m8sha/Desktop/contract_number.csv", stringsAsFactors=FALSE)
rownames(contract)=contract[,1]
contract=contract[,2:3]
#coef=contract[,-1]


# 合约价格
setwd("C:/Users/m8sha/Desktop/DATA/Sig_Price")
f = list.files()

tmp=read.table(f[1],header = TRUE)
tmp = tmp[which(substr(tmp[,1],12,19)=="15:00:00"),]
rownames(tmp)= substr(tmp[,1],1,11)
tmp = tmp[,2:ncol(tmp)]
tmp = tmp[c(1,endpoints(as.xts(tmp),on='month')),]
price = as.matrix(tmp[,1])
rownames(price) = rownames(tmp)

amount = c((ncol(tmp)-1))

for(i in f[-1]){
  tmp = read.table(i,header = TRUE)
  tmp[,1] = as.character(tmp[,1])
  tmp = tmp[which(substr(tmp[,1],12,19)=="15:00:00"),]
  rownames(tmp)= substr(tmp[,1],1,11)
  tmp = tmp[,2:ncol(tmp)]
  tmp = tmp[c(1,endpoints(as.xts(tmp),on='month')),]
  price = cbind(price,tmp[,1])
  amount = c(amount,(ncol(tmp)-1))
}

amount = as.matrix(amount)
names = c('a','hc','i','IF','j','PP','rb','SR','T','ZC')
colnames(price)=rownames(amount)=names

###### 
setwd("C:/Users/m8sha/Desktop/DATA/Strategy")
Strategy.returns <- read.csv("Strategy.returns.txt", 
                             row.names=NULL, sep="", stringsAsFactors=FALSE)
allRt = Strategy.returns[,3:98]
rownames(allRt) = paste(Strategy.returns[,1],Strategy.returns[,2])

total.contract = price

coef.T = coef[which(coef[,1]=='T'),]
amount.T = as.numeric(amount[which(rownames(amount)=='T'),])
mark_p_T = price[,which(colnames(price)=='T')]*as.numeric(coef.T[,2])*amount.T

for (i in names){
  if (i=='T'){total.contract[,which(colnames(total.contract)==i)]=1}
  else{
    coef.i= coef[which(coef[,1]==i),]
    amount.i = as.numeric(amount[which(rownames(amount)==i),])
    num = mark_p_IF/(price[,which(colnames(price)==i)]*as.numeric(coef.i[,2])*amount.i)
    total.contract[,which(colnames(total.contract)==i)]=round(num)
  }
}

# setwd("C:/Users/m8sha/Desktop")
# write.csv(total.contract,'total.contract.csv')

############# return #################
Strategy.Dailylinreturn <- read.csv("C:/Users/m8sha/Desktop/DATA/Strategy/Strategy.Dailylinreturn.txt",
                                    row.names=1, sep="")
rt = matrix(ncol=ncol(Strategy.Dailylinreturn),nrow=nrow(Strategy.Dailylinreturn))
rownames(rt)=substr(rownames(Strategy.Dailylinreturn),1,10)
colnames(rt)=colnames(Strategy.Dailylinreturn)

ind1 = which(substr(rownames(allRt),12,19)=='09:01:00')
ind2 = which(substr(rownames(allRt),12,19)=='15:15:00')
  
for (i in 1:ncol(Strategy.Dailylinreturn)){
  for (j in 1:nrow(Strategy.Dailylinreturn)){
    if (j ==1){day.rt = allRt[ind1[j]:ind2[j],i]}
    else{day.rt = c(allRt[ind1[j]:ind2[j],i],allRt[ind2[j-1]:ind1[j],i])}
    rt[j,i]=sum(day.rt)
  }
}

rt = rt[,!colnames(rt)%in%c("MA_ind30.02.LN_30", "v_LN_indv15_03_15","MA_ind120.07.LN_120","MA_SYL.MA.5min.index.1_5","MA_WS01.GT.15m60m_15","MA_WS07.GT.15m30m60m_60",
                                     "ni_ind10.04.LN_10" ,"ni_LJ.Avanti02PP.5min_60","ni_LJ.MB02PP.5min_30" ,"ni_LJ.ninexiuPP.5min_60","ni_LJ.TW02PP.1H_15","OI_01.4H_240","OI_02.2H_120","OI_LJ.Kelther.Tsi01T.30min_240",
                                     "OI_LJ.multsig03PP.15min_120","OI_LJ.multsig.A.4H_120","OI_LJ.multsig.A.4H_240","OI_LJ.threeswordPP.5min_30")]
ind = endpoints(rt,on='month')
#tmp.rt=rt
rt=tmp.rt

# for (i in 1:ncol(rt)){
#   char = strsplit(colnames(rt)[i],'_')
#   char = char[[1]]
#   k = which(colnames(total.contract)==char[1])
#   #print(k)
#   for (j in 1:(length(ind)-1)){
#     #print(j)
#     if (j ==1){
#       rt[(ind[j]+1):ind[j+1],i]=rt[(ind[j]+1):ind[j+1],i]*total.contract[j,k]*contract$coef[k]}
#     else{
#       rt[(ind[j]+1):ind[j+1],i]=rt[(ind[j]+1):ind[j+1],i]*total.contract[j,k]*contract$coef[k]
#       if (rt[(ind[j]),i]!=0){
#         rt[(ind[j]+1),i]=rt[(ind[j]+1),i]+(total.contract[j,k]-total.contract[j+1,k])*contract$cost[k]
#         }
#     }
#   }
# }



X_lin <- (prices/lag(prices) - 1)[-1]


####################################################################


total.num <- matrix(ncol=96,nrow=19)
colnames(total.num) = colnames(Strategy.returns[,3:ncol(Strategy.returns)])
rownames(total.num) = rownames(price)



