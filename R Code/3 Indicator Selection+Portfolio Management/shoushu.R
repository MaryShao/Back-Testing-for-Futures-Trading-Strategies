w4_1<- read.csv("month19start.csv", row.names=1, stringsAsFactors=FALSE)
w4_2<- read.csv("month19block.csv", row.names=1, stringsAsFactors=FALSE)
w3_1<- read.csv("thirdFRI19start.csv", row.names=1, stringsAsFactors=FALSE)
w3_2<- read.csv("thirdFRI19block.csv", row.names=1, stringsAsFactors=FALSE)
w2_1<- read.csv("secFRI19start.csv", row.names=1, stringsAsFactors=FALSE)
w2_2<- read.csv("secFRI19block.csv", row.names=1, stringsAsFactors=FALSE)
w1_1<- read.csv("firstFRI19start.csv", row.names=1, stringsAsFactors=FALSE)
w1_2<- read.csv("firstFRI19block.csv", row.names=1, stringsAsFactors=FALSE)

monthPrice <- as.matrix(read_xlsx("month.xlsx")[,-1])
rownames(monthPrice) <-c("2018-12-28","2019-01-31","2019-02-28","2019-03-29","2019-04-30","2019-05-31")
thirdPrice <- as.matrix(read_xlsx("third.xlsx")[,-1])
rownames(thirdPrice) <-c("2018-12-28","2019-01-21","2019-02-18","2019-03-18","2019-04-22","2019-05-20","2019-06-24")
secPrice <- as.matrix(read_xlsx("second.xlsx")[,-1])
rownames(secPrice) <-c("2018-12-28","2019-01-14","2019-02-11","2019-03-11","2019-04-15","2019-05-13","2019-06-17")
firstPrice <- as.matrix(read_xlsx("first.xlsx")[,-1])
rownames(firstPrice) <-c("2018-12-28","2019-01-07","2019-02-01","2019-03-04","2019-04-08","2019-05-06","2019-06-10")

coe = matrix(c(10,5,10,100,300,100,5,10,10,10000,100,5,60 ),ncol = 1)
rownames(coe) = c("a_","al","hc","i_","IF","j_","PP","rb","SR","T_","ZC","zn","jm")

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



monthStart <- shoushuStra(w4_1,monthPrice)
monthBlock <- shoushuStra(w4_2,monthPrice)

thirdStart <- shoushuStra(w3_1,thirdPrice)
thirdBlock <- shoushuStra(w3_2,thirdPrice)

secondStart <- shoushuStra(w2_1,secPrice)
secondBlock <- shoushuStra(w2_2,secPrice)

firstStart <- shoushuStra(w1_1,firstPrice)
firstBlock <- shoushuStra(w1_2,firstPrice)


write.csv(round(monthStart$shoushu),"shoushu/monthStartShoushu.csv")
write.csv(round(monthBlock$shoushu),"shoushu/monthBlockShoushu.csv")

write.csv(round(thirdStart$shoushu),"shoushu/thirdStartShoushu.csv")
write.csv(round(thirdBlock$shoushu),"shoushu/thirdBlockShoushu.csv")

write.csv(round(secondStart$shoushu),"shoushu/secondStartShoushu.csv")
write.csv(round(secondBlock$shoushu),"shoushu/secondBlockShoushu.csv")

write.csv(round(firstStart$shoushu),"shoushu/firstStartShoushu.csv")
write.csv(round(firstBlock$shoushu),"shoushu/firstBlockShoushu.csv")

write.csv(monthStart$capital,"CAPmonthStart.csv")
write.csv(monthBlock$capital,"CAPmonthBlock.csv")
write.csv(thirdStart$capital,"CAPthirdStart.csv")
write.csv(thirdBlock$capital,"CAPthirdBlock.csv")
write.csv(secondStart$capital,"CAPsecondStart.csv")
write.csv(secondBlock$capital,"CAPsecondBlock.csv")
write.csv(firstStart$capital,"CAPfirstStart.csv")
write.csv(firstBlock$capital,"CAPfirstBlock.csv")

############################ Contract #######################################
w4_1C<- read.csv("month19startC.csv", row.names=1, stringsAsFactors=FALSE)
w4_2C<- read.csv("month19blockC.csv", row.names=1, stringsAsFactors=FALSE)
w3_1C<- read.csv("thirdFRI19startC.csv", row.names=1, stringsAsFactors=FALSE)
w3_2C<- read.csv("thirdFRI19blockC.csv", row.names=1, stringsAsFactors=FALSE)
w2_1C<- read.csv("secFRI19startC.csv", row.names=1, stringsAsFactors=FALSE)
w2_2C<- read.csv("secFRI19blockC.csv", row.names=1, stringsAsFactors=FALSE)
w1_1C<- read.csv("firstFRI19startC.csv", row.names=1, stringsAsFactors=FALSE)
w1_2C<- read.csv("firstFRI19blockC.csv", row.names=1, stringsAsFactors=FALSE)


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


numT = 11
monthStartC <- shoushuStraC(w4_1C,monthPrice,numT)
monthBlockC <- shoushuStraC(w4_2C,monthPrice,numT)

thirdStartC <- shoushuStraC(w3_1C,thirdPrice,numT)
thirdBlockC <- shoushuStraC(w3_2C,thirdPrice,numT)

secondStartC <- shoushuStraC(w2_1C,secPrice,numT)
secondBlockC <- shoushuStraC(w2_2C,secPrice,numT)

firstStartC <- shoushuStraC(w1_1C,firstPrice,numT)
firstBlockC <- shoushuStraC(w1_2C,firstPrice,numT)


write.csv(round(monthStartC$shoushu),"shoushuC/monthStartShoushuC.csv")
write.csv(round(monthBlockC$shoushu),"shoushuC/monthBlockShoushuC.csv")

write.csv(round(thirdStartC$shoushu),"shoushuC/thirdStartShoushuC.csv")
write.csv(round(thirdBlockC$shoushu),"shoushuC/thirdBlockShoushuC.csv")

write.csv(round(secondStartC$shoushu),"shoushuC/secondStartShoushuC.csv")
write.csv(round(secondBlockC$shoushu),"shoushuC/secondBlockShoushuC.csv")

write.csv(round(firstStartC$shoushu),"shoushuC/firstStartShoushuC.csv")
write.csv(round(firstBlockC$shoushu),"shoushuC/firstBlockShoushuC.csv")



write.csv(monthStartC$capital,"CapC/CapmonthStartShoushuC.csv")
write.csv(monthBlockC$capital,"CapC/CapmonthBlockShoushuC.csv")

write.csv(thirdStartC$capital,"CapC/CapthirdStartShoushuC.csv")
write.csv(thirdBlockC$capital,"CapC/CapthirdBlockShoushuC.csv")

write.csv(secondStartC$capital,"CapC/CapsecondStartShoushuC.csv")
write.csv(secondBlockC$capital,"CapC/CapsecondBlockShoushuC.csv")

write.csv(firstStartC$capital,"CapC/CapfirstStartShoushuC.csv")
write.csv(firstBlockC$capital,"CapC/CapfirstBlockShoushuC.csv")



############################################# devide # in contract into strategies ######################################################
numStra = c(11,6,7,8,10,6,11,7,8,9,11,9,7)
numPf<-7
splitShoushu <- function(file,monthPrice){
  cap<-matrix(nrow = 13)
  
    shoushu <- read.csv(paste("original C/shoushu",file,sep = "/"),row.names=1)
    for (i in 1:nrow(shoushu)) {
      shoushu[i,] <- shoushu[i,]/numStra[i]
    }
    for (j in 1:ncol(shoushu)) {
      if(min(shoushu[shoushu[,j]>0,j])<1){
        shoushu[,j] <- shoushu[,j]*(1/min(shoushu[shoushu[,j]>0,j]))
      }
    }
    for(k in 1:nrow(monthPrice)){
      cap <- cbind(cap,coe[rownames(shoushu),1]*monthPrice[k,]*round(shoushu)[,((k-1)*numPf+1):(k*numPf)]*numStra)
    }
    write.csv(cap[,-1],paste("CAP",file))
    write.csv(round(shoushu),file)

}

splitShoushu("monthBlockShoushuC.csv",monthPrice)
splitShoushu("monthStartShoushuC.csv",monthPrice)
splitShoushu("thirdStartShoushuC.csv",thirdPrice)
splitShoushu("thirdBlockShoushuC.csv",thirdPrice)
splitShoushu("secondBlockShoushuC.csv",secPrice)
splitShoushu("secondStartShoushuC.csv",secPrice)
splitShoushu("firstBlockShoushuC.csv",firstPrice)
splitShoushu("firstStartShoushuC.csv",firstPrice)


nameStra <- rownames(w1_1)
numPf<-7
split2 <- function(file,monthPrice){

  shoushu <- read.csv(paste("original C/shoushu",file,sep = "/"),row.names=1)
  nameCon <- rownames(shoushu)
  ssNew <- matrix(0,nrow = length(nameStra),ncol=ncol(shoushu))
  rownames(ssNew) <- nameStra
  colnames(ssNew) <- colnames(shoushu)
  cap<-ssNew
  
  for(nam in nameCon){
    s <- nameStra[substr(nameStra,1,2)==nam]
    term <- sapply(strsplit(s,"_"),function(x){x[3]})
    n <- length(unique(term))
    for(t in unique(term)){
      rownm <- s[term==t]
      for (i in rownm) {
        ssNew[i,] <- as.matrix(shoushu[nam,]/n/sum(term == t))
      }
    }
  }
  
  for (j in 1:ncol(ssNew)) {
    if(min(ssNew[ssNew[,j]>0,j])<1){
      ssNew[,j] <- ssNew[,j]*(1/min(ssNew[ssNew[,j]>0,j]))
    }
  }
  for(k in 1:nrow(monthPrice)){
    for (c in 1:nrow(shoushu)) {
      cap[substr(rownames(ssNew),1,2)==rownames(shoushu)[c],((k-1)*numPf+1):(k*numPf)] <- 
        coe[rownames(shoushu)[c],1]*monthPrice[k,rownames(shoushu)[c]]*round(ssNew)[substr(rownames(ssNew),1,2)==rownames(shoushu)[c],((k-1)*numPf+1):(k*numPf)]*numStra[c]
    }
  }
  write.csv(cap,paste("CAP",file))
  write.csv(round(ssNew),file)
}

split2("monthBlockShoushuC.csv",monthPrice)
split2("monthBlockShoushuC.csv",monthPrice)
split2("monthStartShoushuC.csv",monthPrice)
split2("thirdStartShoushuC.csv",thirdPrice)
split2("thirdBlockShoushuC.csv",thirdPrice)
split2("secondBlockShoushuC.csv",secPrice)
split2("secondStartShoushuC.csv",secPrice)
split2("firstBlockShoushuC.csv",firstPrice)
split2("firstStartShoushuC.csv",firstPrice)









# price = Signals.T
# write.csv(cbind(price[rownames(price) %in%
# c("2018-12-28 15:15:00","2019-01-31 15:15:00","2019-02-28 15:15:00","2019-03-29 15:15:00","2019-04-30 15:15:00","2019-05-31 15:15:00"),1]
# ,
# price[rownames(price) %in%
#             c("2018-12-28 15:15:00","2019-01-21 15:15:00","2019-02-18 15:15:00","2019-03-18 15:15:00","2019-04-22 15:15:00","2019-05-20 15:15:00","2019-06-24 15:15:00"),1]
# ,
# price[rownames(price) %in%
#             c("2018-12-28 15:15:00","2019-01-14 15:15:00","2019-02-11 15:15:00","2019-03-11 15:15:00","2019-04-15 15:15:00","2019-05-13 15:15:00","2019-06-17 15:15:00"),1]
# ,
# price[rownames(price) %in%
#             c("2018-12-28 15:15:00","2019-01-07 15:15:00","2019-02-01 15:15:00","2019-03-04 15:15:00","2019-04-08 15:15:00","2019-05-06 15:15:00","2019-06-10 15:15:00"),1]
# 
# ),"price.csv")




