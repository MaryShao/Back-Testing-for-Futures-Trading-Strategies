setwd("C:/Users/m8sha/Desktop")

library(xts)
library(quantmod)
library(PerformanceAnalytics)
#install.packages('psych')
library('psych')
#install.packages('clusterSim')
library('clusterSim')
#install.packages('ade4')
library('ade4')
library('PerformanceAnalytics')

PP <- read.table("C:/Users/m8sha/Desktop/PP.txt", quote="\"", comment.char="")
TT <- read.table("C:/Users/m8sha/Desktop/TT.txt", quote="\"", comment.char="")



## JSD
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

## ==============================minute =========================================##
sig = as.matrix(Signals.PP[,3:13])
p = Signals.PP[,2]
retn = matrix(0,nrow=129633,ncol = 11)
colnames(retn) = colnames(sig)
rownames(retn) =Signals.PP[,1]

for(j in 1:11){
  for (i in 1:129632) {
    if ((sig[i,j]==1 && sig[i+1,j]==1)|
        (sig[i,j]==0 && sig[i+1,j]==-1)|
        (sig[i,j]==1 && sig[i+1,j]==0)|
        (sig[i,j]==1 && sig[i+1,j]==-1)){
      retn[i,j]= log(p[i+1]/p[i])
    }
    if((sig[i,j]==-1 && sig[i+1,j]==-1)|
       (sig[i,j]==-1 && sig[i+1,j]==0)|
       (sig[i,j]==0 && sig[i+1,j]== 1)|
       (sig[i,j]==-1 && sig[i+1,j]==1)){
      retn[i,j]= log(p[i]/p[i+1])
    }
  }
}


# Signal.profit<-function(signal.tick,cost){
#   Returns<-numeric()
#   for (i in (1:(nrow(signal.tick)-1))) {
#     if (signal.tick$Signal[i]=="0") {
#       Returns[i]<- 0
#     }
#     if (signal.tick$Signal[i]=="1"&signal.tick$Signal[i+1]=="0"){
#       Returns[i]<- signal.tick$Close[i+1]- signal.tick$Close[i]-2*cost
#     }
#     if (signal.tick$Signal[i]=="1"&signal.tick$Signal[i+1]=="-1"){
#       Returns[i]<- signal.tick$Close[i+1]- signal.tick$Close[i]-2*cost
#     }
#     if (signal.tick$Signal[i]=="-1"&signal.tick$Signal[i+1]=="0"){
#       Returns[i]<- signal.tick$Close[i]-signal.tick$Close[i+1]-2*cost
#     }
#     if (signal.tick$Sgnal[i]=="-1"&signal.tick$Signal[i+1]=="1"){
#       Returns[i]<- signal.tick$Close[i]-signal.tick$Close[i+1]-2*cost
#     }
#     if (signal.tick$Unit[i]=="1"&signal.tick$Unit[i+1]=="1"){
#       Returns[i]<- signal.tick$Close[i+1]- signal.tick$Close[i]
#     }
#     if (signal.tick$Signal[i]=="-1"&signal.tick$Signal[i+1]=="-1"){
#       Returns[i]<- signal.tick$Close[i]- signal.tick$Close[i+1]
#     }
#   }
#   return(Returns)
# }

#################### hourly return ###############################################

retn.hour = retn[substr(Signals.PP[,1], start = 15, stop = 19) == "00:00",]


k<-c()
for(i in 1:129633){
  if (substr(Signals.PP[i,1], start = 15, stop = 19) == "00:00"){
      k<-c(k,i)
  }
}

retn.hour = rbind(retn[1,],matrix(0,nrow=129632,ncol = 11))
for (i in (2:129633)){
  for (j in 1:11) {
    retn.hour[i,j]= retn.hour[i-1,j]+retn[i,j]
  }
}

hourly.retn<-retn.hour[k[1],]

names<-rownames(retn)[k[1]]
for (i in 2:length(k)){
  r<- retn.hour[k[i],]-retn.hour[k[i-1],]
  hourly.retn<-rbind(hourly.retn,r)
  names<-c(names,rownames(retn)[k[i]])
}
rownames(hourly.retn)=names

############# cumulative hourly return ########################################
#c<-as.data.frame(cor(hourly.retn.cum))

hourly.retn.cum<-retn.hour[k[1],]

names<-rownames(retn)[k[1]]
for (i in 2:length(k)){
  r<- retn.hour[k[i],]
  hourly.retn.cum<-rbind(hourly.retn.cum,r)
  names<-c(names,rownames(retn)[k[i]])
}
rownames(hourly.retn.cum)=names



################## Signal & return EMD ####################################


hourly.retn.cum.dist=dist(t(hourly.retn.cum),method ="euclidean",upper =T,diag=T)
plot(hclust(hourly.retn.dist,"average"))

############## RATIO #############################
ratio.retn<- rbind(diag(cov(hourly.retn)),
                   DRatio(hourly.retn),
                   SharpeRatio(hourly.retn)[1,])

ratio.retn.dist=dist(t(ratio.retn),method ="euclidean",upper =T,diag=T)
plot(hclust(ratio.retn.dist,"ward.D"))
#plot(hclust(ratio.retn.dist,"complete"))
#plot(hclust(ratio.retn.dist,"average"))
rt.pam <- pam(ratio.retn.dist, k=2, diss = TRUE)
plot(rt.pam)



