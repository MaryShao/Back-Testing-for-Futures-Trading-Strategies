library(tidyr)
library(clusterSim)


### DATA input ####################################################
setwd("C:/Users/m8sha/Desktop")
Strategy.signals<-read.table('Strategy.signals.txt',header=T,row.names=NULL,sep=' ')

d<-as.character.Date(rbind(Strategy.signals[,1:2]))
u=unite(rbind(d[,1:2]),col='Date',sep=' ',remove = FALSE)
Strategy.signals=Strategy.signals[3:length(colnames(Strategy.signals))]
rownames(Strategy.signals)=u[,1]

## JSD ######################################################################
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
sig = as.matrix(Strategy.signals)

# JSD (2018)
sig.18min = sig[1:87480,]
minProb = matrix(0,nrow=3,ncol = 82)
colnames(minProb) = colnames(Strategy.signals)
rownames(minProb) = c(0,-1,1)
for(i in 1:82){
  for (k in 1:3) {
    minProb[k,i] = sum(sig.18min[,i] == rownames(minProb)[k])/87480
  }
}

sig18Min.JSD = dist.JSD(minProb)
as.matrix(sig18Min.JSD)
plot(hclust(sig18Min.JSD,"average"))
plot(hclust(sig18Min.JSD,"complete"))

rt.pam <- pam(sig18Min.JSD, k=13, diss = TRUE)
plot(rt.pam)

# JSD (2019)
sig.19min = sig[87481:129633,]
minProb = matrix(0,nrow=3,ncol = 82)
colnames(minProb) = colnames(Strategy.signals)
rownames(minProb) = c(0,-1,1)
for(i in 1:82){
  for (k in 1:3) {
    minProb[k,i] = sum(sig.19min[,i] == rownames(minProb)[k])/42153
  }
}

sig19Min.JSD = dist.JSD(minProb)
as.matrix(sig19Min.JSD)
plot(hclust(sig19Min.JSD,"average"))
plot(hclust(sig19Min.JSD,"complete"))

rt.pam <- pam(sig19Min.JSD, k=13, diss = TRUE)
plot(rt.pam)


# JSD (18 19
sig.min = sig
minProb = matrix(0,nrow=3,ncol = 82)
colnames(minProb) = colnames(Strategy.signals)
rownames(minProb) = c(0,-1,1)
for(i in 1:82){
  for (k in 1:3) {
    minProb[k,i] = sum(sig.min[,i] == rownames(minProb)[k])/(dim(sig)[1])
  }
}

sigMin.JSD = dist.JSD(minProb)
as.matrix(sigMin.JSD)
plot(hclust(sigMin.JSD,"average"))
plot(hclust(sigMin.JSD,"complete"))

rt.pam <- pam(sigMin.JSD, k=13, diss = TRUE)
plot(rt.pam)

########====================== HOUR =========================######
sig.hour =Strategy.signals[substr(rownames(Strategy.signals), start = 15, stop = 19) == "00:00",]

# JSD (2018)
sig.18hour = sig.hour[1:1458,]
hour18Prob = matrix(0,nrow=3,ncol = 82)
colnames(hour18Prob) = colnames(Strategy.signals)
rownames(hour18Prob) = c(0,-1,1)
for(i in 1:82){
  for (k in 1:3) {
    hour18Prob[k,i] = sum(sig.18hour[,i] == rownames(hour18Prob)[k])/1458
  }
}

sig18Hour.JSD = dist.JSD(hour18Prob)
as.matrix(sig18Hour.JSD)
plot(hclust(sig18Hour.JSD,"average"))
plot(hclust(sig18Hour.JSD,"complete"))

rt.pam <- pam(sig18Hour.JSD, k=13, diss = TRUE)
plot(rt.pam)

# JSD 2019
sig.2019hour = sig.hour[1459:2160,]
hour19Prob = matrix(0,nrow=3,ncol = 82)
colnames(hour19Prob) = colnames(PP.signal)
rownames(hour19Prob) = c(0,-1,1)
for(i in 1:82){
  for (k in 1:3) {
    hour19Prob[k,i] = sum(sig.2019hour[,i] == rownames(hour19Prob)[k])/702
  }
}

sig19Hour.JSD = dist.JSD(hour19Prob)
as.matrix(sig19Hour.JSD)
plot(hclust(sig19Hour.JSD,"average"))
plot(hclust(sig19Hour.JSD,"complete"))

rt.pam <- pam(sig19Hour.JSD, k=13, diss = TRUE)
plot(rt.pam)

# JSD 18-19
sighour = sig.hour
hourProb = matrix(0,nrow=3,ncol = 82)
colnames(hourProb) = colnames(Strategy.signals)
rownames(hourProb) = c(0,-1,1)
for(i in 1:82){
  for (k in 1:3) {
    hourProb[k,i] = sum(sighour[,i] == rownames(hourProb)[k])/2160
  }
}

sigHour.JSD = dist.JSD(hourProb)
as.matrix(sigHour.JSD)
plot(hclust(sigHour.JSD,"average"))
plot(hclust(sigHour.JSD,"complete"))

rt.pam <- pam(sigHour.JSD, k=13, diss = TRUE)
plot(rt.pam)












######################################################################################
# Times ratio
