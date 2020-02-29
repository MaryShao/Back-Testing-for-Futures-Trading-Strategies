
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


setwd("C:/Users/m8sha/Desktop/1min")
f = list.files()
#f = f[-2]

dat = matrix(nrow = 158146,ncol = 1)

for(i in f){
  tmp = read.csv(i, header = TRUE)
  tmp[,1] = as.character(tmp[,1])
  tmp[,2] = as.character(tmp[,2])
  dat = cbind(dat,tmp[,6])
}


colnames(dat) =c('CFFEX.IF',
                 'CFFEX.T',
                 'CZCE.MA',
                 'CZCE.OI',
                 'CZCE.SR',
                 'CZCE.ZC',
                 'DCE.a',
                 'DCE.i',
                 'DCE.j',
                 'DCE.jm',
                 'DCE.pp',
                 'DCE.v',
                 'SHFE.al',
                 'SHFE.hc',
                 'SHFE.ni',
                 'SHFE.rb',
                 'SHFE.zn')
rownames(dat) = tmp[,1]

#dat.scale = dat
#for (i in 1:17) {
 # minimum = min(dat[,i])
#  maximum = max(dat[,i])
#  dat.scale[,i] = (dat[,i] - minimum)/(maximum - minimum)
#}

########################BENCHMARK###############################
benchmark = matrix(0,nrow = dim(dat)[1]-1, ncol = dim(dat)[2])
colnames(benchmark) =c('CFFEX.IF',
                        'CFFEX.T',
                        'CZCE.MA',
                        'CZCE.OI',
                        'CZCE.SR',
                        'CZCE.ZC',
                        'DCE.a',
                        'DCE.i',
                        'DCE.j',
                        'DCE.jm',
                        'DCE.pp',
                        'DCE.v',
                        'SHFE.al',
                        'SHFE.hc',
                        'SHFE.ni',
                        'SHFE.rb',
                        'SHFE.zn')
rownames(benchmark) = tmp[-1,1]

for (i in 1:17) {
  benchmark[,i] = rep(0.035,341)
}
########################BENCHMARK###############################


dat<-as.xts(dat)
r_log <- CalculateReturns(dat, "log")[-1]
r_lin <- CalculateReturns(dat)[-1]
benchmark<-as.xts(benchmark)



# calculate performance###################################################################

perf= matrix(nrow = 0,ncol = 17)
colnames(perf) =c('CFFEX.IF',
                 'CFFEX.T',
                 'CZCE.MA',
                 'CZCE.OI',
                 'CZCE.SR',
                 'CZCE.ZC',
                 'DCE.a',
                 'DCE.i',
                 'DCE.j',
                 'DCE.jm',
                 'DCE.pp',
                 'DCE.v',
                 'SHFE.al',
                 'SHFE.hc',
                 'SHFE.ni',
                 'SHFE.rb',
                 'SHFE.zn')


mu_lin <- colMeans(r_lin)
sigma_lin <- cov(r_lin)
mu_lin
sigma_lin

mu_log <- colMeans(r_log)
sigma_log <- cov(r_log)
mu_log
sigma_log

# create function for GMVP 
portolioGMVP <- function(Sigma) { 
  ones <- rep(1, nrow(Sigma))  
  Sigma_inv_1 <- solve(Sigma, ones)  #same as: inv(Sigma) %*% ones  
  w <- (1/as.numeric(ones %*% Sigma_inv_1)) * Sigma_inv_1
  return(w) } 
# compute the three versions of GMVP 
portolio_lin <- portolioGMVP(sigma_lin) 
portolio_log <- portolioGMVP(sigma_log)

# performance measures
table.AnnualizedReturns(r_log)

# BernardoLedoitRatio Bernardo and Ledoit ratio of the return distribution
BernardoLedoitRatio(r_log)

# Burke ratio of the return distribution
BurkeRatio(r_log, Rf = 0, modified = FALSE)

# calculate a Calmar or Sterling reward/risk ratio
CalmarRatio(r_log, scale = NA)

# Calculate Uryasev’s proposed Conditional Drawdown at Risk (CDD or CDaR) measure
CDD(r_log, weights = NULL, geometric = TRUE, invert = TRUE, p = 0.95)

#  downside risk (deviation, variance) of the return distribution
DownsideDeviation(r_log, MAR = 0, method = c("full", "subset"),potential = FALSE)
# downside frequency of the return distribution
DownsideFrequency(r_log, MAR = 0)

# d ratio of the return distribution
DRatio(r_log)

# Calculates a standard deviation-type statistic using individual drawdowns.
DrawdownDeviation(r_log)

# InformationRatio = ActivePremium/TrackingError
information_ratio<-c()
for (i in 1:17) {
  information_ratio<-c(information_ratio,InformationRatio(r_log[,i], benchmark[,i], scale = 252))
}
InformationRatio(r_log[,1:17], benchmark, scale = 252)

# calculate Kelly criterion ratio (leverage or bet size) for a strategy
KellyRatio(r_log, Rf = 0, method = "half")

#Martin ratio of the return distribution
MartinRatio(r_log, Rf = 0)

# caclulate the maximum drawdown from peak equity
maxDrawdown(r_log, weights = NULL, geometric = TRUE, invert = TRUE)

#Mean absolute deviation of the return distribution
MeanAbsoluteDeviation(r_log)


#
Omega(r_log, L = 0, method = "simple", output = "point", Rf = 0)

#
PainRatio(r_log, Rf = 0)

#
MAR=0.05
ProspectRatio(r_log,MAR)

#
SharpeRatio(r_log, Rf = 0, p = 0.95, FUN = "StdDev")
#
SkewnessKurtosisRatio(r_log)
#
SmoothingIndex(r_log, neg.thetas = FALSE, MAorder = 2, verbose = FALSE)
#
SortinoRatio(r_log, MAR)
#
UlcerIndex(r_log)
#
UpsideFrequency(r_log, MAR = 0)
#
VaR(r_log, p = 0.95,method = "gaussian", portfolio_method ="single",invert = TRUE)


##3##################################################################################

inf_ratio<- InformationRatio(r_log[,1:17], benchmark, scale = 252)
inf_ratio=c(inf_ratio[2,1],inf_ratio[1,2:17])
inf_ratio
MAR=0.05

perf= rbind(perf,
            colMeans(r_log),
            diag(cov(r_log)),
            table.AnnualizedReturns(r_log),
            BernardoLedoitRatio(r_log),
            BurkeRatio(r_log, Rf = 0, modified = FALSE),
            CalmarRatio(r_log, scale = NA),
            CDD(r_log, weights = NULL, geometric = TRUE, invert = TRUE, p = 0.95),
            DownsideDeviation(r_log, MAR = 0, method = c("full", "subset"),potential = FALSE),
            DRatio(r_log),
            DrawdownDeviation(r_log),
            inf_ratio,
            KellyRatio(r_log, Rf = 0, method = "half"),
            MartinRatio(r_log, Rf = 0),
            maxDrawdown(r_log, weights = NULL, geometric = TRUE, invert = TRUE),
            MeanAbsoluteDeviation(r_log),
            Omega(r_log, L = 0, method = "simple", output = "point", Rf = 0),
            PainRatio(r_log, Rf = 0),
            ProspectRatio(r_log,MAR),
            SharpeRatio(r_log, Rf = 0, p = 0.95, FUN = "StdDev"),
            SkewnessKurtosisRatio(r_log),
            SmoothingIndex(r_log, neg.thetas = FALSE, MAorder = 2, verbose = FALSE),
            SortinoRatio(r_log, MAR),
            UlcerIndex(r_log),
            UpsideFrequency(r_log, MAR = 0),
            VaR(r_log, p = 0.95,method = "gaussian", portfolio_method ="single",invert = TRUE))

rownames(perf)[1:2] = c('mu','cov')
rownames(perf)[13] = c('Information ratio')

# 计算特征值和特征向量##按指标###################################################################

perf.scale = perf
for (i in 1:17) {
  minimum = min(perf[i,])
  maximum = max(perf[i,])
  perf[i,] = (perf[i,] - minimum)/(maximum - minimum)
}

perf.dist=dist(perf,method ="euclidean",upper =T,diag=T)

sigma<-cor(as.matrix(perf.dist))
e<-eigen(sigma)
#特征值
e$values
#特征向量
e$vectors

# 确定主成分个数
#画碎石图
fa.parallel(as.matrix(sigma),n.obs=28,fa='pc')
abline(h=1)

# 提取主成分
fit<-principal(sigma,
               nfactors=3,
               rotate='varimax',  # max variance
               scores=T)

fit

# 绘制主成分的载荷矩阵，查看各个主成分的综合构成变量
fa.diagram(fit,digits=2)

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
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

perf.dist=dist.JSD(t(as.matrix(perf.dist)))

# optimal # of clusters

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

perf.cluster=pam.clustering(t(as.matrix(perf.dist)), k=3)

nclusters = index.G1(t(perf), perf.cluster, t(as.matrix(perf.dist)), centrotypes = "medoids") 
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

perf.cluster=pam.clustering(t(as.matrix(perf.dist)), k=4)
obs.pcoa=dudi.pco(perf.dist, scannf=F, nf=4)
s.class(obs.pcoa$li, fac=as.factor(perf.cluster), grid=F) 
plot(hclust(perf.dist))

s.class(obs.pcoa$li, fac=as.factor(perf.cluster), grid=F, cell=0, label = levels(as.factor(dat.cluster)),cstar=0, col=c(3,2,4)) 
#############################################################################################

# 计算特征值和特征向量##按策略###################################################################
perf.dist2=dist(t(perf),method ="euclidean",upper =T,diag=T)

sigma2<-cor(as.matrix(perf.dist2))
e2<-eigen(sigma2)
#特征值
e2$values
#特征向量
e2$vectors

# 确定主成分个数
#画碎石图
fa.parallel(as.matrix(sigma2),n.obs=17,fa='pc')
abline(h=1)

# 提取主成分
fit<-principal(sigma2,
               nfactors=2,
               rotate='varimax',  # max variance
               scores=T)

fit

# 绘制主成分的载荷矩阵，查看各个主成分的综合构成变量
fa.diagram(fit,digits=2)



perf.dist2=dist.JSD(t(as.matrix(perf.dist2)))

# optimal # of clusters
perf.cluster2=pam.clustering(t(as.matrix(perf.dist2)), k=2)

nclusters = index.G1(t(perf), perf.cluster2, t(as.matrix(perf.dist2)), centrotypes = "medoids") 
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

perf.cluster=pam.clustering(t(as.matrix(perf.dist2)), k=4)
obs.pcoa=dudi.pco(perf.dist2, scannf=F, nf=4)
s.class(obs.pcoa$li, fac=as.factor(perf.cluster2), grid=F) 
plot(hclust(perf.dist2))



######  chart  ########################################################
chart.ACFplus(r_log, maxlag = NULL, elementcolor = "gray", main = NULL)

charts.PerformanceSummary(r_log, Rf = 0, main = NULL, geometric = TRUE, 
                          methods = "none", width = 0, event.labels = NULL, 
                          ylog = FALSE, wealth.index = FALSE, gap = 12, 
                          begin = c("first", "axis"), legend.loc = "topleft", p = 0.95)

hclust(perf.dist)

