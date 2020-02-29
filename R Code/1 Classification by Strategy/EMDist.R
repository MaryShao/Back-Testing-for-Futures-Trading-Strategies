install.packages('emdist')
library(emdist)


returns<-sample(-10:10,18)
o<-numeric()
for (i in 1:length(returns)){
  if (returns[i]==0)
    o[i]<-0
  if (returns[i]>0)
    o[i]<-1
  if (returns[i]<0)
    o[i]<--1
}
l<-rep(1,18)
s<-rep(-1,18)
r<-sample(c(1,-1),18,T)
returns.op<-abs(returns)/sum(abs(returns))

returns.l<-diag(returns%*%t(l))
returns.l[which(returns.l<0)]<-0
returns.L<-returns.l/sum(returns.l)

returns.s<-diag(returns%*%t(s))
returns.s[which(returns.s<0)]<-0
returns.S<-returns.s/sum(returns.s)

returns.r<-diag(returns%*%t(r))
returns.r[which(returns.r<0)]<-0
returns.R<-returns.r/sum(returns.r)


Returns.o<-matrix(c(returns.op,o),ncol=2)
Returns.l<-matrix(c(returns.L,rep(1,18)),ncol=2)
Returns.s<-matrix(c(returns.S,rep(-1,18)),ncol=2)
Returns.r<-matrix(c(returns.R,r),ncol=2)


Returns.O<-matrix(c(returns.op,o),ncol=2)
Returns.L<-matrix(c(returns.L,returns.op),ncol=2)
Returns.S<-matrix(c(returns.S,returns.op),ncol=2)
Returns.R<-matrix(c(returns.R,returns.op),ncol=2)



emdist = matrix(0,nrow=4,ncol = 4)
rownames(emdist)=colnames(emdist)=c('o','l','s','r')
emdist[1,2] = emdist[2,1]=emd(Returns.o, Returns.l, dist="euclidean")
emdist[1,3] = emdist[3,1] = emd(Returns.o, Returns.s, dist="euclidean")
emdist[1,4] = emdist[4,1] = emd(Returns.o, Returns.r, dist="euclidean")
emdist[2,3] = emdist[3,2] = emd(Returns.s, Returns.l, dist="euclidean")
emdist[2,4] = emdist[4,2] = emd(Returns.r, Returns.l, dist="euclidean")
emdist[3,4] = emdist[4,3] = emd(Returns.s, Returns.r, dist="euclidean")
diag(emdist)=0


emdist2 = matrix(0,nrow=4,ncol = 4)
rownames(emdist2)=colnames(emdist2)=c('O','L','S','R')
emdist2[1,2] = emdist2[2,1] = emd(Returns.O, Returns.L, dist="euclidean")
emdist2[1,3] = emdist2[3,1] = emd(Returns.O, Returns.S, dist="euclidean")
emdist2[1,4] = emdist2[4,1] = emd(Returns.O, Returns.R, dist="euclidean")
emdist2[2,3] = emdist2[3,2] = emd(Returns.S, Returns.L, dist="euclidean")
emdist2[2,4] = emdist2[4,2] = emd(Returns.R, Returns.L, dist="euclidean")
emdist2[3,4] = emdist2[4,3] = emd(Returns.S, Returns.R, dist="euclidean")
diag(emdist2)=0

emdist
emdist2


cor(cbind(returns.L,returns.R,returns.S))
cor(cbind(Returns.O,Returns.L,Returns.R,Returns.S))
cov(cbind(l,s,r))


################# Signal test ####################################
l<-matrix(1,nrow=129633,ncol = 11)
w<-matrix(1,nrow=129633,ncol = 11)
r<-matrix(0,nrow=129633,ncol = 11)
t<-matrix(-1,nrow=129633,ncol = 11)
for (i in 1:2160){
  if (i %% 2 ==0){
    l[i,] = -1
    w[i,] = 0
    r[i,] = -1
    t[i,] = 1
  }
}

#sig = as.matrix(Signals.PP[,3:13])
sig=r
p = Signals.PP[,2]

o<-numeric()
for (i in 1:length(p)){
  if (p[i]==0)
    o[i]<-0
  if (p[i]>0)
    o[i]<-1
  if (p[i]<0)
    o[i]<--1
}


# retn = matrix(0,nrow=129633,ncol=11)
retn.o = matrix(0,nrow=129633,ncol=1)
# retn.l = matrix(0,nrow=129633,ncol=1)
# retn.w = matrix(0,nrow=129633,ncol=1)
# retn.r = matrix(0,nrow=129633,ncol=1)
# retn.t = matrix(0,nrow=129633,ncol=1)
colnames(retn) = colnames(sig)
rownames(retn) =Signals.PP[,1]

for(j in 1:11){
  for (i in 1:129632) {
    if ((sig[i,j]==1 && sig[i+1,j]==1)|
        (sig[i,j]==1 && sig[i+1,j]==0)|
        (sig[i,j]==1 && sig[i+1,j]==-1)){
      #retn[i,j]= log(p[i+1]/p[i])
      # retn.l[i]= log(p[i+1]/p[i])
      # retn.w[i]= log(p[i+1]/p[i])
      retn.r[i]= log(p[i+1]/p[i])
      # retn.t[i]= log(p[i+1]/p[i])
    }
    if((sig[i,j]==-1 && sig[i+1,j]==-1)|
       (sig[i,j]==-1 && sig[i+1,j]==0)|
       (sig[i,j]==-1 && sig[i+1,j]==1)){
     # retn[i,j]= log(p[i]/p[i+1])
      # retn.l[i]= log(p[i]/p[i+1])
      # retn.w[i]= log(p[i]/p[i+1])
      retn.r[i]= log(p[i]/p[i+1])
      # retn.t[i]= log(p[i]/p[i+1])
    }
  }
}


################ hourly #######################################################
#retn.hour = retn[substr(Signals.PP[,1], start = 15, stop = 19) == "00:00",]


k<-c()
for(i in 1:129633){
  if (substr(Signals.PP[i,1], start = 15, stop = 19) == "00:00"){
    k<-c(k,i)
  }
}

retn.hour = rbind(retn[1,],matrix(0,nrow=129632,ncol = 11))
retn.hour.l = c(retn.l[1],rep(0,129632))
retn.hour.w = c(retn.w[1],rep(0,129632))
retn.hour.r = c(retn.r[1],rep(0,129632))
retn.hour.t = c(retn.t[1],rep(0,129632))
for (i in (2:129633)){
  for (j in 1:11) {
    retn.hour[i,j]= retn.hour[i-1,j]+retn[i,j]
  }
  retn.hour.l[i] = retn.hour.l[i-1]+retn.l[i]
  retn.hour.w[i] = retn.hour.w[i-1]+retn.w[i]
  retn.hour.r[i] = retn.hour.r[i-1]+retn.r[i]
  retn.hour.t[i] = retn.hour.t[i-1]+retn.t[i]
}

hourly.retn<-retn.hour[k[1],]
hourly.retn.l = retn.hour.l[k[1]]
hourly.retn.w = retn.hour.w[k[1]]
hourly.retn.r = retn.hour.r[k[1]]
hourly.retn.t = retn.hour.t[k[1]]

names<-rownames(retn)[k[1]]
for (i in 2:length(k)){
  r<- retn.hour[k[i],]-retn.hour[k[i-1],]
  hourly.retn<-rbind(hourly.retn,r)
  names<-c(names,rownames(retn)[k[i]])
  hourly.retn.l = c(hourly.retn.l,retn.hour.l[k[i]]-retn.hour.l[k[i-1]])
  hourly.retn.w = c(hourly.retn.w,retn.hour.w[k[i]]-retn.hour.w[k[i-1]])
  hourly.retn.r = c(hourly.retn.r,retn.hour.r[k[i]]-retn.hour.r[k[i-1]])
  hourly.retn.t = c(hourly.retn.t,retn.hour.t[k[i]]-retn.hour.t[k[i-1]])
}
rownames(hourly.retn)=names

############# cumulative hourly return ########################################
#c<-as.data.frame(cor(hourly.retn.cum))

hourly.retn.cum<-retn.hour[k[1],]
hourly.retn.l.cum = retn.hour.l[k[1]]
hourly.retn.w.cum = retn.hour.w[k[1]]
hourly.retn.r.cum = retn.hour.r[k[1]]
hourly.retn.t.cum = retn.hour.t[k[1]]

names<-rownames(retn)[k[1]]
for (i in 2:length(k)){
  r<- retn.hour[k[i],]
  hourly.retn.cum<-rbind(hourly.retn.cum,r)
  names<-c(names,rownames(retn)[k[i]])
  hourly.retn.l.cum = c(hourly.retn.l.cum,retn.hour.l[k[i]])
  hourly.retn.w.cum = c(hourly.retn.w.cum,retn.hour.w[k[i]])
  hourly.retn.r.cum = c(hourly.retn.r.cum,retn.hour.r[k[i]])
  hourly.retn.t.cum = c(hourly.retn.t.cum,retn.hour.t[k[i]])
}
rownames(hourly.retn.cum)=names



################## Signal & return EMD ####################################

### l w r t
RL = as.matrix(cbind(1000*hourly.retn.l,l[1:2160]))
RW = as.matrix(cbind(1000*hourly.retn.w,w[1:2160]))
RR = as.matrix(cbind(1000*hourly.retn.r,r[1:2160]))
RT = as.matrix(cbind(1000*hourly.retn.t,t[1:2160]))

# hourly return
emdist3 = matrix(0,nrow=4,ncol = 4)
rownames(emdist3)=colnames(emdist3)=c('l(1,-1)','w(1,0)','s(0,-1)','t(-1,1)')
emdist3[1,2] = emdist3[2,1] = emd(RL,RW,dist="euclidean")
emdist3[1,3] = emdist3[3,1] = emd(RL,RR, dist="euclidean")
emdist3[1,4] = emdist3[4,1] = emd(RL,RT, dist="euclidean")
emdist3[2,3] = emdist3[3,2] = emd(RW,RR, dist="euclidean")
emdist3[2,4] = emdist3[4,2] = emd(RW,RT, dist="euclidean")
emdist3[3,4] = emdist3[4,3] = emd(RR,RT, dist="euclidean")
diag(emdist3)=0
emdist3

c<-cor(as.matrix(cbind(hourly.retn.l,hourly.retn.w,hourly.retn.r,hourly.retn.t)))

Emd = matrix(0,)



returns.l<-hourly.retn%*%t(matrix.l)[,1]
returns.l[which(returns.l<0)]<-0 
returns.L<-returns.l/sum(returns.l)


returns.s<-hourly.retn%*%t(matrix.s)[,1]
returns.s[which(returns.s<0)]<-0
returns.S<-returns.s/sum(returns.s)

returns.r<-hourly.retn%*%t(matrix.r)[,1]
returns.r[which(returns.r<0)]<-0
returns.R<-returns.r/sum(returns.r)

sig = Signals.PP[substr(Signals.PP[,1], start = 15, stop = 19) == "00:00",][,3:13]
returns.O<-c()
for (i in 1:11){
  returns.O<-c(hourly.retn[,i]%*%t(sig[,i]))
}

returns.o<-matrix(c(returns.op,o),ncol=2)
returns.L<-matrix(c(returns.L,rep(1,18)),ncol=2)
returns.S<-matrix(c(returns.S,rep(-1,18)),ncol=2)
returns.R<-matrix(c(returns.R,r),ncol=2)

emdist3 = matrix(0,nrow=4,ncol = 4)
rownames(emdist3)=colnames(emdist3)=c('O','L','S','R')
emdist3[1,2] = emdist3[2,1] = emd(returns.O, returns.L, dist="euclidean")
emdist3[1,3] = emdist3[3,1] = emd(returns.O, returns.S, dist="euclidean")
emdist3[1,4] = emdist3[4,1] = emd(returns.O, returns.R, dist="euclidean")
emdist3[2,3] = emdist3[3,2] = emd(returns.S, returns.L, dist="euclidean")
emdist3[2,4] = emdist3[4,2] = emd(returns.R, returns.L, dist="euclidean")
emdist3[3,4] = emdist3[4,3] = emd(returns.S, returns.R, dist="euclidean")
diag(emdist3)=0
emdist3
