rm(list=ls())
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(scam)
source('SAEM-MCMC.R')
################################################
#-------------------SIMULATION --------------------
################################################"

A = 10; 
w = 85
m = 30
s = 7 
mu = log(c(A,w,m,s))
theta.true  = list(mu=mu,Omega=diag(c(0.01,0.01,0.02,0.08)),sigma2=5)

#------------------------------- 
lfmc.data <- read.table("data_lfmc.txt", header = TRUE)
myData.real <- lfmc.data[lfmc.data$leaf.type=='Grass E',]
nPlot <- length(unique(myData.real$plot))

#time <-rep(unique(myData.real$time),nPlot) 
time = myData.real$time
#plot <- rep(1:nPlot,each=length(unique(myData.real$time)))
plot <- myData.real$plot


#---------------------------------------- 

Z.sim = rmvnorm(nPlot,theta.true$mu,theta.true$Omega)
reg <- FourParamLogis(time,exp(Z.sim[plot,]))
y <- reg + sqrt(theta.true$sigma2)*rnorm(length(reg))
myData.sim <- data.frame(lfmc=y,reg=reg,plot=plot,plot.fact=as.factor(plot),time=time)
rm(myData.real)

g <- ggplot(data = myData.sim, aes(x = time, y = lfmc)) + 
  geom_point(aes(color=plot.fact))  + geom_smooth(colour="gray", linewidth=0.7)+
  xlab("Time (days)") + #+  facet_wrap(~leaf.type,scales = "free") + 
  ylab("LFMC (%)") # The LFMC temporal dynamics is plotted by "leaf type" 
g 

################################################ 
#########   initialisation 
############################################ 

Z.init <- initZ(myData.sim)
plot(Z.init,Z.sim)
abline(0,1)
theta.init <-initTheta(Z.init,myData.sim,omegaDiag = TRUE)

  
  



########################################################################
#################################### test MCMC 
########################################################################
paramsAlgo = list(nbIterMCMC=1000,rho=c(1/100,1/10,1,2),keep =TRUE,print=TRUE)
test_MCMC <- MCMC_NLME(Z.init,myData.sim,theta.init,paramsAlgo)


#------------------ plot des traces du MCMC
i <- sample(1:nPlot,1)
N = dim(test_MCMC)[1]
burn = round(N*0.25)
seq_iter = seq(burn+1,N,by=4)
par(mfrow=c(2,2))
title = c('A','w',"m","s")
for (k in 1:4){
 plot(seq_iter,test_MCMC[seq_iter,i,k],type='l',main=paste(title[k],'-',i))
 abline(h=Z.sim[i,k],col='red')
}
#---------------- posterior density 
# 
# 
# i <- sample(1:nPlot,1)
# N = dim(test_MCMC)[1]
# burn = round(N*0.25)
# seq_iter = seq(burn+1,N,by=10)
# par(mfrow=c(2,2))
# for (k in 1:4){
#    plot(density(test_MCMC[seq_iter,i,k]),main=paste(title[k],'-',i))
#    abline(v=Z.sim[i,k],col='red')
# }

###################################################""
#------------------- ESSAI SAEM 
######################################################## 
paramsAlgo = list(nbIterMCMC=100,nbIterSAEM = 100,rho=c(1/100,1/10,1,2),keep = FALSE,print=FALSE)
paramsAlgo$gammaSAEM = c(rep(1,10),1/(1:paramsAlgo$nbIterSAEM)^(1/1.1))
paramsAlgo$estim=list(mu=c(1,1,1,1),sigma2=1,Omega=1)

if(paramsAlgo$estim$Omega==0){theta.init$Omega =  theta.true$Omega}
if(paramsAlgo$estim$sigma2==0){theta.init$sigma2 =  theta.true$sigma2}
if(sum(paramsAlgo$estim$mu)<4){theta.init$mu[paramsAlgo$estim$mu==0] = theta.true$mu[paramsAlgo$estim$mu==0]}

res_SAEM_MCMC <- SAEM_MCMC_NLME(myData.sim,paramsAlgo,omegaDiag = TRUE)

names(res_SAEM_MCMC)
par(mfrow=c(2,2))
for (k in 1:4){
  yk <- (c(res_SAEM_MCMC$mu[,k], theta.true$mu[k],theta.init$mu[k]))
  plot((res_SAEM_MCMC$mu[,k]),type='l',ylim = range(yk),main=c('A','w',"m","s")[k])
  abline(h=(theta.true$mu[k]),col='red')
  abline(h=(theta.init$mu[k]),col='green')
  abline(h=colMeans(Z.sim)[k],col='orange')
}

par(mfrow=c(1,1))
ysigma <- (c(res_SAEM_MCMC$sigma2, theta.true$sigma2,theta.init$sigma))
plot(res_SAEM_MCMC$sigma2,type='l',ylim = range(ysigma),main=c('A','w',"m","s")[k])
abline(h=(theta.true$sigma2),col='red')
abline(h=(theta.init$sigma2),col='green')


################# posterior 
theta.estim = list(mu = res_SAEM_MCMC$mu[paramsAlgo$nbIterSAEM+1,])
theta.estim$sigma2 = res_SAEM_MCMC$sigma2[paramsAlgo$nbIterSAEM+1]
theta.estim$Omega = res_SAEM_MCMC$Omega[paramsAlgo$nbIterSAEM+1,,]

paramsAlgo2 <- paramsAlgo
paramsAlgo2$keep = TRUE
paramsAlgo2$nbIterMCMC<- 5000
paramsAlgo2$print=TRUE
post_Ind_param_MCMC <- MCMC_NLME(res_SAEM_MCMC$Z,myData.sim,theta.estim,paramsAlgo2)
 
seq_iter = seq(paramsAlgo2$nbIterMCMC/2, paramsAlgo2$nbIterMCMC,by=2)
i = sample(1:nPlot,1)
Data.i <- myData.sim%>%filter(plot == i)
seq_time <- seq(min(Data.i$time),max(Data.i$time),len=100)
Z.post.i <- post_Ind_param_MCMC[seq_iter,i,]
F.i <- sapply(1:nrow(Z.post.i),function(m){
    FourParamLogis(seq_time,exp(Z.post.i[m,]))}) 
d1 = dim(F.i)[1]
d2 = dim(F.i)[2]
F.i = F.i + sqrt(theta.estim$sigma2)*matrix(rnorm(d1*d2,0,1),d1,d2)
q.i <-  apply(F.i, 1, quantile, probs = c(0.05, 0.95), na.rm = TRUE)
Data_simpost.i=data.frame(time=seq_time )
Data_simpost.i$q05<- q.i[1,]
Data_simpost.i$q95<- q.i[2,]
  
ggplot() +
    geom_ribbon(data = Data_simpost.i,
                aes(x = time, ymin = q05, ymax = q95),
                fill = "grey80") +
    geom_point(data = Data.i,
               aes(x =  time, y =lfmc ))


 

