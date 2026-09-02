library(ggplot2)
library(parallel)
library(viridis)
library(dplyr)
#----------------------------------------- 
source('Functions_EM.R')
#---------------------------------------- 



#---------------------------------------- 
#-------------- SIMULATION ---------------- 
#---------------------------------------- 

#p = 0.3
#mu = c(2,5)
#sigma = c(1,1) 

#n=1000 
#U = rbinom(n,1,p)
#Y = rnorm(n,mu[1],1)*U + (1-U)*rnorm(n,mu[2],1)
#sampleY <- as.data.frame(Y)
#names(sampleY)='y'
#gg_hist <- ggplot(sampleY,aes(x=y)) + geom_histogram()
#params = list(mu=mu,p=p,sigma=sigma,n=n)

#params$sigma = c(1,1)
#save(sampleY,gg_hist, params,file='datasim_Mixture2.Rdata')    
load(file='datasim_Mixture2.Rdata')  
params_true <- params
rm(params)
y <- sampleY$y
#--------------------------------------------- 
#------------    Figure for likelihood
#-------------------------------------------- 
N = 500
res_like = data.frame(mu1=double(),
                      mu2 = double(), 
                      ll = double())
grid_mu1 = seq(0,7,len=N)
grid_mu2 = seq(0,7,len=N)
grid_params <- as.data.frame(matrix(0,N*N,2))
names(grid_params) <- c('mu1', 'mu2') 
grid_params[,1]<-  rep(grid_mu1,N)
grid_params[,2] <- rep(grid_mu2,each=N)
res_loglik <- mclapply(1:N^2,function(i){
  params.i = params_true
  params.i$mu <- c(grid_params[i,1],grid_params[i,2])
  logLikelihood(params.i,y=sampleY$y)},mc.cores = 4)
res_like = grid_params
res_like[,3]  = unlist(res_loglik)
names(res_like)[3] <- 'loglik'


##########################################################
#------------  EM on only MU 
##########################################################

#---------------- RUNs  EM from various starting points

params0 = list(list(mu=c(1,2),sigma=c(1,1),p=0.3),list(mu=c(0.5,4),sigma=c(1,1),p=0.3),list(mu=c(6,1),sigma=c(1,1),p=0.3),list(mu=c(5,6.5),sigma=c(1,1),p=0.3))
params0[[5]] = list(mu=c(6,6),sigma=c(1,1),p=0.3)
res_all_run <-do.call(rbind,lapply(1:length(params0),function(run){
  res_EM_run  <- EM_2Mixture(params0[[run]],y,tol=10^{-5},estim_only_mu = TRUE)
  niter <- nrow(res_EM_run)
  res_EM_run <- cbind(res_EM_run,rep(run,niter))
  res_EM_run <- cbind(res_EM_run, c('init',rep('iter',nrow(res_EM_run)-1)))
  names(res_EM_run)[7]='Init'
  names(res_EM_run)[8]='numIter'
  
  return(res_EM_run)}))
  


#------------ PLOT ------------------
res_all_run$Init <- as.factor(res_all_run$Init)
res_all_run$numIter <- as.factor(res_all_run$numIter)
true_mu <- as.data.frame(matrix(params_true$mu,nrow=1))
names(true_mu)  =c('mu1','mu2')
gg <- ggplot(res_like,aes(x=mu1,y=mu2,z=loglik)) + geom_tile(aes(fill=loglik)) + scale_fill_viridis() + stat_contour(color="white", size=0.25)
gg <- gg + geom_point(data=res_all_run, aes(x=mu1, y=mu2,z=NULL,colour=Init,shape=numIter),size=3) + geom_line(data=res_all_run, aes(x=mu1, y=mu2,z=NULL,col=Init),linewidth=1.3) + scale_shape_manual(values=c(15, 20))
gg <- gg + geom_point(data=true_mu, aes(x=mu1, y=mu2,z=NULL),size=4,color='black') 
#gg <- gg+  geom_line(linewidth=1.3) + geom_point(size=1.5) 
gg<- gg + theme(axis.text=element_text(size=15),axis.title=element_text(size=15),legend.text = element_text(size=10))
gg <- gg  + theme(legend.title = element_text(size=15))+ guides(shape="none") + labs(x=expression(mu[1]),y=expression(mu[2]))
gg

#scale_fill_hue()$palette(5)  
##########################################################
#------------  EM with MU 
##########################################################


params0 = list(list(mu=c(1,2),sigma=c(1,1),p=0.5),list(mu=c(0.5,4),sigma=c(1,1),p=0.5),list(mu=c(6,1),sigma=c(2,2),p=0.5),list(mu=c(5,6.5),sigma=c(2,2),p=0.5))
params0[[5]] = list(mu=c(4,4),sigma=c(1,1),p=0.5)
res_all_run_all <-do.call(rbind,lapply(1:length(params0),function(run){
  res_EM_run  <- EM_2Mixture(params0[[run]],y,tol=10^{-7},estim_only_mu = FALSE)
  niter <- nrow(res_EM_run)
  res_EM_run <- cbind(res_EM_run,rep(run,niter))
  res_EM_run <- cbind(res_EM_run, c(1:nrow(res_EM_run)))
  names(res_EM_run)[7]='Init'
  names(res_EM_run)[8]='numIter'
  
  return(res_EM_run)}))

#------------------------------------------------
tabl <- matrix(unlist(params0),ncol=5,byrow=TRUE)
tabl <- cbind(tabl,rep(0,5))
row.names(tabl) = 1:5
for (init in 1:5){
  w <- which(res_all_run_all$Init==init)
  tabl[init,6]<- max(res_all_run_all$loglik[w])
}
lltrue <- logLikelihood(params_true[1:5],y)
tabl <- rbind(matrix(c(unlist(params_true)[1:5],lltrue),nrow=1),tabl)
xtable(tabl)


res_all_run_all$Init <- as.factor(res_all_run_all$Init)
res_all_run_all <- res_all_run_all %>% filter(numIter>3)
gg <- ggplot(res_all_run_all,aes(x=numIter,y=loglik,colour = Init)) + geom_line(linewidth=1.3) + geom_point(size=1.5) 
gg<- gg + theme(axis.text=element_text(size=15),axis.title=element_text(size=15),legend.text = element_text(size=15))
gg  + theme(legend.title = element_text(size=15))

