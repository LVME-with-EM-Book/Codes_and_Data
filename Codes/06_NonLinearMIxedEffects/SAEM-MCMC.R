library(mvtnorm)
#####################################################
FourParamLogis = function(time,phi){
  
  #  phi = [A,w,m,s]
  nt = length(time)
  if(is.list(phi)){phi  = unlist(phi)}
  
  #phi = vecteur de taille $4$ ou matrice length(t) lines and $4$ colonnes 
  if(is.vector(phi)){ # n = 1 seul individu
      phi = matrix(phi,nt,4,byrow = TRUE)
  }
  A <- phi[,1]
  w  <- phi[,2]
  m <- phi[,3]
  s <- phi[,4]
  r <- (m-time)/s
  res <- (A-w)/(1+exp(r))+w
  return(res)
}
####################################################### 

completeLogLikNLME <- function(myData,Z,theta){
  
  #myData$lfmc = vecteur de lfmc{ij}
  #myData$time = vecteur de time{ij}
  #myData$indice plot  = vecteur de placette{ij}
  mu <- theta$mu 
  Omega <- theta$Omega
   
  phi <- exp(Z[myData$plot,])
  ly <- sum(dnorm(myData$lfmc,FourParamLogis(myData$time,phi),sqrt(theta$sigma2),log = TRUE))
  lZ <- sum(dmvnorm(Z,theta$mu,theta$Omega,log  = TRUE))
  return(ly+lZ)
  
}

############################################################
#----------------------------- INITIALISATION à l'oeil sur 
############################################################# 

initZ <- function(myData){
  nPlot = max(myData$plot)
  Z.init = matrix(0,nPlot,4)
  for (i in 1:nPlot){
    Data.i <- myData%>%filter(plot == i)
    tmin.i <- min(Data.i$time)  
    tmax.i <- max(Data.i$time)  
    
    A.i = mean(Data.i$lfmc[Data.i$time==tmax.i]) # A.i valeur minimale de lfmc
    w.i = mean(Data.i$lfmc[Data.i$time==tmin.i]) # w.i valeur minimale de lfmc
    
    # Selon le model #log((A.i-w.i)/(Data.i$lfmc-w.i)-1) ~ m.i/s.i - time /s.i
    U.i  =  suppressWarnings(log((A.i-w.i)/(Data.i$lfmc-w.i)-1))
    r.i <- which(!is.finite(U.i) | is.na(U.i)) 
    res_lm.i <- lm(U.i[-r.i] ~ Data.i$time[-r.i])
    s.i <- -1/res_lm.i$coefficients[2]
    m.i <- s.i*res_lm.i$coefficients[1]
    if(m.i<0){m.i = mean(Data.i$time)}
    Z.init[i,] <- log(c(A.i,w.i,m.i,s.i))
  }
  return(Z.init)
}

#-------- if init by optimisation
fitForInit <- function(param,y,time){
  phi= c(param)
  res <- sum((y-FourParamLogis(time,phi))^2)
  return(res)
}

#------------------------- 

initTheta  = function(Z.init,myData,omegaDiag = TRUE){
  ## From Z.init compute mu.init, Omega.init et sigma2.init
  mu.init = colMeans(Z.init)
  if(omegaDiag){
    K = ncol(Z.init)
    Omega.init <- matrix(0,K,K)
    diag(Omega.init) <- apply(Z.init,2,var)
  }else{
    Omega.init = var(Z.init)
  }
  res.init <- myData$lfmc-FourParamLogis(myData$time,exp(Z.init[myData$plot,]))
  sigma2.init <- mean(res.init^2,rm.na = TRUE)
  theta.init  = list(mu=mu.init,Omega=Omega.init,sigma2=sigma2.init)
  return(theta.init)
}


############################################## 
MCMC_NLME = function(Z.init,myData,theta,paramsAlgo){
  
  #------------- params de l'algo 
  nbIterMCMC <-paramsAlgo$nbIterMCMC
  keep <- paramsAlgo$keep
  rho  <- paramsAlgo$rho
  
  #--------------- Initialisation 
  Z <- Z.init 
  nPlot <- nrow(Z.init)
  nParam <- ncol(Z.init)
  ll <- completeLogLikNLME(myData,Z,theta)
  
  
  #------------ stockage ou non de la chaine de markov
  if(keep){
    traceZ <- array(0,c(1+nbIterMCMC,nPlot,nParam))
    traceZ[1,,] <-Z
  }
  
  #---------------------- iterations de
  place  = 1
  for (iter in 1:(nbIterMCMC)){
    if(paramsAlgo$print){print(iter)}
    for (i in 1:nPlot){
        for (k in 1:nParam){
          rho_iter <- sample(rho,1)
          #---------------- 
          Zc = Z 
          Zc[i,k] <- Z[i,k] +  rnorm(1,0,rho_iter)
          #----------------- 
          llc <- completeLogLikNLME(myData,Zc,theta)
          if(llc == -Inf){llc = -10000}
          u = log(runif(1))
          #----------------- 
          if(u < (llc-ll)){
            Z <- Zc
            ll <- llc}
        }# fin k 
    }# fin i 
    if(keep){place= place+1; traceZ[place,,] <-Z}
  }# fin iter
  
  
  #print(keep)
  if(keep){return(traceZ)}else{return(Z)}
}

########################### SAEM MCMC

#paramsAlgo = list(nbIterMCMC=100,nbIterSAEM=2000,rho=c(1/10,1,2),keep = FALSE)
#paramsAlgo$gammaSAEM = c(rep(1,10),1/(1:paramsAlgo$nbIterSAEM))



SAEM_MCMC_NLME <- function(myData,paramsAlgo,omegaDiag = FALSE){
  
  paramsAlgo$keep = FALSE
  if(is.null(paramsAlgo$estim)){
    paramsAlgo$estim=list(mu=c(1,1,1,1),sigma2=1,Omega=1)
  }
  #---------------- params algo 
  nSAEM <-  paramsAlgo$nbIterSAEM
  which_estim_mu <- which(paramsAlgo$estim$mu==1)
  #---------------- Initialisation 
  y <- myData$lfmc
  nPlot <- length(unique(myData$plot))
  n <- length(myData$lfmc)
  
  Z.init <- initZ(myData)
  theta.init <-initTheta(Z.init,myData,omegaDiag)
  
  
  # Stat 0 
  
  S1 <- colSums(Z.init)
  S2  <- t(Z.init) %*% Z.init 
  S3 <- sum((y-FourParamLogis(myData$time,exp(Z.init[myData$plot,])))^2)

  theta = theta.init
  Z = Z.init
  ##-------------- keep trace
  traceMu <- matrix(0,paramsAlgo$nbIterSAEM+1,4)
  traceMu[1,] <- theta$mu
  traceOmega <- array(0,c(paramsAlgo$nbIterSAEM+1,4,4))
  traceOmega[1,,] <- theta$Omega
  traceSigma2 <- rep(0,paramsAlgo$nbIterSAEM)
  traceSigma2[1] <- theta$sigma2
  
  ##------------- 
  for (i in 1:paramsAlgo$nbIterSAEM){
    print(i)
    gamma.i <- paramsAlgo$gammaSAEM[i]
    # Etape Simul 
    Z <- MCMC_NLME(Z,myData,theta,paramsAlgo)
    
    # Etape SA 
    sumZ <- colSums(Z)
    S1 <- S1 + gamma.i*(sumZ-S1)
   
    tZZ =  t(Z) %*% Z
    if(omegaDiag){
      S2 = S2 + gamma.i*(tZZ*diag(1,4)-S2)
    }else{
      S2 <- S2 + gamma.i*(tZZ-S2)
    } 
    
    
    sumRes <- sum((y-FourParamLogis(myData$time,exp(Z[myData$plot,])))^2)
    S3  <- S3 +  gamma.i* (sumRes-S3)
    
    # Etape M 
    theta$mu[which_estim_mu] <- 1/nPlot*S1[which_estim_mu]
    
    if(paramsAlgo$estim$Omega==1){
      if(omegaDiag){
      theta$Omega <- S2/nPlot -  tcrossprod(theta$mu)*diag(1,4)
      }else{
       theta$Omega <- S2/nPlot -  tcrossprod(theta$mu) 
      }
    }
    if(paramsAlgo$estim$sigma2==1){
      theta$sigma2 <-  S3/(n)             
    }
    
    #---------------------------Stockage 
    traceMu[i+1,] <- theta$mu
    traceOmega[i+1,,] <- theta$Omega
    traceSigma2[i+1] <- theta$sigma2  
    
  }
   
  return(list(mu = traceMu,Omega=traceOmega,sigma2 = traceSigma2,Z=Z))
  
}
