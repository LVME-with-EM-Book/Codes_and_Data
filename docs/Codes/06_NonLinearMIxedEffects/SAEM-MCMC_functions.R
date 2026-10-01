library(mvtnorm)

#####################################################
# LFMC response function: declining four-parameter logistic
#####################################################
# Returns  lfmc(t) = (A - w) / (1 + exp((m - t)/s)) + w
# phi = [A, w, m, s], either a vector of length 4 (one individual)
#       or a matrix with length(time) rows and 4 columns.
FourParamLogis = function(time,phi){

  #  phi = [A,w,m,s]
  nt = length(time)
  if(is.list(phi)){phi  = unlist(phi)}

  # phi = vector of size 4, or a matrix of length(t) rows and 4 columns
  if(is.vector(phi)){ # n = 1 individual only
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

#####################################################
# Complete-data log-likelihood of the NLME model
#####################################################
# log p(y, Z ; theta) = log p(y | Z ; theta) + log p(Z ; theta)
#   - myData$lfmc : vector of lfmc{ij} (observations)
#   - myData$time : vector of time{ij}
#   - myData$plot : vector of plot index for each observation
#   - Z           : matrix (nPlot x 4) of individual log-parameters
#   - theta       : list(mu, Omega, sigma2)
# The individual parameters phi_i = exp(Z_i) enter the mean response.
completeLogLikNLME <- function(myData,Z,theta){

  #myData$lfmc = vector of lfmc{ij}
  #myData$time = vector of time{ij}
  #myData$plot index = vector of plot{ij}
  mu <- theta$mu
  Omega <- theta$Omega

  phi <- exp(Z[myData$plot,])
  # Observation log-likelihood: sum over i,j of log N(lfmc{ij}; f(t_ij,phi_i), sigma2)
  ly <- sum(dnorm(myData$lfmc,FourParamLogis(myData$time,phi),sqrt(theta$sigma2),log = TRUE))
  # Random-effects log-density: sum over i of log N(phi_i; mu, Omega)
  lZ <- sum(dmvnorm(Z,theta$mu,theta$Omega,log  = TRUE))
  return(ly+lZ)

}

############################################################
#----------------------------- INITIALISATION  
############################################################# 

initZ <- function(myData){
  nPlot = max(myData$plot)
  Z.init = matrix(0,nPlot,4)
  for (i in 1:nPlot){
    Data.i <- myData%>%filter(plot == i)
    tmin.i <- min(Data.i$time)
    tmax.i <- max(Data.i$time)

    A.i = mean(Data.i$lfmc[Data.i$time==tmax.i]) # A.i minimum value of lfmc
    w.i = mean(Data.i$lfmc[Data.i$time==tmin.i]) # w.i minimum value of lfmc

    # According to the model #log((A.i-w.i)/(Data.i$lfmc-w.i)-1) ~ m.i/s.i - time /s.i
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

# #-------- if init by optimisation
# fitForInit <- function(param,y,time){
#   phi= c(param)
#   res <- sum((y-FourParamLogis(time,phi))^2)
#   return(res)
# }

#-------------------------

initTheta  = function(Z.init,myData,omegaDiag = TRUE){
  ## From Z.init compute mu.init, Omega.init and sigma2.init
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
# Metropolis-within-Gibbs sampler of the individual parameters
##############################################
# Samples (approximately) from the posterior p(Z | y ; theta) for fixed theta.
# Each coordinate Z[i,k] is updated in turn with a random-walk Metropolis step:
#   proposal  Zc[i,k] = Z[i,k] + N(0, rho)   (rho drawn at random in paramsAlgo$rho)
#   acceptance with probability exp(min(0, loglik(Zc) - loglik(Z)))
# If keep = TRUE, the whole chain traceZ (nbIterMCMC+1) is returned,
# otherwise only the last state Z is returned.
MCMC_NLME = function(Z.init,myData,theta,paramsAlgo){

  #------------- algorithm parameters
  nbIterMCMC <-paramsAlgo$nbIterMCMC
  keep <- paramsAlgo$keep
  rho  <- paramsAlgo$rho

  #--------------- Initialisation of the chain
  Z <- Z.init
  nPlot <- nrow(Z.init)
  nParam <- ncol(Z.init)
  ll <- completeLogLikNLME(myData,Z,theta)


  #------------ store the whole chain or not
  if(keep){
    traceZ <- array(0,c(1+nbIterMCMC,nPlot,nParam))
    traceZ[1,,] <-Z
  }

  #---------------------- MCMC iterations
  place  = 1
  for (iter in 1:(nbIterMCMC)){
    if(paramsAlgo$print){print(iter)}
    for (i in 1:nPlot){
        for (k in 1:nParam){
          # random-walk proposal scale
          rho_iter <- sample(rho,1)
          #----------------
          Zc = Z
          Zc[i,k] <- Z[i,k] +  rnorm(1,0,rho_iter)
          #-----------------
          # Metropolis acceptance ratio (in log scale)
          llc <- completeLogLikNLME(myData,Zc,theta)
          if(llc == -Inf){llc = -10000}
          u = log(runif(1))
          #-----------------
          if(u < (llc-ll)){
            Z <- Zc
            ll <- llc}
        }# end k (parameter)
    }# end i (plot)
    if(keep){place= place+1; traceZ[place,,] <-Z}
  }# end iter (MCMC iterations)


  #print(keep)
  if(keep){return(traceZ)}else{return(Z)}
}

###########################
# SAEM-MCMC: Stochastic Approximation EM with an MCMC step
###########################
# Estimates theta = (mu, Omega, sigma2) by maximum likelihood, where the
# latent individual parameters Z are sampled with an MCMC step (instead of an
# exact simulation). At iteration i, with step size gamma_i:
#   Simulation (E-step):  sample Z ~ p(Z | y ; theta) via MCMC_NLME
#   Stochastic approx. :  update the sufficient statistics S1, S2, S3
#                         S <- S + gamma_i (s(Z) - S)
#   Maximisation (M)   :  update theta from the statistics
#
# paramsAlgo = list(nbIterMCMC, nbIterSAEM, rho, keep, print, gammaSAEM)
#   gammaSAEM is the vector of step sizes gamma_i.
# omegaDiag = TRUE forces Omega to be diagonal.
#
# e.g. paramsAlgo = list(nbIterMCMC=100, nbIterSAEM=2000, rho=c(1/10,1,2), keep = FALSE)
#      paramsAlgo$gammaSAEM = c(rep(1,10), 1/(1:paramsAlgo$nbIterSAEM))

SAEM_MCMC_NLME <- function(myData,paramsAlgo,omegaDiag = FALSE){

  paramsAlgo$keep = FALSE
  # By default all parameters are estimated (subset estimation is possible
  # by setting paramsAlgo$estim$mu / $Omega / $sigma2 to 0 for fixed ones)
  if(is.null(paramsAlgo$estim)){
    paramsAlgo$estim=list(mu=c(1,1,1,1),sigma2=1,Omega=1)
  }
  #---------------- algorithm parameters
  nSAEM <-  paramsAlgo$nbIterSAEM
  which_estim_mu <- which(paramsAlgo$estim$mu==1)   # mu components to estimate
  #---------------- Initialisation from the optimisation-based init
  y <- myData$lfmc
  nPlot <- length(unique(myData$plot))
  n <- length(myData$lfmc)

  Z.init <- initZ(myData)
  theta.init <-initTheta(Z.init,myData,omegaDiag)


  # Initial sufficient statistics (at Z = Z.init)
  #   S1 = sum_i Z_i            (sufficient stat for mu)
  #   S2 = sum_i Z_i Z_i'       (sufficient stat for Omega)
  #   S3 = sum_ij residual^2    (sufficient stat for sigma2)
  S1 <- colSums(Z.init)
  S2  <- t(Z.init) %*% Z.init
  S3 <- sum((y-FourParamLogis(myData$time,exp(Z.init[myData$plot,])))^2)

  theta = theta.init
  Z = Z.init
  ##-------------- keep a trace of the estimates across iterations
  traceMu <- matrix(0,paramsAlgo$nbIterSAEM+1,4)
  traceMu[1,] <- theta$mu
  traceOmega <- array(0,c(paramsAlgo$nbIterSAEM+1,4,4))
  traceOmega[1,,] <- theta$Omega
  traceSigma2 <- rep(0,paramsAlgo$nbIterSAEM)
  traceSigma2[1] <- theta$sigma2

  ##------------- main SAEM-MCMC loop
  for (i in 1:paramsAlgo$nbIterSAEM){
    print(i)
    gamma.i <- paramsAlgo$gammaSAEM[i]      # current step size
    # Simulation step: sample the latent individual parameters with MCMC
    Z <- MCMC_NLME(Z,myData,theta,paramsAlgo)

    # Stochastic approximation of the sufficient statistics
    sumZ <- colSums(Z)
    S1 <- S1 + gamma.i*(sumZ-S1)

    tZZ =  t(Z) %*% Z
    if(omegaDiag){
      S2 = S2 + gamma.i*(tZZ*diag(1,4)-S2)   # keep Omega diagonal
    }else{
      S2 <- S2 + gamma.i*(tZZ-S2)
    }


    sumRes <- sum((y-FourParamLogis(myData$time,exp(Z[myData$plot,])))^2)
    S3  <- S3 +  gamma.i* (sumRes-S3)

    # Maximisation step: update theta from the sufficient statistics
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

    #--------------------------- store the current estimates
    traceMu[i+1,] <- theta$mu
    traceOmega[i+1,,] <- theta$Omega
    traceSigma2[i+1] <- theta$sigma2

  }

  return(list(mu = traceMu,Omega=traceOmega,sigma2 = traceSigma2,Z=Z))

}
