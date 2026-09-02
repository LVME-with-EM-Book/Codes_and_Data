# Fit genotype mixture & admixuture for Taita thrush data (GLS00)

rm(list=ls()); palette('R3'); par(pch=20)

# -> Book: TaitaThrush-k5-m8-5-6-3-4-10-8-seed1-try20.Rdata

# Dirs
library(extraDistr); library(tensor); library(MSCquartets); library(sirt); 
library(dirmult)
library(ade4); library(mclust)
source('../Functions/FunctionsGenotypeMixture.R')
source('../Functions/FunctionsLatex.R')
dataDir <- '../../Data/TaitaThrush/'
dataName <- 'TaitaThrush'
figDir <- '../../Figures/'
exportFig <- FALSE

# Data 
# dataName <- 'TaitaThrush-m22-5-6-3-4-10-8'
dataName <- 'TaitaThrush-m8-5-6-3-4-10-8'
load(paste0(dataDir, dataName, '.Rdata'))
n <- nrow(genotypes); p <- ncol(genotypes);

# Algo parms
seed <- 1; set.seed(seed)
kMax <- 5; tryNb <- 20; tol <- 1e-6
algoParms <- paste0('-seed', seed, '-try', tryNb, '-tol', tol)

################################################################################
# Fit admixture (with Dirichlet prob_i=W_i)
resName <- paste0(dataName, '-admixture-k', kMax)
resFile <- paste0(dataDir, resName, algoParms, '.Rdata')
if(!file.exists(resFile)){
  emList <- initList <- list(); par(mfrow=c(ceiling(sqrt(2*kMax)), round(sqrt(2*kMax))))
  for(k in 1:kMax){
    # Init one tau / marker
    cat('\n k =', k, 'each:')
    tauInit <- array(dim=c(n, p, k))
    for(j in 1:p){tauInit[, j, ] <- InitEM(genotypes=genotypes[, j, drop=FALSE], k=k)$tau}
    if(k > 1){for(j in 1:p){
      tauInit[, j, ] <- tauInit[, j, ] + 1/k
      tauInit[, j, ] <- tauInit[, j, ] / rowSums(tauInit[, j, ])
    }}
    emList[[k]] <- VEMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                          k=k, tauInit=tauInit, tol=tol)
    cat(' ->', emList[[k]]$elbo)
    plot(emList[[k]]$elboPath, type='b', ylab='', xlab='', main=paste('k =', k), 
         ylim=quantile(emList[[k]]$elboPath, prob=c(.1, 1)))
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    # Init same tau
    cat('\n k =', k, 'same:')
    tauInit <- array(dim=c(n, p, k))
    tauInit[, 1, ] <- InitEM(genotypes=genotypes[, j, drop=FALSE], k=k)$tau
    for(j in 2:p){tauInit[, j, ] <- tauInit[, 1, ]}
    if(k > 1){for(j in 1:p){
      tauInit[, j, ] <- tauInit[, j, ] + 1/k
      tauInit[, j, ] <- tauInit[, j, ] / rowSums(tauInit[, j, ])
    }}
    emTmp <- VEMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                    k=k, tauInit=tauInit, tol=tol)
    cat(' ->', emTmp$elbo)
    if(emTmp$elbo > emList[[k]]$elbo){cat(' !!!'); emList[[k]] <- emTmp}
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    # Random init one tau / marker
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try each', tt, ':')
      tauTry <- array(dim=c(n, p, k))
      for(j in 1:p){tauTry[, j, ] <- rdirichlet(n, rep(1, k))}
      emTmp <- VEMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                      k=k, tauInit=tauTry, tol=tol)
      cat(' ->', emTmp$elbo)
      if(emTmp$elbo > emList[[k]]$elbo){cat(' !!!'); emList[[k]] <- emTmp}
    }}
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    # Random init same tau / marker
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try same', tt, ':')
      tauTry <- array(dim=c(n, p, k))
      tauTry[, 1, ] <- rdirichlet(n, rep(1, k))
      for(j in 2:p){tauTry[, j, ] <- tauTry[, 1, ]}
      emTmp <- VEMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                      k=k, tauInit=tauTry, tol=tol)
      cat(' ->', emTmp$elbo)
      if(emTmp$elbo > emList[[k]]$elbo){cat(' !!!'); emList[[k]] <- emTmp}
    }}
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    plot(emList[[k]]$elboPath, type='b', ylab='', xlab='', main=paste('k =', k), 
         ylim=quantile(emList[[k]]$elboPath, prob=c(.1, 1)))
  }
  save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
}else{load(resFile)}

################################################################################
# Fit mixtures
resName <- paste0(dataName, '-mixture-k', kMax)
resFile <- paste0(dataDir, resName, algoParms, '.Rdata')
if(!file.exists(resFile)){
  emList <- initList <- list(); par(mfrow=c(ceiling(sqrt(kMax)), round(sqrt(kMax))))
  for(k in 1:kMax){
    cat('\n k =', k, ':')
    # Init tau
    initList[[k]] <- InitEM(genotypes=genotypes, k=k)
    if(k > 1){
      initList[[k]]$tau <- initList[[k]]$tau + 1/k
      initList[[k]]$tau <- initList[[k]]$tau / rowSums(initList[[k]]$tau)
    }
    emList[[k]] <- EMgenotypeMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                       k=k, tauInit=initList[[k]]$tau, tol=tol)
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    apply(emList[[k]]$gamma, c(1, 2), sum, na.rm=TRUE)
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try', tt, ':')
      emTmp <- EMgenotypeMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                   k=k, tauInit=rdirichlet(n, rep(1, k)), tol=tol)
      if(emTmp$logLik > emList[[k]]$logLik){cat('!!!'); emList[[k]] <- emTmp}
    }}
    plot(emList[[k]]$logLikPath, type='b', ylab='', xlab='', main=paste('k =', k))
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
  }
}else{load(resFile)}

################################################################################
# Fit fake admixture (with fixed prob_i)
resName <- paste0(dataName, '-admixtureFake-k', kMax)
resFile <- paste0(dataDir, resName, algoParms, '.Rdata')
if(!file.exists(resFile)){
  emList <- initList <- list(); par(mfrow=c(ceiling(sqrt(2*kMax)), round(sqrt(2*kMax))))
  for(k in 1:kMax){
    # Init one tau / marker
    cat('\n k =', k, 'each:')
    tauInit <- array(dim=c(n, p, k))
    for(j in 1:p){tauInit[, j, ] <- InitEM(genotypes=genotypes[, j, drop=FALSE], k=k)$tau}
    if(k > 1){for(j in 1:p){
      tauInit[, j, ] <- tauInit[, j, ] + 1/k
      tauInit[, j, ] <- tauInit[, j, ] / rowSums(tauInit[, j, ])
    }}
    emList[[k]] <- EMgenotypeFakeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                       k=k, tauInit=tauInit, tol=tol)
    cat(' ->', emList[[k]]$logLik)
    plot(emList[[k]]$logLikPath, type='b', ylab='', xlab='', main=paste('k =', k),
         ylim=quantile(emList[[k]]$logLikPath, prob=c(.1, 1)))
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    # Init same tau
    cat('\n k =', k, 'same:')
    tauInit <- array(dim=c(n, p, k))
    tauInit[, 1, ] <- InitEM(genotypes=genotypes[, j, drop=FALSE], k=k)$tau
    for(j in 2:p){tauInit[, j, ] <- tauInit[, 1, ]}
    if(k > 1){for(j in 1:p){
      tauInit[, j, ] <- tauInit[, j, ] + 1/k
      tauInit[, j, ] <- tauInit[, j, ] / rowSums(tauInit[, j, ])
    }}
    emTmp <- EMgenotypeFakeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                         k=k, tauInit=tauInit, tol=tol)
    cat(' ->', emTmp$logLik)
    if(emTmp$logLik > emList[[k]]$logLik){cat(' !!!'); emList[[k]] <- emTmp}
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    # Random init one tau / marker
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try each', tt, ':')
      tauTry <- array(dim=c(n, p, k))
      for(j in 1:p){tauTry[, j, ] <- rdirichlet(n, rep(1, k))}
      emTmp <- EMgenotypeFakeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                     k=k, tauInit=tauTry, tol=tol)
      cat(' ->', emTmp$logLik)
      if(emTmp$logLik > emList[[k]]$logLik){cat(' !!!'); emList[[k]] <- emTmp}
    }}
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    # Random init same tau / marker
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try same', tt, ':')
      tauTry <- array(dim=c(n, p, k))
      tauTry[, 1, ] <- rdirichlet(n, rep(1, k))
      for(j in 2:p){tauTry[, j, ] <- tauTry[, 1, ]}
      emTmp <- EMgenotypeFakeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                     k=k, tauInit=tauTry, tol=tol)
      cat(' ->', emTmp$logLik)
      if(emTmp$logLik > emList[[k]]$logLik){cat(' !!!'); emList[[k]] <- emTmp}
    }}
    plot(emList[[k]]$logLikPath, type='b', ylab='', xlab='', main=paste('k =', k),
         ylim=quantile(emList[[k]]$logLikPath, prob=c(.1, 1)))
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
  }
}else{load(resFile)}

