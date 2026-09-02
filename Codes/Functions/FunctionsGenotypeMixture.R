LogDirichlet <-function(alpha){sum(lgamma(alpha))-lgamma(sum(alpha))} 

EntropyDirichlet <-function(alpha){LogDirichlet(alpha) + 
    (sum(alpha)-length(alpha))*digamma(sum(alpha)) -
    sum((alpha-1)*digamma(alpha))} 

Geno2Genotype <- function(geno){
  if(sum(is.na(geno)) > 0){geno <- geno[-which(rowSums(is.na(geno))>0), ]}
  p2 <- ncol(geno); genotypes <- matrix('', nrow(geno), p2/2)
  miss <- geno <- matrix(as.numeric(geno), nrow=nrow(geno), ncol=p2)
  for(i in 1:nrow(geno)){for(j in 1:(p2/2)){
    tmp <- sort(c(geno[i, 2*(j-1)+1], geno[i, 2*(j-1)+2]))
    if(tmp[1]==0){miss[i, j] <- 1}else{miss[i, j] <- 0}
    genotypes[i, j] <- paste0(tmp[1], tmp[2]) 
  }}
  return(list(genotypes=genotypes, miss=miss))
}

MakeGenoArray <- function(markers){
  n <- nrow(markers); p <- ncol(markers)
  alleles <- lapply(1:p, function(j){sort(unique(as.factor(markers[, j])))})
  m <- sapply(1:p, function(j){length(unique(as.factor(markers[, j])))})
  mMax <- max(m)
  markArray <- array(0, dim=c(n, p, mMax))
  for(i in 1:n){for(j in 1:p){
    markArray[i, j, which(alleles[[j]]==markers[i, j])] <- 1
    }}
  return(markArray)
}

InitEM <- function(genotypes, k, tauEps=1e-4){
  n <- nrow(genotypes); p <- ncol(genotypes)
  genoFrame <- list()
  for(j in 1:p){genoFrame[[j]] <- as.factor(genotypes[, j])}
  mca <- dudi.acm(genoFrame, scannf=FALSE, nf=10)
  cluster <- kmeans(mca$li, centers=k, nstart=100, iter.max=100)$cluster
  tau <- matrix(tauEps, n, k)
  for(i in 1:n){tau[i, cluster[i]] <- tau[i, cluster[i]] + 1}
  tau <- tau / rowSums(tau)
  return(list(mca=mca, cluster=cluster, tau=tau))
}

LogPhiHW <- function(haploTab, logGammaHaplo){
  sum(haploTab[1, , ]*logGammaHaplo, na.rm=TRUE) +
    sum(haploTab[2, , ]*logGammaHaplo, na.rm=TRUE) +
    log(2)*sum(haploTab[1, , ]*haploTab[2, , ], na.rm=TRUE)
  }

DirLogLik <- function(alpha, x, p=mean(rowSums(x))){
  k <- length(alpha); n <- nrow(x); alphat <- rep(1, n)%o%alpha + x; 
  return(- n*LogDirichlet(alpha) + sum(log(x)%*%(alpha-1)))
}

DirGradLogLik <- function(alpha, x, p=mean(rowSums(x))){
  k <- length(alpha); n <- nrow(x); alphat <- rep(1, n)%o%alpha + x; 
  return(n*rep(digamma(sum(alpha)), k) - n*digamma(alpha) + colSums(log(x)))
}

DirMLE <- function(alpha, x){
  fit <- optim(par=alpha, fn=DirLogLik, gr=DirGradLogLik, x=x, 
               control=list(fnscale=-1), method='L-BFGS-B', lower=rep(1e-5, k))
  return(list(alpha=fit$par, logLik=fit$value))
}

DirMultLogLik <- function(alpha, x, p=mean(rowSums(x))){
  k <- length(alpha); n <- nrow(x); alphat <- rep(1, n)%o%alpha + x; 
  return(- n*lgamma(sum(alpha)+p) + n*lgamma(sum(alpha)) + n*lgamma(p+1) + 
           sum(lgamma(alphat)) - n*sum(lgamma(alpha)) - sum(lgamma(x+1)))
}

DirMultGradLogLik <- function(alpha, x, p=mean(rowSums(x))){
  k <- length(alpha); n <- nrow(x); alphat <- rep(1, n)%o%alpha + x; 
  return(- n*rep(digamma(sum(alpha)+p), k) + n*rep(digamma(sum(alpha)), k) + 
           colSums(digamma(alphat)) - n*digamma(alpha))
}

DirMultMLE <- function(alpha, x){
  fit <- optim(par=alpha, fn=DirMultLogLik, gr=DirMultGradLogLik, x=x, 
               control=list(fnscale=-1), method='L-BFGS-B', lower=rep(1e-5, k))
  return(list(alpha=fit$par, logLik=fit$value))
}

EMgenotypeMixtureHW <- function(haplotypes, haploArray, k=ncol(tauInit), tauInit, tol=1e-6, iterMax=1e3, tauEps=1e-4){
  # Genotype mixture under Hardy-Weinberg equilibrium
  # k=3; tauInit=initList[[k]]$tau; tol=1e-6; iterMax=1e1; tauEps=1e-4
  n <- nrow(haplotypes)/2; p <- ncol(haplotypes)
  geno2haplo <- matrix(0, n, 2*n)
  for(i in 1:n){geno2haplo[i, 2*(i-1)+1] <- geno2haplo[i, 2*(i-1)+2] <- 1}
  m <- sapply(1:p, function(j){length(unique(haplotypes[, j]))})
  mMax <- max(m)
  tauGeno <- tauInit; 
  diff <- 2*tol; iter <- 0
  logLik <- rep(NA, iterMax)
  while((diff > tol) & (iter < iterMax)){
    iter <- iter+1; if(iter%%round(sqrt(iterMax))==0){cat('', iter)}
    # M step
    nGroup <- colSums(tauGeno)
    pi <- nGroup/n
    logPi <- log(pi)
    tauHaplo <- t(geno2haplo) %*% tauGeno  
    gamma <- array(dim=c(k, p, mMax))
    for(g in 1:k){for(j in 1:p){for(a in 1:m[j]){
      gamma[g, j, a] <- sum(haploArray[, j, a]*tauHaplo[, g])/(sum((!is.na(haplotypes[, j]))*tauHaplo[, g]))
    }}}
    logGamma <- log(gamma)
    logPhi <- matrix(0, n, k)
    for(i in 1:n){for(g in 1:k){
      logPhi[i, g] <- LogPhiHW(haploArray[2*(i-1)+(1:2), , ], logGamma[g, , ])}}
    
    # E step
    logTau <- logPhi + rep(1, n)%o%logPi
    logTau <- logTau - apply(logTau, 1, max)
    tauNew <- exp(logTau); tauNew <- tauNew / rowSums(tauNew)
    tauNew <- tauNew + tauEps; tauNew <- tauNew / rowSums(tauNew)
    
    # Test
    diff <- max(abs(tauNew - tauGeno))
    tauGeno <- tauNew
    logLik[iter] <- sum(tauGeno%*%logPi) + sum(tauGeno*logPhi) - sum(tauGeno*log(tauGeno))
    # cat(logLik[iter], diff, '\n')
  }
  return(list(tau=tauGeno, logPhi=logPhi, gamma=gamma, pi=pi, m=m, 
              iter=iter, logLik=logLik[iter], logLikPath=logLik[1:iter]))
}

EMgenotypeFakeAdMixtureHW <- function(haplotypes, haploArray, k=ncol(tauInit), tauInit, tol=1e-6, iterMax=1e3, tauEps=1e-4){
  # Genotype mixture under Hardy-Weinberg equilibrium
  # tol=1e-6; iterMax=1e1; tauEps=1e-4
  n <- nrow(haplotypes)/2; p <- ncol(haplotypes)
  geno2haplo <- matrix(0, n, 2*n)
  for(i in 1:n){geno2haplo[i, 2*(i-1)+1] <- geno2haplo[i, 2*(i-1)+2] <- 1}
  m <- sapply(1:p, function(j){length(unique(haplotypes[, j]))})
  mMax <- max(m)
  tauGeno <- tauInit; 
  diff <- 2*tol; iter <- 0
  logLik <- rep(NA, iterMax)
  while((diff > tol) & (iter < iterMax)){
    iter <- iter+1; if(iter%%round(sqrt(iterMax))==0){cat('', iter)}
    # M step
    pi <- apply(tauGeno, c(1, 3), sum)
    pi <- pi/rowSums(pi)
    logPi <- log(pi)
    tauHaplo <- tensor(t(geno2haplo), tauGeno, alongA=2, alongB=1)  
    nGroup <- apply(tauHaplo, c(2, 3), sum)
    gamma <- array(dim=c(k, p, mMax))
    for(g in 1:k){for(j in 1:p){for(a in 1:m[j]){
      gamma[g, j, a] <- sum(haploArray[, j, a]*tauHaplo[, j, g])/(sum((!is.na(haplotypes[, j]))*tauHaplo[, j, g]))
    }}}
    # for(g in 1:k){for(j in 1:p){for(a in 1:m[j]){
    #   gamma[g, j, a] <- sum(haploArray[, j, a]*tauHaplo[, j, g])/nGroup[j, g]
    # }}}
    logGamma <- log(gamma)
    logPhi <- array(dim=c(n, p, k))
    for(i in 1:n){for(j in 1:p){for(g in 1:k){
      logPhi[i, j, g] <- LogPhiHW(haploArray[2*(i-1)+(1:2), j, , drop=FALSE], logGamma[g, j, , drop=FALSE])
    }}}
    
    # E step
    tauNew <- logPhi
    for(j in 1:p){
      tauNew[, j, ] <- tauNew[, j, ] + logPi
      tauNew[, j, ] <- tauNew[, j, ] - apply(tauNew[, j, , drop=FALSE], 1, max)
      tauNew[, j, ] <- exp(tauNew[, j, ])
      tauNew[, j, ] <- tauNew[, j, ] / rowSums(tauNew[, j,  , drop=FALSE])
      tauNew[, j, ] <- tauNew[, j, ] + tauEps
      tauNew[, j, ] <- tauNew[, j, ] / rowSums(tauNew[, j,  , drop=FALSE])
    }

    # Test
    diff <- max(abs(tauNew - tauGeno))
    tauGeno <- tauNew
    logLik[iter] <- sum(sapply(1:p, function(j){tauGeno[, j, ]*logPi})) + 
                          sum(tauGeno*logPhi) - sum(tauGeno*log(tauGeno))
    # cat(logLik[iter], diff, '\n')
  }
  return(list(tau=tauGeno, logPhi=logPhi, gamma=gamma, pi=pi, m=m, 
              iter=iter, logLik=logLik[iter], logLikPath=logLik[1:iter]))
}

VEMgenotypeAdMixtureHW <- function(haplotypes, haploArray, k=ncol(tauInit), tauInit, tol=1e-4, iterMax=1e3, tauEps=1e-4){
  # Genotype admixture under Hardy-Weinberg equilibrium
  # k=3; tauInit=emList[[k]]$tau; tol=1e-6; iterMax=1e3; tauEps=1e-4
  n <- nrow(haplotypes)/2; p <- ncol(haplotypes)
  geno2haplo <- matrix(0, n, 2*n)
  for(i in 1:n){geno2haplo[i, 2*(i-1)+1] <- geno2haplo[i, 2*(i-1)+2] <- 1}
  m <- sapply(1:p, function(j){length(unique(haplotypes[, j]))})
  mMax <- max(m)
  tauGeno <- tauInit;
  alpha <- rep(1, k)
  diff <- 2*tol; iter <- 0
  elbo <- rep(NA, iterMax)
  while((diff > tol) & (iter < iterMax)){
    iter <- iter+1; if(iter%%round(sqrt(iterMax))==0){cat('', iter)}
    # M step
    if(iter > 1){alpha <- DirMLE(alpha=alpha, x=exp(espLogW))$alpha}
    tauHaplo <- tensor(t(geno2haplo), tauGeno, alongA=2, alongB=1)  
    nGroup <- apply(tauHaplo, c(2, 3), sum)
    gamma <- array(dim=c(k, p, mMax))
    for(g in 1:k){for(j in 1:p){for(a in 1:m[j]){
      gamma[g, j, a] <- sum(haploArray[, j, a]*tauHaplo[, j, g])/(sum((!is.na(haplotypes[, j]))*tauHaplo[, j, g]))
    }}}
    logGamma <- log(gamma)
    logPhi <- array(dim=c(n, p, k))
    for(i in 1:n){for(j in 1:p){for(g in 1:k){
      logPhi[i, j, g] <- LogPhiHW(haploArray[2*(i-1)+(1:2), j, , drop=FALSE], logGamma[g, j, , drop=FALSE])
    }}}
    
    # VE step
    alphaTilde <- rep(1, n)%o%alpha + apply(tauGeno, c(1, 3), sum)
    espLogW <- digamma(alphaTilde) - digamma(rowSums(alphaTilde))%o%rep(1, k)
    tauNew <- array(dim=c(n, p, k)); 
    for(j in 1:p){
      tauNew[, j, ] <- espLogW + logPhi[, j, ] 
      tauNew[, j, ] <- tauNew[, j, ] - apply(tauNew[, j, , drop=FALSE], 1, max)
      tauNew[, j, ] <- exp(tauNew[, j, ])
      tauNew[, j, ] <- tauNew[, j, ] / rowSums(tauNew[, j,  , drop=FALSE])
      tauNew[, j, ] <- tauNew[, j, ] + tauEps
      tauNew[, j, ] <- tauNew[, j, ] / rowSums(tauNew[, j,  , drop=FALSE])
    }
    
    # Test
    diff <- max(abs(tauNew - tauGeno))
    tauGeno <- tauNew
    elbo[iter] <- -n*LogDirichlet(alpha) + 
      sum(espLogW %*% (alpha-1)) +
      sum(sapply(1:p, function(j){tauGeno[, j, ]*espLogW})) + 
      sum(tauGeno*logPhi) +
      sum(sapply(1:n, function(i){EntropyDirichlet(alphaTilde[i, ])})) - 
      sum(tauGeno*log(tauGeno))
    # cat(iter, elbo[iter], diff, '\n')
  }
  # plot(elbo[1:iter], type='b')
  return(list(tau=tauGeno, logPhi=logPhi, gamma=gamma, alpha=alpha, 
              alphaTilde=alphaTilde, espLogW=espLogW, 
              m=m, iter=iter, elbo=elbo[iter], elboPath=elbo[1:iter]))
}
