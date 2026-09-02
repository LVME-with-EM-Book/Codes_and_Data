#-------------------------------------------------------------------------------
# PPCA 
#-------------------------------------------------------------------------------
EMpPCA <- function(Y, q=2, init=NULL, tol=1e-5){
  n <- nrow(Y); p <- ncol(Y)
  S <- cov(Y)*(n-1)/n; eigS <- eigen(S)
  if(q==p){
    B <- eigS$vectors; sigma2 <- 0; Sigma <- S
  }else{
    if(is.null(init)){
      nppca <- NaivePPCA(Y, q)
      init <- list(B=nppca$B, sigma2=nppca$sigma2)
    }
    B <- init$B; sigma2 <- init$sigma2
    diff <- 2*tol
    while(diff > tol){
      Gamma <- t(B) %*% B + sigma2*diag(q)
      invGamma <- solve(Gamma)
      Bnew <- S %*% B %*% solve (sigma2*diag(q) + invGamma%*%t(B)%*%S%*%B)
      sigma2new <- mean(diag(S - S%*%Bnew%*%invGamma%*%t(Bnew)))
      diff <- max(max(abs(Bnew - B)), abs(sigma2 - sigma2new))
      B <- Bnew; sigma2 <- sigma2new
    }
    Sigma <- B%*%t(B) + sigma2*diag(p)
  }
  logL <- sum(dmvnorm(Y, mean=colMeans(Y), sigma=Sigma, log=TRUE))
  parmNb <- 1 + q*(2*p - q + 1)/2 + p # SR 14fev25 : ajout de p suite Rq PG
  if(scale){parmNb <- parmNb - 1}
  aic <- logL - parmNb
  bic <- logL - parmNb*log(n)/2
  return(list(mu=colMeans(Y), B=B, sigma2=sigma2, Sigma=Sigma, Gamma=Gamma, mu=colMeans(Y), 
              logL=logL, parmNb=parmNb, bic=bic, aic=aic, init=init))
}

NaivePPCA <- function(Y, q=2, scale=FALSE){
  Y <- as.matrix(Y)
  n <- nrow(Y); p <- ncol(Y)
  mu <- colMeans(Y); 
  if(scale){
    Ysd <- sqrt(apply(Y, 2, var)*(n-1)/n)
    Y <- (Y - rep(1, n)%o%mu) %*% diag(1/Ysd)
  }
  S <- cov(Y)*(n-1)/n; eigS <- eigen(S)
  if(q < p){sigma2 <- mean(eigS$values[q+(1:(p-q))])}else{sigma2 <- 0}
  if(q==1){
    B <- matrix(eigS$vectors[, 1] * sqrt(eigS$values[1] - sigma2), p, 1)
  }else{
    B <- eigS$vectors[, 1:q] %*% sqrt(diag(eigS$values[1:q] - sigma2))
  }
  Gamma <- t(B) %*% B # + sigma2*diag(q)
  Sigma <- B%*%t(B) + sigma2*diag(p)
  logL <- sum(dmvnorm(Y, mean=colMeans(Y), sigma=Sigma, log=TRUE))
  parmNb <- 1 + q*(2*p - q + 1)/2 + p # SR 14fev25 : ajout de p suite Rq PG
  if(scale){parmNb <- parmNb - 1}
  aic <- logL - parmNb
  bic <- logL - parmNb*log(n)/2
  return(list(mu=mu, B=B, sigma2=sigma2, Sigma=Sigma, Gamma=Gamma, mu=mu, 
              logL=logL, parmNb=parmNb, bic=bic, aic=aic))
}

#-------------------------------------------------------------------------------
# PPCA + missing
#-------------------------------------------------------------------------------
EMpPCAmiss <- function(Yobs, R, q=2, tol=1e-5, iterMax=1e4, init=NULL){
  # tol=1e-5; iterMax=500
  n <- nrow(Yobs); p <- ncol(Yobs); nObs <- sum(R)
  # Missing structure
  Obs <- lapply(1:n, function(i){which(R[i, ]==1)})
  A <- lapply(1:p, function(j){which(R[, j]==1)})
  # Init
  if(is.null(init)){
    Ymean <- rep(1, n) %o% colMeans(Yobs, na.rm=TRUE)
    Ytmp <- Yobs; Ytmp[which(R==0)] <- Ymean[which(R==0)]
    init <- NaivePPCA(Y=Ytmp, q=q, scale=FALSE)
  }
  B <- init$B; sigma2 <- init$sigma2; Sigma <- init$Sigma; mu <- init$mu
  # EM loop
  logLpath <- rep(NA, iterMax); diff <- 2*tol; iter <- 0
  while((iter < iterMax) & (diff > tol)){
    iter <- iter+1
    # E-step
    if(q==1){
      M <- matrix(sapply(1:n, function(i){
        t(B[Obs[[i]], ]) %*% solve(Sigma[Obs[[i]], Obs[[i]]]) %*% (Yobs[i, Obs[[i]]] - mu[Obs[[i]]])
      }), n, 1)
    }else{
      M <- t(sapply(1:n, function(i){
        t(B[Obs[[i]], ]) %*% solve(Sigma[Obs[[i]], Obs[[i]]]) %*% (Yobs[i, Obs[[i]]] - mu[Obs[[i]]])
      }))
    }
    Q <- array(dim=c(n, q, q))
    for(i in 1:n){
      Q[i, , ] <- diag(rep(1, q)) - t(B[Obs[[i]], ]) %*% solve(Sigma[Obs[[i]], Obs[[i]]]) %*% B[Obs[[i]], ]
    }
    # M-step
    muNew <- colMeans(Yobs - M%*%t(B), na.rm=TRUE)
    if(q==1){
      sigma2new <- sum(sapply(1:n, function(i){sum((Yobs[i, Obs[[i]]] - muNew[Obs[[i]]] - B[Obs[[i]], ]%*%M[i, , drop=FALSE])^2) + 
          sum(diag(B[Obs[[i]], ]%*%(Q[i, , ]%*%t(B[Obs[[i]], ]))))})) / nObs
    }else{
      sigma2new <- sum(sapply(1:n, function(i){sum((Yobs[i, Obs[[i]]] - muNew[Obs[[i]]] - B[Obs[[i]], ]%*%M[i, ])^2) + 
          sum(diag(B[Obs[[i]], ]%*%(Q[i, , ]%*%t(B[Obs[[i]], ]))))})) / nObs
    }
    Bnew <- t(sapply(1:p, function(j){
      invQ <- matrix(0, q, q)
      for(i in A[[j]]){invQ <- invQ + (M[i, ]%o%M[i, ]) + Q[i, , ]}
      invQ <- solve(invQ)
      if(q==1){
        return(sum(sapply(A[[j]], function(i){(Yobs[i, j] - muNew[j])*M[i, ]})) %*% invQ)
      }else{
        return(rowSums(sapply(A[[j]], function(i){(Yobs[i, j] - muNew[j])*M[i, ]})) %*% invQ)
      }
    }))
    if(q==1){Bnew <- t(Bnew)}
    # Test & update
    diff <- max(max(abs(Bnew - B)), abs(sigma2 - sigma2new), abs(mu - muNew))
    B <- Bnew; sigma2 <- sigma2new; mu <- muNew
    if(q==1){
      Sigma <- as.vector(B)%o%as.vector(B) + sigma2*diag(p)
    }else{Sigma <- B%*%t(B) + sigma2*diag(p)}
    logLpath[iter] <- sum(sapply(1:n, function(i){
      dmvnorm(Yobs[i, Obs[[i]], drop=FALSE], mean=mu[Obs[[i]]], sigma=Sigma[Obs[[i]], Obs[[i]]], log=TRUE)
    }))
    # print(c(iter, diff, logLpath[iter])); plot(logLpath, type='b', pch=20)
  }
  logLpath[1] <- NA
  logL <- logLpath[iter]; logLpath <- logLpath[1:iter]
  parmNb <- 1 + q*(2*p - q + 1)/2 + p # SR 14fev25 : ajout de p suite Rq PG
  if(scale){parmNb <- parmNb - 1}
  aic <- logL - parmNb
  bic <- logL - parmNb*log(n)/2
  return(list(mu=mu, B=B, sigma2=sigma2, Sigma=Sigma, Gamma=Gamma, logLpath=logLpath, 
              logL=logL, parmNb=parmNb, bic=bic, aic=aic, init=init))
}

#-------------------------------------------------------------------------------
# Utils
#-------------------------------------------------------------------------------
CondMomentsPPCA <- function(Y, ppca){
  # Conditional moments of the latent variables in probabilistic PCA
  # In absence of missing data
  # ppca is a list resulting from EMpPCA or NaivePPCA
  n <- nrow(Y); p <- ncol(Y); q <- ncol(ppca$B); 
  mu <- colMeans(Y); invSigma <- solve(ppca$Sigma)
  M <- (Y - rep(1, n)%o%mu) %*% invSigma %*% ppca$B
  Q <- diag(q) - t(ppca$B) %*% invSigma %*% ppca$B
  return(list(mu=mu, M=M, Q=Q))
}

ImputePPCAall <- function(Yobs, R, ppca){
  # Imputation of all missing values at once
  Obs <- lapply(1:n, function(i){which(R[i, ]==1)})
  Miss <- lapply(1:n, function(i){which(R[i, ]==0)})
  Yimp <- matrix(NA, n, p); Yvar <- array(dim=c(n, p, p))
  for(i in 1:n){
    if(length(Miss[[i]]) > 0){
      invSigmai <- solve(ppca$Sigma[Obs[[i]], Obs[[i]]])
      Yimp[i, Miss[[i]]] <- ppca$mu[Miss[[i]]] + 
        (ppca$Sigma[Miss[[i]], Obs[[i]]] %*%  invSigmai) %*% (Yobs[i, Obs[[i]]] - ppca$mu[Obs[[i]]])
      Yvar[i, Miss[[i]], Miss[[i]]] <- ppca$Sigma[Miss[[i]], Miss[[i]]] - 
        ppca$Sigma[Miss[[i]], Obs[[i]]] %*% invSigmai %*% ppca$Sigma[Obs[[i]], Miss[[i]]]
    }
  }
  return(list(Yimp=Yimp, Yvar=Yvar))
}

ImputePPCAeach <- function(Yfull, R, ppca){
  # Imputation of each missing values one by one
  n <- nrow(Yfull); p <- ncol(Yfull); nObs <- sum(R)
  missList <- c()
  for(i in 1:n){for(j in 1:p){if(R[i, j]==0){missList <- rbind(missList, c(i, j))}}}
  invSigmai <- 1/diag(ppca$Sigma)
  Yimp <- Yvar <- matrix(NA, n, p)
  for(i in 1:n){for(j in 1:p){if(R[i, j]==0){
    invSigmai <- solve(ppca$Sigma[-j, -j])
    Yimp[i, j]<- ppca$mu[j] + 
      (ppca$Sigma[j, -j] %*%  invSigmai) %*% (Yfull[i, -j] - ppca$mu[-j])
    Yvar[i, j] <- ppca$Sigma[j, j] - ppca$Sigma[j, -j] %*% invSigmai %*% ppca$Sigma[-j, j]
  }}}
  return(list(Yimp=Yimp, Yvar=Yvar))
}
