# Loglik and gradients
LogLikBernoulli <- function(alpha, X, Y, W=rep(1, length(Y))){sum(W * (Y*(X%*%alpha) - log(1+exp(X%*%alpha))))}
# LogLikBernoulli <- function(alpha, X, Y, W=rep(1, length(Y))){sum(W * dbinom(Y, rep(1, length(Y)), plogis(X%*%alpha), log=TRUE))}
DerLogLikBernoulli <- function(alpha, X, Y, W=rep(1, length(Y))){as.vector(t(X) %*% (W*(Y - plogis(X%*%alpha))))}
LogLikPoisson <- function(beta, X, Y, W=rep(1, length(Y))){sum(W * dpois(Y, exp(X%*%beta), log=TRUE))}
DerLogLikPoisson <- function(beta, X, Y, W=rep(1, length(Y))){as.vector(t(X) %*% (W*(Y - exp(X%*%beta))))}
LogLikZIP <- function(theta, X, Y, W=rep(1, length(Y))){
  # Bernoulli part = presence
  alpha <- theta[1:round(length(theta)/2)]; beta <- theta[-(1:round(length(theta)/2))]
  pi <- plogis(X%*%alpha); lambda <- exp(X%*%beta)
  return(sum(W * (log((1-pi)*(Y==0) + pi*dpois(Y, lambda)))))
}

# Init
InitZIP <- function(X, Y){
  alpha <- as.vector(glm(1*(Y>0) ~ -1 + X, family='binomial')$coef)
  beta <- as.vector(glm(Y ~ -1 + X, family='poisson')$coef)
  return(list(alpha=alpha, beta=beta))
}

# EM
EMZIPreg <- function(X, Y, thetaInit=NULL, iterMax=1e3, tol=1e-6, tolTau=1e-4){
  # Bernoulli part = presence
  # X=Covariates; Y=Abundance; init=NULL; tolTau=1e-4; tol=1e-6; iterMax=1e3
  # X=matrix(1, length(Abundance), 1); Y=Abundance; thetaInit=NULL; tolTau=1e-4; tol=1e-6; iterMax=1e3
  Y0 <- 1*(Y==0); Yp <- 1 - Y0
  if(is.null(thetaInit)){thetaInit <- InitZIP(X, Y)}
  alpha <- thetaInit$alpha; beta <- thetaInit$beta
  logL <- rep(NA, iterMax); diff <- 2*tol; iter <- 1
  logL[iter] <- LogLikZIP(theta=c(alpha, beta), X=X, Y=Y)
  while((diff > tol) & (iter < iterMax)){
    iter <- iter+1; cat('iter', iter, '\n')
    pi <- plogis(X%*%alpha); lambda <- exp(X%*%beta)
    tau <- Yp + Y0*pi*exp(-lambda)/((1-pi) + pi*exp(-lambda))
    tau <- cbind((1-tau), tau); 
    tau <- tau+tolTau; tau <- tau / rowSums(tau)
    alphaNew <- optim(par=alpha, f=LogLikBernoulli, g=DerLogLikBernoulli, Y=tau[, 2], X=X, W=rep(1, length(Y)), 
                      control=list(fnscale=-1))$par
    betaNew <- optim(par=beta, f=LogLikPoisson, g=DerLogLikPoisson, Y=Y, X=X, W=tau[, 2], 
                     control=list(fnscale=-1))$par
    diff <- max(abs(c(alphaNew, betaNew)-c(alpha, beta)))
    alpha <- alphaNew; beta <- betaNew
    logL[iter] <- LogLikZIP(theta=c(alpha, beta), X=X, Y=Y)
    cat(alpha, beta, diff, logL[iter], '\n')
  }
  theta <- list(alpha=alpha, beta=beta)  
  # logLvec <- logL[1:iter]; logLvec; plot(logLvec, type='b')
  return(list(theta=theta, tau=tau, 
              pi=pi, lambda=lambda, fitted.values=pi*lambda,
              logL=logL[iter], logLvec=logL[1:iter]))
}

