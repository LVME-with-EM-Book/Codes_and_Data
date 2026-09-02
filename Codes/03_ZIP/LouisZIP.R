# Check Louis formulas for the ZIP(pi, lambda) distribution

rm(list=ls())
library(pscl)
library(ellipse)
source('../Functions/FunctionsZIPreg.R')

# Parms
n <- 100; pi = .5; lambda <- 10
thetaTrue <- list(pi=pi, lambda=lambda)

# Functions
SimZIP <- function(n, pi, lambda){
  Z <- rbinom(n, 1, pi)
  Y <- Z*rpois(n, lambda)
  return(list(Z=Z, Y=Y))
  }
FitZIP <- function(Y){
  fit <- zeroinfl(as.matrix(Y, n, 1) ~ 1, dist="poisson")$coef
  # return(list(pi=1-plogis(fit$zero), lambda=exp(fit$count)))
  return(c(pi=1-plogis(fit$zero), lambda=exp(fit$count)))
}
Pi2Tau <- function(pi, lambda = 1){
  pi*exp(-lambda) / ((1 - pi) + pi*exp(-lambda))
}

################################################################################
# tau  = f(pi)
plot(curve(Pi2Tau(x)), type='l', col=0, lwd=2, xlim=c(0, 1), ylim=c(0, 1), 
     xlab=expression(pi), ylab=expression(tau), cex=1.5, cex.axis=1.5, cex.lab=1.5)
abline(0, 1, lwd=2, col=8, lty=2)
lambdaList <- c(0.2, 1, 5)
for(lambda in lambdaList){
  curve(Pi2Tau(x, lambda=lambda), col=which(lambdaList==lambda), lwd=2, add=TRUE)
}
# DerZIP <- function(Y, Z, pi, lambda){
#   Zp <- sum(Z); Yp <- sum(Y); P <- sum(Y > 0)
#   gamma <- 1/pi/(1-pi)
#   return(list(grad = c(gamma*Zp - n/(1-pi), -Zp + Yp/lambda), 
#               hess = matrix(c(-Zp/pi^2 + (n-Zp)/(1- pi)^2, 0, 
#                               0, -Yp / lambda^2), 2, 2)))
# }
VarZIP <- function(Y, pi, lambda){
  n <- length(Y); Yp <- sum(Y); P <- sum(Y > 0)
  alpha <- pi*exp(-lambda)/((1-pi) + pi*exp(-lambda))
  espZ <- (Y > 0) + (Y == 0)*alpha
  espZp <- sum(espZ)
  varZ <- espZ * (1 - espZ)
  varZp <- sum(varZ)
  # espZp <- P + (n-P)*alpha
  # varZp <- -(n-P)*alpha*(1-alpha) # SR 29nov2024: pourquoi '-' ?
  gamma <- 1/pi/(1-pi)
  varGradient <- varZp * matrix(c(gamma^2, -gamma, 
                                  -gamma, 1), 2, 2)
  # eigen(varGradient)$values
  espHessian <- matrix(c(- espZp/pi^2 - (n - espZp)/(1-pi)^2, 0, 
                         0, -Yp/lambda^2), 2, 2)
  # eigen(espHessian)$values
  # eigen(espHessian + varGradient)$values
  return(list(espZp=espZp, varZp=varZp, varGradient=varGradient, 
              espHessian=espHessian, asympVar=-solve(espHessian + varGradient)))
}

# Simulations
B <- 1000
theta <- matrix(NA, B, 2)
# grad <- matrix(NA, B, 2); hess <- array(NA, dim=c(B, 2, 2)) 
varGradient <- espHessian <- asymVar <- array(NA, dim=c(B, 2, 2))
stat <- matrix(NA, B, 2)
chi2 <- rep(NA, B)
for(b in 1:B){ # b <- 1
  if(b %% round(sqrt(B))==0){cat(b, '')}
  sim <- SimZIP(n, thetaTrue$pi, thetaTrue$lambda)
  # der <- DerZIP(sim$Y, sim$Z, pi, lambda)
  # grad[b, ] <- der$grad; hess[b, , ] <- der$hess
  theta[b, ] <- unlist(FitZIP(sim$Y))
  varZIP <- VarZIP(sim$Y, theta[b, 1], theta[b, 2])
  varGradient[b, , ] <- varZIP$varGradient
  espHessian[b, , ] <- varZIP$espHessian
  asymVar[b, , ] <- varZIP$asympVar
  stat[b, ] <- (theta[b, ] - unlist(thetaTrue))/sqrt(diag(asymVar[b, , ]))
  chi2[b] <- t(theta[b, ] - unlist(thetaTrue)) %*% solve(asymVar[b, , ]) %*%
    (theta[b, ] - unlist(thetaTrue))
}

par(mfrow=c(2, 2))
hist(stat[, 1], freq=FALSE, breaks=sqrt(B))
curve(dnorm(x), from=min(stat[, 1]), to=max(stat[, 1]), add=TRUE, col=2, lwd=2)

hist(stat[, 2], freq=FALSE, breaks=sqrt(B))
curve(dnorm(x), from=min(stat[, 2]), to=max(stat[, 2]), add=TRUE, col=2, lwd=2)

plot(stat, pch=20)
lines(ellipse(diag(2), scale = c(1, 1), centre = c(0, 0), level = 0.95,
        t = sqrt(qchisq(.95, 2)), which = c(1, 2), npoints = 100), lwd=2, col=2)
abline(h=0, v=0)

statMax <- max(chi2)
hist(chi2, freq=FALSE, breaks=sqrt(B), xlim=c(0, statMax))
curve(dchisq(x, df=2), from=0, to=statMax, add=TRUE, col=2, lwd=2)

solve(cov(theta - rep(1, B)%o%unlist(thetaTrue)))
apply(-varGradient-espHessian, c(2, 3), mean)

cov(theta - rep(1, B)%o%unlist(thetaTrue))
apply(asymVar, c(2, 3), mean)

################################################################################
# Data
dataName <- 'CodBarents'
data <- read.table('../../Data/03_ZIP/BarentsFish.csv', sep=';', head=TRUE)
Covariates <- as.matrix(data[, (2:5)])
Counts <- data[, 6:ncol(data)]
# load('BarentsFish.Rdata')
j <- 21 # Species Ga_mo
Abundance <- Counts[, j]
n <- length(Abundance); p <- 1+ncol(Covariates)
fit <- zeroinfl(Abundance ~  1, dist="poisson")
fit$alpha <- -fit$coefficients$zero; fit$beta <- fit$coefficients$count
# fit <- EMZIPreg(X=matrix(1, n, 1), Y=Abundance, presence=TRUE)
cat(signif(c(fit$alpha, fit$beta), 4))

pi <- plogis(fit$alpha)
lambda <- exp(fit$beta)
cat(signif(c(pi, lambda), 4))

louis <- VarZIP(Y=Abundance, pi=pi, lambda=lambda)
cat(signif(sqrt(diag(louis$asympVar)), 4))

piCI <- pi + sqrt(louis$asympVar[1, 1])*qnorm(p=c(.025, 0.975))
signif(piCI, 4)
lambdaCI <- lambda + sqrt(louis$asympVar[2, 2])*qnorm(p=c(.025, 0.975))
signif(lambdaCI, 4)

