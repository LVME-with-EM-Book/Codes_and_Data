# Barents fish (Cod) for the ZIP section

rm(list=ls()); palette('R3'); par(pch=20, mfrow=c(1, 1))

# Dirs
dirFig <- '../../Figures/'
exportFig <- FALSE
source('../Functions/FunctionsLatex.R')
source('../Functions/FunctionsZIPreg.R')

################################################################################
# Data
dataName <- 'CodBarents'
data <- read.table('../../Data/03_ZIP/BarentsFish.csv', sep=';', head=TRUE)
Covariates <- as.matrix(data[, (2:5)])
Counts <- data[, 6:ncol(data)]
# load('BarentsFish.Rdata')
j <- 21 # Species Ga_mo
Abundance <- Counts[, j]; Presence <- (Abundance>0)
n <- length(Abundance); p <- 1+ncol(Covariates)
Covariates <- as.matrix(scale(Covariates))
Barents <- as.data.frame(cbind(Covariates, Abundance))
names(Barents)[5] <- 'Abundance'
Covariates <- cbind(rep(1, n), Covariates)
colnames(Covariates)[1] <- 'Intercept'

################################################################################
# Description
F_Matrix2Tabular(head(Barents))
if(exportFig){png(paste0(dirFig, dataName, '-Hist.png'))}
hist(Barents$Abundance, breaks=2*sqrt(nrow(Barents)), main='', xlab='Abundance', cex=1.5, cex.axis=1.5, cex.lab=1.5)
if(exportFig){dev.off()}
sum(Abundance==0)

################################################################################
# Poisson regression model for the Abundance
regP <- glm(Abundance ~ Latitude + Longitude + Depth + Temperature, data=Barents, family='poisson')
summary(regP)
regP$beta <- regP$coefficients
# # Check
# opt <- optim(par=lm(log(1+Abundance) ~ -1+Covariates)$coef, fn=LogLikPoisson, X=Covariates, Y=Abundance, W=rep(1, n), control=list(fnscale=-1))
# c(logLik((regP)), opt$value)
# rbind(regP$coef, opt$par)

################################################################################
# Logistic regression model for the Abundance
regL <- glm(Presence ~  Latitude + Longitude + Depth + Temperature, data=Barents, family='binomial')
summary(regL)
# plot(regL$fitted.values, Presence, pch=20); 
regL$alpha <- regL$coefficients
regL$linPred <- as.vector(Covariates%*%regL$alpha)
logLik(regL)
# # Check
# opt <- optim(par=rep(0, ncol(Covariates)), fn=LogLikBernoulli, X=Covariates, Y=Presence, W=rep(1, n), control=list(fnscale=-1))
# c(logLik((regL)), opt$value)
# rbind(regL$coef, opt$par)

################################################################################
# Zero inflated Poisson regression model for the Abundance
thetaInit <- InitZIP(X=Covariates, Y=Abundance)
# Direct optimization
opt <- optim(unlist(thetaInit), f=LogLikZIP, Y=Abundance, X=Covariates, control=list(fnscale=-1))
opt$alpha <- opt$par[1:ncol(Covariates)]; opt$beta <- opt$par[-(1:ncol(Covariates))]
opt$theta <- list(alpha=opt$alpha, beta=opt$beta)
# Package
library(pscl)
fit <- zeroinfl(Abundance ~ -1 + Covariates, dist="poisson")
fit$theta <- list(alpha=fit$coefficients$zero, beta=fit$coefficients$count)
# EM
em <- EMZIPreg(X=Covariates, Y=Abundance, thetaInit=thetaInit)
c(LogLikZIP(unlist(thetaInit), X=Covariates, Y=Abundance), 
  LogLikZIP(unlist(opt$theta), X=Covariates, Y=Abundance), 
  LogLikZIP(unlist(fit$theta), X=Covariates, Y=Abundance),
  LogLikZIP(unlist(em$theta), X=Covariates, Y=Abundance))
rbind(unlist(thetaInit), unlist(fit$theta), opt$par, unlist(em$theta))
regZIP <- em
regZIP$linPredPresence <- Covariates%*%regZIP$theta$alpha
regZIP$linPredAbundance <- Covariates%*%regZIP$theta$beta
regZIP$tau <- regZIP$tau[, 2]
# regZIP$fitted.values <- regZIP$pi*regZIP$lambda

################################################################################
# Parameter estimates
# regPcoef <- as.vector(regP$coefficients)
# regLcoef <- as.vector(regL$coefficients)
# regZIPcoef <- regZIP$coefficients
# regZIPcoef$presence <- regZIPcoef$alpha
regCoef <- round(cbind(c(regL$alpha, rep(NA, length(regP$beta))), 
                       c(rep(NA, length(regL$alpha)), regP$beta), 
                       as.vector(c(regZIP$theta$alpha, regZIP$theta$beta))), 3)
colnames(regCoef) <- c('Logistic', 'Poisson', 'ZIP')
regCoef <- rbind(c(names(regP$coefficients), names(regL$coefficients)),
                 t(regCoef))
F_Matrix2Tabular(regCoef, labelRowNames='covariate')

################################################################################
# Probability of presence
xGrid <- seq(-4, 4, by=.001)
library(extraDistr)
# Logistic regression
orderRegL <- order(regL$linPred)
if(exportFig){png(paste0(dirFig, dataName, '-PresProb-Logistic.png'))}
plot(regL$linPred, Presence, xlab='Logistic regression', ylab='presence probability',
     pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5, col=2-Presence); 
lines(xGrid, plogis(xGrid), col=4, lwd=2); 
# points(regL$linPred[orderRegL], regZIP$pi[orderRegL], col=4); 
# points(regL$linPred[orderRegL], regZIP$tau[orderRegL], col=3); 
if(exportFig){dev.off()}
# ZIP
orderZIP <- order(regZIP$linPredPresence)
orderZIP0 <- order(regZIP$linPredPresence[which(!Presence)])
if(exportFig){png(paste0(dirFig, dataName, '-PresProb-ZIP_Logistic.png'))}
plot(regZIP$linPredPresence, Presence, xlab='ZIP regression', ylab='presence probability',
     pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5, col=2-Presence); 
lines(xGrid, plogis(xGrid), lwd=2); 
points(regZIP$linPredPresence[orderZIP], regL$fitted.values[orderZIP], pch='+', col=4, cex=1.5); 
if(exportFig){dev.off()}
# Comparison 
if(exportFig){png(paste0(dirFig, dataName, '-PresPred-ZIP-logistic.png'))}
plot(regZIP$tau, regL$fitted.values, 
     xlab='ZIP regression',  ylab='Logistic regression', 
     pch=20, col=1+(Abundance==0), cex=1.5, cex.axis=1.5, cex.lab=1.5)
abline(v=.5, h=.5, col=4, lty=2, lwd=2)
if(exportFig){dev.off()}

# plot(regZIP$pi, regZIP$tau, 
#      xlab=expression(pi),  
#      ylab=expression(tau), 
#      pch=20, col=1+(Abundance==0), cex=1.5, cex.axis=1.5, cex.lab=1.5)
# abline(0, 1, col=8)
# abline(v=.5, h=.5, col=4, lty=2, lwd=2)
# table(regZIP$tau>.5, Presence)
# table(regZIP$pi>.5, Presence)
# table(regL$fitted.values>.5, Presence)

################################################################################
# Expected abundances
# Poisson regression
if(exportFig){png(paste0(dirFig, dataName, '-Abundance-Poisson.png'))}
plot(1+regP$fitted.values, 1+Abundance, xlab='1 + Poisson regression', 
     pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5, log='xy'); 
abline(0, 1, lwd=2)
lines(1+sort(regP$fitted.values), 1+qpois(.025, sort(regP$fitted.values)), col=2, lty=2, lwd=2)
lines(1+sort(regP$fitted.values), 1+qpois(.975, sort(regP$fitted.values)), col=2, lty=2, lwd=2)
if(exportFig){dev.off()}
# ZIP
orderZIP <- order(regZIP$fitted.values)
if(exportFig){png(paste0(dirFig, dataName, '-Abundance-ZIP.png'))}
plot(1+regZIP$fitted.values, 1+Abundance, xlab='1 + ZIP regression', 
     pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5, log='xy'); 
abline(0, 1, lwd=2)
lines(1+regZIP$fitted.values[orderZIP], 1+qzip(.025, regZIP$lambda, 1-regZIP$pi)[orderZIP], col=2, lty=2, lwd=2)
lines(1+regZIP$fitted.values[orderZIP], 1+qzip(.925, regZIP$lambda, 1-regZIP$pi)[orderZIP], col=2, lty=2, lwd=2)
if(exportFig){dev.off()}

if(exportFig){png(paste0(dirFig, dataName, '-Presence-Abundance-ZIP.png'))}
plot(regZIP$pi, 1+regZIP$lambda, xlab=expression(pi[i]), ylab=expression(1+lambda[i]),
     pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5, log='y', col=1+(Abundance==0)); 
if(exportFig){dev.off()}

# ZIP: Poisson part
orderZIP <- order(regZIP$linPredAbundance)
orderZIP <- order(regZIP$fitted.values)

par(mfrow=c(3, 1))
plot(regZIP$fitted.values, regZIP$pi, log='x')
plot(regZIP$fitted.values, regZIP$lambda, log='xy'); abline(0, 1)
plot(1+regZIP$fitted.values, 1+regZIP$fitted.values, log='xy', col=0); abline(0, 1)
lines(1+regZIP$fitted.values[order(regZIP$fitted.values)], 1+qzip(.025, regZIP$lambda[order(regZIP$fitted.values)], 1-regZIP$pi[order(regZIP$fitted.values)]), col=4)
lines(1+regZIP$fitted.values[order(regZIP$fitted.values)], 1+qzip(.975, regZIP$lambda[order(regZIP$fitted.values)], 1-regZIP$pi[order(regZIP$fitted.values)]), col=4)
lines(1+regZIP$fitted.values[order(regZIP$fitted.values)], 1+qzip(.025, regZIP$lambda, 1-regZIP$pi)[order(regZIP$fitted.values)], col=2)
lines(1+regZIP$fitted.values[order(regZIP$fitted.values)], 1+qzip(.975, regZIP$lambda, 1-regZIP$pi)[order(regZIP$fitted.values)], col=2)
points(1+regZIP$fitted.values, 1+Abundance, col=1); abline(0, 1)

################################################################################
# # Comparison with Poisson
# if(exportFig){png(paste0(dirFig, dataName, '-ExpectedAbundance.png'))}
# plot(Abundance, regP$fitted.values, pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5)
# points(Abundance, regZIP$fitted.values, pch=20, col=2)
# abline(0, 1, lwd=2)
# if(exportFig){dev.off()}
# # Classification
# if(exportFig){png(paste0(dirFig, dataName, '-PresenceProbZIP.png'))}
# boxplot(presProbZIP ~ (Abundance > 0), ylab='Presence probability', 
#         pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5); 
# if(exportFig){dev.off()}

################################################################################
## Model comparison
logLik(regL)
logL <- c(logLik(regP), regZIP$logL)
df <- c(p, 2*p)
aic <- logL - df
bic <- logL - 0.5*df*log(n)
compModel <- cbind(df, round(logL, 1), round(bic, 1))
rownames(compModel) <- c('Poisson', 'ZIP')
colnames(compModel) <- c('nb parms', 'logL', 'BIC')
F_Matrix2Tabular(compModel, labelRowNames='model')

LRT <- 2*(regZIP$logL - logLik(regP))
pval <- pchisq(LRT, df=p, lower.tail=FALSE)
c(LRT, pval)
