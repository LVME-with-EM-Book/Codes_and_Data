# PCA for Villeger data analysis: fish species in Terminos Lagoon (Gulf of Mexico)
# GrP13-Book-MultiAnalEcolData, chap 12

rm(list=ls()); par(pch=20); palette('R3')
seed <- 1; set.seed(seed)
library(mvtnorm);
source('../Functions/FunctionsProbabilisticPCA.R')
source('../Functions/FunctionsProbabilisticPCA.R')
figDir <- '../../Figures/'
dataDir <- '../../Data/03_PPCA/'
exportFig <- FALSE
scale <- FALSE

# Data
dataName <- 'Villeger2012c' # Villeger2012c
traits <- read.table(paste0(dataDir, dataName, '_AJ-traits.csv'), sep=';', dec=',', header=TRUE)
rownames(traits) <- traits[, 1]
Y <- as.matrix(traits[, -1]); n <- nrow(Y); p <- ncol(Y)
if(scale){Y <- scale(Y)}

# Missing values
missRate <- .2
R <- matrix(rbinom(n*p, 1, 1-missRate), n, p)
Yobs <- Y; Yobs[which(R==0)] <- NA

# pPCA
qMax <- p; 
resFile <- paste0(dataDir, dataName, '-ppcaMiss-miss', 100*missRate, '-qMax', qMax, '.Rdata')
if(!file.exists(resFile)){
  par(mfrow=c(ceiling(sqrt(qMax)), round(sqrt(qMax))))
  ppcaList <- list()
  for(q in 1:qMax){
    ppcaList[[q]] <- EMpPCAmiss(Yobs=Yobs, R=R, q=q)
    plot(ppcaList[[q]]$logLpath, main=q, ylim=quantile(ppcaList[[q]]$logLpath, probs=c(0.1, 1), na.rm=TRUE), xlab='')
  }
  # par(mfrow=c(ceiling(sqrt(qMax)), round(sqrt(qMax))))
  # for(q in (qMax-1):1){
  #   init <- ppcaList[[q+1]]; init$B <- init$B[, -(q+1), drop=FALSE]
  #   ppcaTmp <- EMpPCAmiss(Yobs=Yobs, R=R, q=q, init=init)
  #   plot(ppcaTmp$logLpath, main=q, ylim=quantile(ppcaTmp$logLpath, probs=c(0.1, 1), na.rm=TRUE), xlab='')
  #   if(ppcaTmp$logL > ppcaList[[q]]$logL){ppcaList[[q]] <- ppcaTmp}
  # }
  save(ppcaList, file=resFile)
}else{load(resFile)}
par(mfrow=c(ceiling(sqrt(qMax)), round(sqrt(qMax))))
for(q in 1:qMax){
  plot(ppcaList[[q]]$logLpath, main=q, ylim=quantile(ppcaList[[q]]$logLpath, probs=c(0.1, 1), na.rm=TRUE), xlab='')
}

# Eigenvalues
par(mfrow=c(3, 2))
eigSigma <- t(sapply(1:qMax, function(q){eigen(ppcaList[[q]]$Sigma)$values}))
plot(eigSigma[1, ], col=1, type='b', ylim=c(0, max(eigSigma)))
for(h in 1:qMax){points(eigSigma[h, ], col=h, type='b')}
points(2:qMax, eigSigma[1:(qMax-1), p], type='b', lty=2, col=8)
points(1:qMax, eigSigma[1:qMax, 1], type='b', lty=2, col=8)
abline(h=0)

# eigSigma <- t(sapply(1:qMax, function(q){eigen(ppcaList[[q]]$Sigma)$values}))
# plot(eigSigma[, 1], col=1, type='b', ylim=c(0, max(eigSigma)))
# for(h in 1:qMax){points(eigSigma[, h], col=h, type='b')}
# points(2:qMax, eigSigma[1:(qMax-1), p], type='b', lty=2, col=8)
# points(1:qMax, eigSigma[1:qMax, 1], type='b', lty=2, col=8)
# abline(h=0)

# Selection
logL <- unlist(lapply(ppcaList, function(ppca){ppca$logL}))
aic <- unlist(lapply(ppcaList, function(ppca){ppca$aic}))
qAic <- which.max(aic)
bic <- unlist(lapply(ppcaList, function(ppca){ppca$bic}))
qBic <- which.max(bic)
plot(logL, type='b', ylim=range(c(logL, aic, bic)))
lines(aic, type='b', col=4); abline(v=qAic, col=4, lty=2)
lines(bic, type='b', col=2); abline(v=qBic, col=2, lty=2)

# Imputation
q <- qBic; ppca <- ppcaList[[q]]; 
imputeBic <- ImputePPCAall(Yobs=Yobs, R=R, ppca=ppca)
Yfull <- Yobs; Yfull[which(R==0)] <- imputeBic$Yimp[which(R==0)]
plot(Y[which(R==0)], imputeBic$Yimp[which(R==0)]); abline(h=0, v=0, a=0, b=1)

# Imputation with increasing q
imputeList <- lapply(1:qMax, function(q){ImputePPCAall(Yobs=Yobs, R=R, ppca=ppcaList[[q]])})
plot(Y[which(R==0)], imputeList[[1]]$Yimp[which(R==0)]); abline(h=0, v=0, a=0, b=1)
imputeMat <- as.vector(imputeList[[1]]$Yimp[which(R==0)])
for(q in 2:qMax){
  points(Y[which(R==0)], imputeList[[q]]$Yimp[which(R==0)], col=q)
  imputeMat <- cbind(imputeMat, as.vector(imputeList[[q]]$Yimp[which(R==0)]))
}
boxplot(imputeMat - as.vector(Y[which(R==0)])); abline(h=0)
points(apply(imputeMat, 2, sd), type='b', col=2)
points(apply(imputeMat - as.vector(Y[which(R==0)]), 2, quantile, prob=0.75), type='b', col=4)
points(apply(imputeMat - as.vector(Y[which(R==0)]), 2, quantile, prob=0.25), type='b', col=4)
points(-apply(imputeMat, 2, sd), type='b', col=2)

# # Imputation Iterative
# YfullInit <- Yfull
# plot(Y[which(R==0)], imputeBic$Yimp[which(R==0)]); abline(h=0, v=0, a=0, b=1)
# imputeMat <- as.vector(Yfull[which(R==0)])
# hMax <- 10
# for(h in 2:hMax){
#   ppcaFull <- NaivePPCA(Y=Yfull, q=q)
#   imputeEach <- ImputePPCAeach(Yfull=Yfull, R=R, ppca=ppcaFull)
#   points(Y[which(R==0)], imputeEach$Yimp[which(R==0)], col=h)
#   Yfull <- Yobs; Yfull[which(R==0)] <- imputeEach$Yimp[which(R==0)]
#   imputeMat <- cbind(imputeMat, as.vector(Yfull[which(R==0)]))
# }
# boxplot(imputeMat - as.vector(Y[which(R==0)])); abline(h=0)
# 
