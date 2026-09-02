# PCA for Villeger data analysis: fish species in Terminos Lagoon (Gulf of Mexico)
# GrP13-Book-MultiAnalEcolData, chap 12

rm(list=ls()); par(mfrow=c(4, 3), pch=20); palette('R3')
library(mvtnorm); library(mclust); library(fields)
source('../Functions/FunctionsProbabilisticPCA.R')
figDir <- '../../Figures/'
dataDir <- '../../Data/03_PPCA/'
exportFig <- FALSE
scale <- TRUE

# Data
dataName <- 'Villeger2012c' # Villeger2012c
traits <- read.table(paste0(dataDir, dataName, '_AJ-traits.csv'), sep=';', dec=',', header=TRUE)
rownames(traits) <- traits[, 1]
Y <- as.matrix(traits[, -1])
head(Y); n <- nrow(Y); p <- ncol(Y)

# Description
# if(exportFig){png(paste0(figDir, dataName, '-scatter.png'))}
# plot(as.data.frame(Y), pch=20, cex=.5)
# if(exportFig){dev.off()}

library(GGally)
graph_width <- 12.5
graph_height <- 12.5
graph_unit <- "in"
ggpair_plot <- ggpairs(Y) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

if(exportFig){
  ggsave(paste0(figDir, dataName, '-scatter.png'),
         plot = ggpair_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit) 
}

# Scaling
dataName <- paste0(dataName, '-scale', scale)


if(scale){
  Y <- scale(Y); S <- cor(Y)}else{
  Y <- Y - rep(1, n)%o%colMeans(Y); S <- cov(Y)*(n-1)/n; 
}

eigS <- eigen(S)

# PCA
pca <- list(sdev=sqrt(eigS$values), x=(Y%*%eigS$vectors), rotation=eigS$vectors, scale=scale)
plot(pca$x, col=0)
text(pca$x, col=1, labels=rownames(Y))
abline(v=0, h=0)

# pPCA
ppcaList <- list()
for(q in 1:ncol(Y)){
  # Mettre include_mu_pen = TRUE pour mettre le $p$ de la pénalité de mu
  ppcaList[[q]] <- NaivePPCA(Y, q=q, scale=scale)
}
logL <- unlist(lapply(ppcaList, function(ppca){ppca$logL}))
aic <- unlist(lapply(ppcaList, function(ppca){ppca$aic}))
qAic <- which.max(aic)
bic <- unlist(lapply(ppcaList, function(ppca){ppca$bic}))
qBic <- which.max(bic)

if(exportFig){png(paste0(figDir, dataName, '-selection.png'))}
yLim <- range(c(logL, bic, aic))
plot(logL, type='b', col=1, ylim=yLim, lwd=2,  
     xlab='q', ylab='log-likelihood', cex=1.5, cex.axis=1.5, cex.lab=1.5, pch=20)
points(aic, type='b', col=4, lwd=2, lty=2, pch=21)
abline(v=qAic, col=4, lty=2, lwd=2)
text(qAic, min(yLim)+.9*diff(yLim), label=expression(q[AIC]), col=4, cex=2)
points(bic, type='b', col=2, lwd=2, lty=3, pch=24)
abline(v=qBic, col=2, lty=3, lwd=2)
text(qBic, min(yLim)+.9*diff(yLim), label=expression(q[BIC]), col=2, cex=2)
if(exportFig){dev.off()}

# Estimation of Sigma
if(exportFig){png(paste0(figDir, dataName, '-eigen.png'))}
par(pch=20, mex=1)
plot(pca$sdev^2, type='b', lwd=2, xlab='q', ylab='Eigenvalues', 
     cex=1.5, cex.axis=1.5, cex.lab=1.5, pch=20)
points(eigen(ppcaList[[qAic]]$Sigma)$values, type='b', col=4, lwd=2, pch=21, lty=2)
abline(v=qAic, col=4, lty=2, lwd=2)
text(qAic, .9*max(pca$sdev^2), label=expression(q[AIC]), col=4, cex=2)
points(eigen(ppcaList[[qBic]]$Sigma)$values, type='b', col=2, lwd=2, pch=24, lty=3)
abline(v=qBic, col=2, lty=3, lwd=2)
text(qBic, .9*max(pca$sdev^2), label=expression(q[BIC]), col=2, cex=2)
if(exportFig){dev.off()}

# Cumulated
if(exportFig){png(paste0(figDir, dataName, '-cumEigen.png'))}
par(pch=20, mex=1)
plot(cumsum(pca$sdev^2), type='b', lwd=2, xlab='q', ylab='Cumulated eigenvalues', 
     cex=1.5, cex.axis=1.5, cex.lab=1.5, pch=20)
points(cumsum(eigen(ppcaList[[qAic]]$Sigma)$values), type='b', col=4, lwd=2, pch=21, lty=2)
abline(v=qAic, col=4, lty=2, lwd=2)
text(qAic, .9*sum(pca$sdev^2), label=expression(q[AIC]), col=4, cex=2)
points(cumsum(eigen(ppcaList[[qBic]]$Sigma)$values), type='b', col=2, lwd=2, pch=24, lty=3)
abline(v=qBic, col=2, lty=3, lwd=2)
text(qBic, .9*sum(pca$sdev^2), label=expression(q[BIC]), col=2, cex=2)
if(exportFig){dev.off()}

# Shrinkage
q <- qBic
pca$m <- pca$x[, 1:q] %*% diag(1/pca$sdev[1:q])
colnames(pca$m) <- paste0('dim ', (1:q))
ppca <- NaivePPCA(Y, q=q, scale=scale)
M <- CondMomentsPPCA(Y, ppca)$M
for(j in 1:q){M[, j] <- sign(cor(M[, j], pca$x[, j]))*M[, j]}
gamma <- sqrt(eigS$values[1:q] - ppca$sigma2)
shrinkage <- gamma/sqrt(eigS$values[1:q])
print(ppca$sigma2)
print(round(rbind(sqrt(eigS$values[1:q]), gamma, shrinkage), 3))

PlotShrinkage <- function(pca, M, pcList){
  if(exportFig){png(paste0(figDir, dataName, '-shrinkage-dim', pcList[1], '-dim', pcList[2], '.png'))}
  plot(pca$m[, pcList], pch=20, cex.lab=1.5, cex=1.5)
  points(M[, pcList], col=2, pch=20, cex=1.5)
  for(i in 1:n){lines(rbind(pca$m[i, pcList], M[i, pcList]), col=8, lwd=2)}
  abline(h=0, v=0)
  if(exportFig){dev.off()}
}

PlotShrinkage(pca=pca, M=M, pcList=1:2)
PlotShrinkage(pca=pca, M=M, pcList=3:4)

# Correlation circle
condMoments <- CondMomentsPPCA(Y, ppca)
image.plot(1:(p+qBic), 1:(p+qBic), cor(cbind(Y, condMoments$M)))

CorCircle <- function(Y, M, pcList){
  if(exportFig){png(paste0(figDir, dataName, '-corCircle-dim', pcList[1], '-dim', pcList[2], '.png'))}
  corCircle <- cor(Y, M)
  plot(0, 0, xlim=c(-1, 1), ylim=c(-1, 1), xlab=paste('dim', pcList[1]), ylab=paste('dim', pcList[2]), col=0)
  for(j in 1:p){
    lines(c(0, corCircle[j, pcList[1]]), c(0, corCircle[j, pcList[2]]))
    text(corCircle[j, pcList[1]], corCircle[j, pcList[2]], labels=colnames(Y)[j])
  }
  angle <- seq(0, 2*acos(-1), length.out=200)
  lines(cos(angle), sin(angle), col=8); abline(h=0, v=0)
  if(exportFig){dev.off()}
}

CorCircle(Y=Y, M=condMoments$M, pcList=1:2)
CorCircle(Y=Y, M=condMoments$M, pcList=3:4)

Contribution <- function(ppca, pcList){
  if(exportFig){png(paste0(figDir, dataName, '-corCircle-dim', pcList[1], '-dim', pcList[2], '.png'))}
  plot(0, 0, xlab=paste('dim', pcList[1]), ylab=paste('dim', pcList[2]), col=0, 
       xlim=range(ppca$B[, pcList[1]]), ylim=range(ppca$B[, pcList[2]]))
  for(j in 1:p){
    lines(c(0, ppca$B[j, pcList[1]]), c(0, ppca$B[j, pcList[2]]))
    text(ppca$B[j, pcList[1]], ppca$B[j, pcList[2]], labels=colnames(Y)[j])
  }
  abline(h=0, v=0)
  if(exportFig){dev.off()}
}

Contribution(ppca=ppca, pcList=1:2)
Contribution(ppca=ppca, pcList=3:4)

