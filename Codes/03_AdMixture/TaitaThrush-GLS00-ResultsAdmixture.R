# Results of genotype admixture for Taita thrush data (GLS00)

rm(list=ls()); palette('R3'); par(mfrow=c(1, 1), pch=20)

# -> Book: TaitaThrush-k5-m8-5-6-3-4-10-8-seed1-try20.Rdata

# Dirs
library(extraDistr); library(tensor); library(MSCquartets); library(sirt); 
library(ade4); library(mclust)
source('../Functions/FunctionsGenotypeMixture.R')
source('../Functions/FunctionsLatex.R')
dataDir <- '../../Data/TaitaThrush/'
dataName <- 'TaitaThrush'
figDir <- '../../Figures/'
exportFig <- FALSE

# Data & fit
# dataName <- 'TaitaThrush-m22-5-6-3-4-10-8'
dataName <- 'TaitaThrush-m8-5-6-3-4-10-8'
resParms <- '-admixture-k5-seed1-try0-tol1e-06'
load(paste0(dataDir, dataName, '.Rdata'))
load(paste0(dataDir, dataName, resParms, '.Rdata'))
kMax <- length(emList)
n <- nrow(genotypes); p <- ncol(genotypes);
popNb <- length(unique(pop))

# ################################################################################
# initParms <- '-admixtureFake-k5-seed1-try0-tol1e-06'
# load(paste0(dataDir, dataName, initParms, '.Rdata'))
# initList <- emList; emList <- c()
# kMax <- length(initList)
# # Fast fit from admixtureFake
# resFile <- paste0(dataDir, dataName, algoParms, '.Rdata')
# if(!file.exists(resFile)){
#   emList <- list()
#   for(k in 1:kMax){
#     emList[[k]] <- VEMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
#                                           k=k, tauInit=initList[[k]]$tau, tol=tol)
#   }
#   save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
# }

################################################################################
# Admixture
alleles <- sapply(1:p, function(j){sort(unique(haplotypes[, j]))})
m <- sapply(alleles, length)
elbo <- sapply(1:kMax, function(k){emList[[k]]$elbo})
ent <- sapply(1:kMax, function(k){-sum(emList[[k]]$tau*log(emList[[k]]$tau))})
penDim <- 1:kMax + (1:kMax)*(sum(m-1))
bic <- elbo - penDim*log(n)/2
bicMax <- max(bic)
bicProb <- exp(bic-bicMax) / sum(exp(bic-bicMax))
icl <- bic - ent
rbind(elbo, penDim, ent, bic, icl, bicProb)
if(exportFig){png(paste0(figDir, dataName, '-mixture-BIC.png'))}
plot(elbo, type='b', ylim=range(c(elbo, bic, icl)), pch=0, 
     ylab='log-likelihood', xlab=expression(K), cex=1.5, cex.axis=1.5, cex.lab=1.5)
lines(bic, type='b', col=2, pch=1, cex=1.5); abline(v=which.max(bic), col=2)
lines(icl, type='b', col=4, pch=2, cex=1.5); abline(v=which.max(icl), col=4)
if(exportFig){dev.off()}

# Fit
k <- which.max(bic)
k <- 3
em <- emList[[k]]
par(mfrow=c(ceiling(sqrt(p)), round(sqrt(p))), mex=.6)
dimnames(em$gamma)[[1]] <- 1:k
for(j in 1:p){
  if(exportFig){png(paste0(figDir, dataName, '-mixture-marker', j, '-gamma.png'))}
  barplot(t(em$gamma[, j, 1:m[j]]), main=colnames(genotypes)[j], 
          col=1+(1:m[j]),
          # col=gray.colors(m[j]), 
          cex=1.5, cex.axis=1.5)
  if(exportFig){dev.off()}
}

# Clustering
par(mfrow=c(1, 1), pch=20)
# bp <- barplot(t(em$pi[order(pop), ]), border=1+(1:k), col=0)
if(exportFig){png(paste0(figDir, dataName, '-mixture-clustering.png'))}
bp <- barplot(t(em$alphaTilde[order(pop), ]), col=4+(1:k), border=4+(1:k)) 
bpCoef <- mean(diff(bp))
# popPch <- c(0, 1, 2, 6); 
# points(bp, rep(.5, n), col=sort(pop), pch=popPch[sort(pop)], cex=1.5)
abline(v=bpCoef*(cumsum(table(pop))[-popNb]), col=1, lwd=3, lty=2)
popPos <- c(0, cumsum(table(pop)))
popPos <- popPos[-(popNb+1)] + diff(popPos)/2
for(h in 1:popNb){text(bpCoef*popPos[h], .5, label=popName[h], cex=1.5)}
if(exportFig){dev.off()}

if(k==3){
  simplexPrepare()
  for(i in 1:n){simplexPoint(em$alphaTilde[i, ]/sum(em$alphaTilde[i, ]), col=pop[i], pch=14+pop[i], cex=1.5)}
  simplexPoint(em$alpha/sum(em$alpha), cex=3, pch='+', col=8)
}
Zhat <- apply(em$alphaTilde, 1, which.max)
table(pop, Zhat)

