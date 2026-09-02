# Results of genotype fake admixture (fixed proba prob_i) for Taita thrush data (GLS00)

rm(list=ls()); palette('R3'); par(pch=20)

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
kMax <- 5
dataName <- 'TaitaThrush'
fitName <- '-admixture-k5-m7-5-6-3-3-10-8-seed1-try20'
load(paste0(dataDir, dataName, fitName, '.Rdata'))
n <- nrow(genotypes); p <- ncol(genotypes);

################################################################################
# Fake admixture (with fixed prob_i)

# Dirichlet for the pi_i
dirList <- list()
dirList[[1]] <- list(alpha=1, alpha0=1, xsi=1, logLik=0)
for(k in 2:kMax){
  dirList[[k]] <- dirichlet.mle(emList[[k]]$pi)
  dirList[[k]]$logLik <- sum(sapply(1:n, function(i){ddirichlet(emList[[k]]$pi[i, ], dirList[[k]]$alpha, log=TRUE)}))
  }

# Results
par(mfrow=c(1, 1))
alleles <- sapply(1:p, function(j){sort(unique(haplotypes[, j]))})
m <- sapply(alleles, length)
logL <- sapply(1:kMax, function(k){emList[[k]]$logLik})
penDim <- (0:(kMax-1)) + (1:kMax)*(sum(m-1))
# logL <- sapply(1:kMax, function(k){emList[[k]]$logLik + dirList[[k]]$logLik})
# penDim <- (1:kMax) + (1:kMax)*(sum(m-1))
ent <- sapply(1:kMax, function(k){-sum(emList[[k]]$tau*log(emList[[k]]$tau))})
bic <- logL - penDim*log(n)/2
bicMax <- max(bic)
bicProb <- exp(bic-bicMax) / sum(exp(bic-bicMax))
icl <- bic - ent
rbind(logL, penDim, ent, bic, icl, bicProb)
plot(logL, type='b', ylim=range(c(logL, bic, icl)))
lines(bic, type='b', col=2); abline(v=which.max(bic), col=2)
lines(icl, type='b', col=4); abline(v=which.max(icl), col=4)

# Fit
k <- which.max(bic)
k <- 3
em <- emList[[k]]
par(mfrow=c(ceiling(sqrt(p)), round(sqrt(p))), mex=.6)
for(j in 1:p){
  barplot(t(em$gamma[, j, 1:m[j]]), main=colnames(genotypes)[j], col=1+(1:m[j]))
}

# Clustering
par(mfrow=c(1, 1), pch=20)
# bp <- barplot(t(em$pi[order(pop), ]), border=1+(1:k), col=0)
bp <- barplot(t(em$pi[order(pop), ]), col=4+(1:k), border=4+(1:k))
points(bp, rep(-.005, n), col=sort(pop), pch=14+sort(pop))
abline(v=.5*mean(diff(bp))+bp[cumsum(table(pop))][-4], col=1, lwd=3, lty=2)
simplexPrepare()
for(i in 1:n){simplexPoint(em$pi[i, ], col=pop[i], pch=14+pop[i], cex=1.5)}
simplexPoint(dirList[[k]]$xsi, cex=3, pch='+', col=8)
Zhat <- apply(em$pi, 1, which.max)
table(pop, Zhat)

