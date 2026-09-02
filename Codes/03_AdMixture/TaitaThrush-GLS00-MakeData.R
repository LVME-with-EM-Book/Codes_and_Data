# Fit genotype mixture & admixuture for Taita thrush data (GLS00)

rm(list=ls()); palette('R3'); par(pch=20)
seed <- 1; set.seed(seed)

# -> Book: TaitaThrush-k5-m8-5-6-3-4-10-8-seed1-try20.Rdata

# Dirs
# library(extraDistr); library(tensor); library(MSCquartets); library(sirt); 
# library(ade4); library(mclust)
source('../Functions/FunctionsGenotypeMixture.R')
source('../Functions/FunctionsLatex.R')
dataDir <- '../../Data/03_Admixture/'
dataName <- 'TaitaThrush'
figDir <- '../../Figures/'
exportFig <- FALSE

# Data 
haplotypes <- as.matrix(read.table(paste0(dataDir, 'thrush-data.str')))
haplotypes[which(haplotypes==-9)] <- NA
covar <- haplotypes[, 1:2]; haplotypes <- haplotypes[, -(1:2)]; 
id <- covar[, 1]; pop <- covar[, 2]
n2 <- nrow(haplotypes); n <- n2/2; p <- ncol(haplotypes)
colnames(haplotypes) <- paste0('M', 1:p)
popName <- c('Chawia', 'Ngangao', 'Mbololo', 'Yale')

# Description 
sapply(1:p, function(j){table(haplotypes[, j])})
m <- sapply(1:p, function(j){length(unique(haplotypes[, j], na.rm=TRUE))})
haploArray <- MakeGenoArray(haplotypes)

# Reduce nb of alleles
mMax <- 10; redList <- which(m > mMax)
for(j in redList){
  sel <- which(!is.na(haplotypes[, j]))
  gmm <- Mclust(haplotypes[sel, j], modelNames='E', G=mMax)
  haplotypes[sel, j] <- gmm$classification
}
sapply(1:p, function(j){table(haplotypes[, j])})
m <- sapply(1:p, function(j){length(unique(haplotypes[, j], na.rm=TRUE))})
haploArray <- MakeGenoArray(haplotypes)
F_Matrix2Tabular(cbind(rep(1:3, each=2), haplotypes[1:6, ], popName[pop[1:6]]))

# Genotypes for init
genotypes <- matrix(0, n, p)
for(i in 1:n){for(j in 1:p){
  tmp <- sort(c(haplotypes[2*(i-1)+1, j], haplotypes[2*(i-1)+2, j]))
  genotypes[i, j] <- paste0(tmp[1], '-', tmp[2])
}}
pop <- covar[2*(1:nrow(genotypes)), 2]
for(j in 1:ncol(genotypes)){print(j); print(unique(genotypes[, j]))}

dataName <- paste0(dataName, '-m')
for(j in 1:(p-1)){dataName <- paste0(dataName, m[j], '-')}
dataName <- paste0(dataName, m[p])
dataFile <- paste0(dataDir, dataName, '.RData')
save(covar, pop, id, haplotypes, haploArray, genotypes, popName, file=dataFile)
