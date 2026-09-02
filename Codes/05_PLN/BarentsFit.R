# Illustration of the PLN model on Barents data

rm(list=ls())
library(PLNmodels)
library(combinat)

# Dir
dataName <- 'BarentsFish'
dataDir <- '../../Data/Barents/'

# Data
data <- read.csv(paste0(dataDir, dataName, '.csv'), sep=';')
X <- as.matrix(data[, 2:5]); Y <- as.matrix(data[, 6:ncol(data)])
Xscale <- scale(X)
save(X, Y, Xscale, file=paste0(dataDir, dataName, '.Rdata'))
# load(paste0(dataDir, dataName, '.Rdata'))

# PLN
plnNull <- PLN(Y ~ 1, control=PLN_param(config_post=list(jackknife=TRUE, variational_var=TRUE)))
plnFull <- PLN(Y ~ Xscale, control=PLN_param(config_post=list(jackknife=TRUE, variational_var=TRUE)))
dim(vcov(plnFull))
save(plnNull, plnFull, file=paste0(dataDir, dataName, '-All.Rdata'))
# load(paste0(dataDir, dataName, '-All.Rdata'))

# Model list
d <- ncol(Xscale); n <- nrow(Xscale); intercept <- matrix(1, n, 1); colnames(intercept) <- 'Intercept'
modelList <- list(); modelNum <- 1
modelList[[modelNum]] <- list(model=c(0), X=intercept)
for(c in 1:d){
  combin <- t(combn(1:d, c))
  for(m in 1:nrow(combin)){
    modelNum <- modelNum+1;
    Xtmp <- cbind(intercept, Xscale[, combin[m, ]])
    colnames(Xtmp) <- c('Intercept', colnames(X)[combin[m, ]])
    modelList[[modelNum]] <- list(model=c(0, combin[m, ]), X=Xtmp)}
}
save(modelList, file=paste0(dataDir, dataName, '-modelList.Rdata'))
# load(paste0(dataDir, dataName, '-modelList.Rdata'))

# Fit list
plnList <- list()
for(m in 1:length(modelList)){plnList[[m]] <- PLN(Y ~ -1 + modelList[[m]]$X)}
save(plnList, file=paste0(dataDir, dataName, '-plnList.Rdata'))
# load(paste0(dataDir, dataName, '-plnList.Rdata'))

