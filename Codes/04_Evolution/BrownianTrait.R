# Figure for Brownian evolution of a trait

rm(list=ls())
# seed <- 4, 12, 19
seed <- 4; set.seed(seed)
figDir <- '../../Figures/'
figName <- 'BrownianTrait'
exportFig <- FALSE

# Tree + branch lengths
mu <- 0; sigma <- 1
n <- 5
d <- c(1, 1, 2.25, 1.5, 1.5, 1.25, 1.5, .75)
pa <- c(6, 6, 8, 7, 7, 8, 9, 9)
off <- lapply(n+(1:(n-1)),function(j){which(pa==j)})
h <- rep(0, 2*n-1)
for(j in n+(1:(n-1))){h[j] <- h[off[[j-n]][1]] + d[off[[j-n]][1]]}
hTot <- max(h)

# Function
BrownPath <- function(mu, sigma=1, d, step=.001){
  t <- seq(0, d, by=step)
  w <- mu +  sqrt(step)*sigma*cumsum(rnorm(length(t)))
  return(list(t=t, w=w))
}

# Branching path
z <- rep(NA, 2*n-1); z[2*n-1] <- mu
pathList <- list(); yLim <- mu+sigma*c(-1, 1)
for(j in (n-1):1){# j <- n-1
  pathList[[j]] <- list()
  for(i in off[[j]]){
    pathList[[j]][[i]] <- BrownPath(mu=z[n+j], sigma=sigma, d=d[i])
    yLim <- range(c(yLim, pathList[[j]][[i]]$w))
    # lines(hTot-h[n+j]+path$t, path$w, col=j)
    z[i] <- pathList[[j]][[i]]$w[length(pathList[[j]][[i]]$w)]
    # points(hTot-h[n+j]+path$t[length(path$t)], z[i], pch=20, cex=1.5)
    # text(hTot-h[n+j]+path$t[length(path$t)], z[i], label=paste0('Z', i), col=2)
  }
}

# Plot (horizontal)
if(exportFig){png(paste0(figDir, figName, '-h-seed', seed, '.png'), height=360)}
plot(0, 0, col=0, xlim=c(0, 1.1*hTot), ylim=yLim, xlab='', ylab='', 
     axes=0, cex.lab=1.5)
for(j in (n-1):1){# j <- n-1
  for(i in off[[j]]){
    lines(hTot-h[n+j]+pathList[[j]][[i]]$t, pathList[[j]][[i]]$w, col=8)
    points(hTot-h[i], z[i], pch=20, cex=1.5)
    # text(hTot-h[n+j]+path$t[length(path$t)], z[i], label=paste0('Z', i), col=2)
  }
}
points(0, z[2*n-1], pch=20, cex=1.5)
lines(c(0, 0), yLim, lwd=2)
# text(hTot/40, mean(yLim), label='trait', cex=1.5, srt=90)
text(0, max(yLim), label='trait', cex=1.5, pos=4)
lines(c(0, hTot), min(yLim)*c(1, 1), lwd=2)
text(hTot/2, min(yLim)+diff(yLim)/30, label='time', cex=1.5)
text(hTot-h[1], z[1], label=expression(Y[1]), cex=1.5, pos=4)
text(hTot-h[2], z[2], label=expression(Y[2]), cex=1.5, pos=4)
text(hTot-h[3], z[3], label=expression(Y[3]), cex=1.5, pos=4)
text(hTot-h[4], z[4], label=expression(Y[4]), cex=1.5, pos=4)
text(hTot-h[5], z[5], label=expression(Y[5]), cex=1.5, pos=4)
text(hTot-h[6], z[6], label=expression(Z[6]), cex=1.5, pos=4)
text(hTot-h[7], z[7], label=expression(Z[7]), cex=1.5, pos=4)
text(hTot-h[8], z[8], label=expression(Z[8]), cex=1.5, pos=4)
text(hTot-h[9], z[9], label=expression(Z[9]), cex=1.5, pos=4)
if(exportFig){dev.off()}

# Plot (vertical)
if(exportFig){png(paste0(figDir, figName, '-v-seed', seed, '.png'), height=360)}
plot(0, 0, col=0, xlim=yLim, ylim=c(0, 1.1*hTot), xlab='', ylab='', 
     axes=0, cex.lab=1.5)
for(j in (n-1):1){# j <- n-1
  for(i in off[[j]]){
    lines(pathList[[j]][[i]]$w, h[n+j]-pathList[[j]][[i]]$t, col=8)
    points(z[i], h[i], pch=20, cex=1.5)
    # text(hTot-h[n+j]+path$t[length(path$t)], z[i], label=paste0('Z', i), col=2)
  }
}
points(z[2*n-1], 0, pch=20, cex=1.5)
lines(yLim, c(0, 0), lwd=2)
# text(hTot/40, mean(yLim), label='trait', cex=1.5, srt=90)
text(.9*max(yLim), hTot/40, label='trait', cex=1.5, pos=4)
lines(min(yLim)*c(1, 1), c(0, hTot), lwd=2)
text(min(yLim)+diff(yLim)/40, hTot/2, label='time', cex=1.5, srt=90)
text(z[1], h[1], label=expression(Y[1]), cex=1.5, pos=3)
text(z[2], h[2], label=expression(Y[2]), cex=1.5, pos=3)
text(z[3], h[3], label=expression(Y[3]), cex=1.5, pos=3)
text(z[4], h[4], label=expression(Y[4]), cex=1.5, pos=3)
text(z[5], h[5], label=expression(Y[5]), cex=1.5, pos=3)
text(z[6], h[6], label=expression(Z[6]), cex=1.5, pos=3)
text(z[7], h[7], label=expression(Z[7]), cex=1.5, pos=3)
text(z[8], h[8], label=expression(Z[8]), cex=1.5, pos=3)
text(z[9], h[9], label=expression(Z[9]), cex=1.5, pos=3)
if(exportFig){dev.off()}
