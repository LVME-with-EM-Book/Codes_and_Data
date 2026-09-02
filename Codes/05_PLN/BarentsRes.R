# Illustration of the PLN model on Barents data

rm(list=ls()); par(lwd=2, pch=20, mfrow=c(1, 1))
library(PLNmodels); library(fields); library(RColorBrewer)
library(tidyverse)
# See http://www.sthda.com/english/wiki/colors-in-r for RColorBrewer
exportFig <- FALSE
source("../utils_figure_parameters.R")
source("../utils_ggplot_theme.R")

#-------------------------------------------------------------------------------
# Dir
dataName <- 'BarentsFish'
dataDir <- '../../Data/Barents/'
figDir <- '../../Figures/'
alpha <- .05
palette <- "PiYG"

# Data & fit
load(paste0(dataDir, dataName, '.Rdata'))
load(paste0(dataDir, dataName, '-All.Rdata'))
load(paste0(dataDir, dataName, '-modelList.Rdata'))
load(paste0(dataDir, dataName, '-plnList.Rdata'))
n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

#-------------------------------------------------------------------------------
# Full model: covariates'effect
#-------------------------------------------------------------------------------
beta <- matrix(as.numeric(plnFull$model_par$B), d+1, p)
varBeta <- matrix(diag(vcov(plnFull)), d+1, p); 
rownames(beta) <- rownames(varBeta) <- colnames(modelList[[length(modelList)]]$X)
colnames(beta) <- colnames(varBeta) <- colnames(Y)
# image.plot(1:d, 1:p, beta[-1, ], xlab='covariates', ylab='species', cex.lab=1.5)

betaBrk <- seq(-8, 8, by=2)
# speciesOrder <- order(prcomp(t(beta[-1, ]))$x[, 1])
# heatmap(t(beta[-1, speciesOrder]), Rowv=NA, Colv=NA, 
#         breaks=betaBrk, col=brewer.pal(n=length(betaBrk)-1, name = palette))

if(exportFig){png(paste0(figDir, dataName, '-plnFull-beta.png'))}
image.plot(1:d, 1:p, beta[-1, ], xlab='', ylab='', axes=0, 
           breaks=betaBrk, col=brewer.pal(n=length(betaBrk)-1, name = palette))
axis(1, labels=colnames(Xscale), at=1:d, cex.axis=1.25)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

est_coeff_plot <- beta[-1,] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Species") %>% 
  pivot_longer(cols = -c("Species"), 
               names_to = "Covariate", values_to = "Estimate") %>% 
  ggplot(aes(x = Covariate,
             y = Species, 
             fill = Estimate)) +
  geom_raster() +
  scale_fill_gradient2(low = "#C51B7D", mid = "white", high = "#4D9221") +
  coord_cartesian(expand = FALSE) +
  labs(x = "", y = "") +
  theme(axis.ticks = element_blank(),
        text = element_text(family = "LM Roman 10", size = 13))
if(exportFig){
  ggsave(filename = paste0(figDir, dataName, '-plnAll-beta.png'),
         unit = graph_unit, height = graph_height, 
         width = graph_width, plot = est_coeff_plot)
}
# Pseudo-stats
stat <- (beta / sqrt(varBeta))[-1, ] # Remove the intercept for the tests
maxBrk <- 6; statBrk <- seq(-maxBrk, maxBrk, by=2)
statPlot <- stat; 
statPlot[which(abs(stat) > maxBrk)] <- maxBrk*sign(stat[which(abs(stat) > maxBrk)])
statBrkLegends <- statBrk; statBrkLegends[1] <- '-Inf'; statBrkLegends[length(statBrk)] <- '+Inf'

if(exportFig){png(paste0(figDir, dataName, '-plnFull-stat.png'))}
image.plot(1:d, 1:p, statPlot, xlab='', ylab='', axes=0, 
           breaks=statBrk, lab.breaks=statBrkLegends, col=brewer.pal(n=length(statBrk)-1, name = palette))
axis(1, labels=colnames(Xscale), at=1:d, cex.axis=1.25)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

# Pseudo-p values
pVal <- 2*pnorm(abs(stat), lower.tail=FALSE)
pAdj <- matrix(p.adjust(pVal, method="BH"), d, p)
signif <- sign(stat) * (pAdj < alpha)

if(exportFig){png(paste0(figDir, dataName, '-plnFull-signif.png'))}
image.plot(1:d, 1:p, statPlot, xlab='', ylab='', axes=0, 
           breaks=statBrk, lab.breaks=statBrkLegends, col=brewer.pal(n=length(statBrk)-1, name = palette))
axis(1, labels=colnames(Xscale), at=1:d, cex.axis=1.25)
axis(2, labels=colnames(Y), at=1:p, las=2)
for(c in 1:d){for(j in 1:p){
  if(signif[c, j]==1){points(c+.1, j, pch=24, cex=1.5, bg=3)}
  if(signif[c, j]==-1){points(c-.1, j, pch=25, cex=1.5, bg=2)}
}}
if(exportFig){dev.off()}

#-------------------------------------------------------------------------------
# Model selection
#-------------------------------------------------------------------------------
dList <- unlist(lapply(modelList, function(model){ncol(model$X)}))
dimList <- p*(p+1)/2 + p*dList
penList <- round(dimList * log(n)/2)
elboList <- round(unlist(lapply(plnList, function(fit){fit$criteria$loglik})))
bicList <- round(unlist(lapply(plnList, function(fit){fit$criteria$BIC})))
iclList <- round(unlist(lapply(plnList, function(fit){fit$criteria$ICL})))
entropyList <- round(unlist(lapply(plnList, function(fit){fit$entropy})))
bicOrder <- order(bicList, decreasing=TRUE)
bicBest <- bicOrder[1]
plnBest <- plnList[[bicBest]]

df_data <- cbind(d = dList - 1, elbo = elboList, bic = bicList) %>% 
  as.data.frame() %>% 
  group_by(d) %>% 
  arrange(elbo) %>% 
  ungroup() %>% 
  mutate(model = 1:nrow(.))

rg_bic <- range(df_data$bic)
rg_elbo <- range(df_data$elbo)

criterions_plot <- ggplot(df_data %>% 
         mutate(trans_bic = bic - rg_bic[1] + rg_elbo[2]), 
       aes(x = model, y = elbo)) +
  geom_point(aes(color = "ELBO", shape = factor(d)),
             size = 3) + 
  geom_line(aes(color = "ELBO")) +
  geom_point(aes(y = trans_bic, color = "Var. BIC", shape = factor(d)),
             size = 3) + 
  geom_line(aes(y = trans_bic, color = "Var. BIC")) +
  theme(text = element_text(family = "LM Roman 10"),
        legend.position = "inside",
        legend.position.inside = c(0.65, 0.18),
        legend.background = element_blank(),
        legend.box = "horizontal",          # Arrange legends in a horizontal box
        legend.direction = "vertical",      # Arrange items within each legend vertically
        legend.spacing.x = unit(1, "cm"), 
        axis.text.y = element_text(color = "#ed7953"),
        axis.text.y.right = element_text(color = "#9c179e")) +
  scale_x_continuous(breaks = c(1, 2, 6, 12, 16)) +
  labs(x = "Model number", color = "",
       shape = "Nbr. covariates") +
  scale_y_continuous(
    name = "ELBO",
    sec.axis = sec_axis(transform = function(x) x + rg_bic[1] - rg_elbo[2], 
                        name="Variational BIC")
  ) +
  scale_color_manual(values = c("#ed7953", "#9c179e")) +
  scale_shape_manual(values = c(1, 3, 2, 8, 19)) +
  geom_vline(xintercept = 11, linetype = 3) +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(shape = NA)
    ),
    shape = guide_legend(ncol = 1)
  )
if(exportFig){
  ggsave(filename = paste0(figDir, dataName, '-pln-bic.png'),
         unit = graph_unit, height = graph_height, 
         width = graph_width, plot = criterions_plot)
}

if(exportFig){png(paste0(figDir, dataName, '-pln-bic.png'))}
par(new=FALSE)
plot(elboList, type='b', ylim=range(elboList, bicList+elboList[1]-bicList[1]), 
     pch=20, lwd=2, xlab='models', ylab='elbo', cex.lab=1.5)
for(c in 1:(d+1)){
  text(mean(which(dList==c)), min(elboList)+.6*diff(range(elboList)), paste0('d=', c), cex=1.25)
}
abline(v=.5+which(diff(dList)>0), col=1, lwd=1)
par(new=TRUE)
plot(bicList, type='b', xlab='', ylab='', col=2, pch=20, lwd=2, axes=0, 
     ylim=range(elboList)+(max(bicList)-max(elboList)))
axis(side=4, col=2, col.ticks=2, col.axis=2)
par(new=FALSE)
abline(v=which.max(bicList), lty=2, col=2, lwd=2)
if(exportFig){dev.off()}

# Three best models
for(m in bicOrder[1:3]){
  cat(m, '& & ', colnames(Xscale)[modelList[[m]]$model[-1]], 
      '&', elboList[m], '&', dimList[m], '&', penList[m], '&', bicList[m], '\\\\ \n')
}
# Null model
m <- 1
cat(m, '& (null) &', colnames(Xscale)[modelList[[m]]$model[-1]], 
    '&', elboList[m], '&', dimList[m], '&', penList[m], '&', bicList[m], '\\\\ \n')
# Full model
m <- 16
cat(m, '& (full) &', colnames(Xscale)[modelList[[m]]$model[-1]], 
    '&', elboList[m], '&', dimList[m], '&', penList[m], '&', bicList[m], '\\\\ \n')

#-------------------------------------------------------------------------------
# Covariance structure
#-------------------------------------------------------------------------------
summary(diag(plnNull$model_par$Sigma))
summary(diag(plnFull$model_par$Sigma))
summary(diag(plnBest$model_par$Sigma))

cov2cor(plnBest$model_par$Sigma) %>% 
  corrplot::corrplot(order = "hclust", hclust.method = "ward.D")
hclust(cov2cor(plnBest$model_par$Sigma), method = "ward.D2")
library(corrplot)
orig_names <- colnames(plnBest$model_par$Sigma)
## Si on veut une visu plus compacte
# levels_plot <- corrMatOrder(cov2cor(plnBest$model_par$Sigma), 
#                             order = "hclust", 
#                             hclust.method = "ward.D2") %>% 
#   orig_names[.]
levels_plot <- sort(orig_names)
  

all_est_cors <- imap_dfr(list(Null = plnNull,
              Full = plnFull,
              Best = plnBest),
         function(x, nm){
           x$model_par$Sigma %>% 
             cov2cor() %>% 
             as.data.frame() %>% 
             rownames_to_column(var = "row") %>% 
             pivot_longer(cols = -c("row"),
                          names_to = "col",
                          values_to = "cor") %>% 
             mutate(model = nm)
         }) %>% 
  mutate(model = factor(model, levels = c("Null", "Best", "Full")),
         col = factor(col, levels = rev(levels_plot)),
         row = factor(row, levels = levels_plot))
all_cors_plot <- ggplot(all_est_cors) +
  aes(x = row, y = col, fill = cor) +
  geom_raster() +
  facet_wrap(~model) +
  scale_fill_gradient2(high = "#800000",
                       mid = "white", 
                       low = "#00008f",
                       limits = c(-1, 1)) +
  labs(x = "", y = "", fill = "Corr.") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        text = element_text(size = 10),
        axis.ticks = element_blank()) +
  coord_equal()
if(exportFig){
  ggsave(filename = paste0(figDir, dataName, '-plnAll-cors.png'),
         unit = graph_unit, height = graph_height, 
         width = 2 * graph_width, plot = all_cors_plot)
}


corBrk <- seq(-1, 1, length.out=11)

if(exportFig){png(paste0(figDir, dataName, '-plnNull-cor.png'))}
image.plot(1:p ,1:p, cov2cor(plnNull$model_par$Sigma), xlab='', ylab='', 
           breaks=corBrk, col=tim.colors(length(corBrk)-1), axes=0)
axis(1, labels=colnames(Y), at=1:p, las=2)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, dataName, '-plnFull-cor.png'))}
image.plot(1:p ,1:p, cov2cor(plnFull$model_par$Sigma), xlab='', ylab='', 
           breaks=corBrk, col=tim.colors(length(corBrk)-1), axes=0)
axis(1, labels=colnames(Y), at=1:p, las=2)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

# # Prediction
# XbetaAll <- cbind(rep(1, nrow(Y)), Xscale)%*%plnFull$model_par$B
# image.plot(cov2cor(plnNull$model_par$Sigma))
# image.plot(cor(XbetaAll))
# plot(cor(XbetaAll), cov2cor(plnNull$model_par$Sigma))

if(exportFig){png(paste0(figDir, dataName, '-plnBest-cor.png'))}
image.plot(1:p ,1:p, cov2cor(plnBest$model_par$Sigma), xlab='', ylab='', 
           breaks=corBrk, col=tim.colors(length(corBrk)-1), axes=0)
axis(1, labels=colnames(Y), at=1:p, las=2)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

plot(plnFull$model_par$Sigma, plnBest$model_par$Sigma, ylim=(range(plnNull$model_par$Sigma)))
points(plnFull$model_par$Sigma, plnNull$model_par$Sigma, col=2)
abline(h=0, v=0, a=0, b=1)

