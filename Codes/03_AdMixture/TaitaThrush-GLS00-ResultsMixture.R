# Results of genotype mixture for Taita thrush data (GLS00)

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

# Librairies --------------------------------------------------------------

library(tidyverse) # For data manipulation
# library(traitor)

# Graph parameters --------------------------------------------------------

source("../utils_figure_parameters.R")
source("../utils_ggplot_theme.R")

# Data & fit
kMax <- 5
dataName <- 'TaitaThrush-m8-5-6-3-4-10-8'
# algoParms <- '-mixture-k5-seed1-try20'
algoParms <- '-mixture-k5-seed1-try20-tol1e-06'
load(paste0(dataDir, dataName, '.Rdata'))
load(paste0(dataDir, dataName, algoParms, '.Rdata'))
n <- nrow(genotypes); p <- ncol(genotypes);
popNb <- length(unique(pop))

################################################################################
# Mixture
alleles <- sapply(1:p, function(j){sort(unique(haplotypes[, j]))})
m <- sapply(alleles, length)
logL <- sapply(1:kMax, function(k){emList[[k]]$logLik})
ent <- sapply(1:kMax, function(k){-sum(emList[[k]]$tau*log(emList[[k]]$tau))})
penDim <- (0:(kMax-1)) + (1:kMax)*(sum(m-1))
aic <- logL - penDim
bic <- logL - penDim*log(n)/2
bicMax <- max(bic)
bicProb <- exp(bic-bicMax) / sum(exp(bic-bicMax))
icl <- bic - ent
rbind(logL, penDim, ent, bic, icl, bicProb)

mixture_est_multiple_ks <- data_frame(k = 1:5,
                                      ICL = icl,
                                      AIC = aic,
                                      BIC = bic,
                                      "Log-lik" = logL)

criterions_plot <- mixture_est_multiple_ks %>% 
  pivot_longer(-k) %>% 
  ggplot() + 
  aes(x = k, y = value, color = name) + 
  geom_line() + 
  geom_point(aes(shape = name), size = 3) +
  scale_x_continuous(breaks = 1:5) +
  labs(x = "Number of components",
       y = "Criterion value",
       color =  "",
       shape = "") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.15, 0.82)) +
  geom_vline(xintercept = 3, linetype = 2) +
  scale_color_viridis_d(option = "C") +
  scale_shape_manual(values = c(4, 2, 6, 0)) +
  guides(
    color = guide_legend(
      override.aes = list(linetype = c("dashed"))
    )
  )
if(exportFig){
  ggsave("../../Figures/TaitaThrush-mixture-BIC.png",
         plot = criterions_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)
}

# Fit
k <- which.max(bic)
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
Copy

custom_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                    "#ebe221", "#e377c2", "white", "#bcbd22", "#17becf")

plot_list <- map(1:7, function(j) {
  t(em$gamma[, j, ]) %>%
    as.data.frame() %>%
    mutate_all(function(x) replace_na(x, 0)) %>%
    mutate(Marker = j) %>%
    rowid_to_column(var = "Allele") %>%
    pivot_longer(cols = -c("Allele", "Marker"),
                 names_to = "Population", values_to = "Percentage") %>%
    mutate(Allele = factor(Allele, levels = 1:10)) %>%
    ggplot(aes(x = Population, y = Percentage, fill = Allele)) +
    geom_col(color = "lightgray", width = .95) +
    coord_cartesian(expand = FALSE) +
    theme(legend.position = "none") +
    labs(title = paste("Marker", j)) +
    scale_fill_manual(values = custom_palette) +
    theme(plot.title = element_text(hjust = 0.5, family = "LM Roman 10")) +
    labs(x = "", y = "")
})
design <- "
AAA#BBB#CCC#DDD
##EEE#FFF#GGG##
"
my_plot <- patchwork::wrap_plots(plot_list, design = design)
my_plot
if(exportFig)
  ggsave("../../Figures/TaitaThrush-all-markers.png",
       plot = my_plot,
       width = 2*graph_width,
       height = graph_height,
       units = graph_unit)

# Clustering
par(mfrow=c(1, 1), pch=20)
if(exportFig){png(paste0(figDir, dataName, '-mixture-clustering.png'))}

bp <- barplot(t(em$tau[order(pop), ]), col=4+(1:k), border=4+(1:k)) 
bpCoef <- mean(diff(bp))
abline(v=bpCoef*(cumsum(table(pop))[-popNb]), col=1, lwd=3, lty=2)
popPos <- c(0, cumsum(table(pop)))
popPos <- popPos[-(popNb+1)] + diff(popPos)/2
for(h in 1:popNb){text(bpCoef*popPos[h], .5, label=popName[h], cex=1.5)}
if(exportFig){dev.off()}

cluster_individuals <- em$tau[order(pop), ] %>% 
  as.data.frame() %>% 
  mutate(Origin = factor(popName[pop[order(pop)]], 
                         levels = c("Ngangao", "Chawia", 
                                    "Mbololo", "Yale"))) %>% 
  group_by(Origin) %>% 
  mutate(Individual = 1:n()) %>% 
  ungroup() %>% 
  pivot_longer(cols = - c("Individual", "Origin"), 
                 names_to = "Population", 
               values_to = "Probability",
               names_prefix = "V")

plot_cluster_list <- group_by(cluster_individuals, Origin) |> 
  group_map(function(.x, .g){
    origin_ <- as.character(.g$Origin[1])
    xlab = ""
    ylab = ""
    if(origin_ %in% c("Ngangao", "Chawia"))
      ylab = "Probability"
    if(origin_ %in% c("Mbololo"))
      xlab = "Individual"
    ggplot(mutate(.x, Origin = origin_), aes(x = Individual, y = Probability)) + 
      geom_col(aes(fill = Population), width = 1, color = "lightgray", linewidth = 0.1) +
      facet_grid(~Origin) +
      scale_fill_viridis_d(option = "C") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) + 
      theme(legend.position = "none") +
      labs(x = xlab, y = ylab) +
      theme(text = element_text(family = "LM Roman 10"))
  })
design <- "
AAAAAAAAAAAAAAAAAAAA
BBBBBCCCCCCCCCCCCCDD
"
plot_cluster_individuals <- patchwork::wrap_plots(plot_cluster_list, design = design)

if(exportFig){
  ggsave(paste0(figDir, dataName, '-mixture-clustering.png'),
         plot = plot_cluster_individuals,
         width = graph_width,
         height = graph_height / 2,
         units = graph_unit)
}
Zhat <- apply(em$tau, 1, which.max)
table(pop, Zhat)

