rm(list=ls())
library(ggplot2)
library(tidyverse)
library(parallel)
library(viridis)
library(dplyr)
library(xtable)
library(cowplot)
#----------------------------------------- 
source('Functions_EM.R')
source("../utils_ggplot_theme.R")
source("../utils_figure_parameters.R")
#---------------------------------------- 

fig_path <- "../../Figures/EMToyExample/"
exportFig <- FALSE

#=======================================
#-------------- Data ---------------- 
#==========================================
library(palmerpenguins)

bill_length_distribution <- ggplot(penguins %>% na.omit(),
                                   aes(x=bill_length_mm)) +
  geom_histogram(color = "black",fill="cyan4", breaks = seq(32, 60, by = 1)) +
  labs(x = "Bill length (mm)",
       y = "Counts") +
  scale_x_continuous(breaks = seq(35, 55, by = 10)) +
  # Needed to fit side by side with penguin photo
  theme(plot.margin = margin(t = 1, #Top margin
                             r = 0.1,  # Right margin
                             l = 0, 
                             b = 0,
                             unit = "cm"))
if(exportFig)
  ggsave(filename = paste0(fig_path, "pinguinsBillLength.png"),
         plot = bill_length_distribution,
         width = graph_width, height = graph_height, units = graph_unit)

##########################################################
#------------  EM for all parameters
##########################################################

y <- penguins$bill_length_mm
y <- y[!is.na(y)]
#---------------- RUNs  EM from various starting points

KM <- kmeans(y, 2)

params0 = list(list(mu=c(40,50),var=c(5,5),p=0.5),
               list(mu=c(20,50),var=c(5,5),p=0.5),
               list(mu=c(35,70),var=c(5,5),p=0.6),
               list(mu=c(50,40),var=c(10,10),p=0.4))
params0[[5]] = list(mu=c(40,50),var=c(1,1),p=0.5)
params0[[6]] = list(mu=KM$centers,var=c(3,3),p=0.5)
res_all_run_all <-do.call(rbind,lapply(1:length(params0),function(run){
  res_EM_run  <- EM_2Mixture(params0[[run]],y,tol=10^{-7},estim_only_mu = FALSE)
  niter <- nrow(res_EM_run)
  res_EM_run <- cbind(res_EM_run,rep(run,niter))
  res_EM_run <- cbind(res_EM_run, c(1:nrow(res_EM_run)))
  names(res_EM_run)[7]='Init'
  names(res_EM_run)[8]='numIter'
  
  return(res_EM_run)}))
#------------------------------------------------
tabl <- matrix(unlist(params0),ncol=5,byrow=TRUE)
tabl <- cbind(tabl,rep(0,length(params0)))
row.names(tabl) = 1:length(params0)
for (init in 1:length(params0)){
  w <- which(res_all_run_all$Init==init)
  tabl[init,6]<- max(res_all_run_all$loglik[w])
}
xtable(tabl)


res_all_run_all$Init <- as.factor(res_all_run_all$Init)
res_all_run_all <- res_all_run_all %>% filter(numIter>3)
ll_trajectory <- ggplot(res_all_run_all,aes(x=numIter - 4,
                                            y=loglik,colour = Init)) + 
  geom_line() + 
  geom_point(size = .5) +
  # theme(axis.text=element_text(size=15),
  #       axis.title=element_text(size=15),
  #       legend.text = element_text(size=15)) +
  labs(y = "Log likelihood",
       x = "Iterations")
if(exportFig)
  ggsave(filename = paste0(fig_path, "em_Iter_allparams.png"),
         plot = ll_trajectory,
         width = graph_width, height = graph_height, units = graph_unit)


#scale_fill_hue()$palette(5)  
##########################################################
#------------  EM on only MU 
##########################################################
theta_true = list(mu = c(35, 45))
theta_true$var = c(9,9)
theta_true$p = 0.36
Zsim = sample(c(1,2),length(y),prob=c(theta_true$p,1-theta_true$p),replace=TRUE)
ysim = theta_true$mu[Zsim] +sqrt(theta_true$var[Zsim])*rnorm(length(y))
hist(ysim)
N = 600

res_like = data.frame(mu1=double(),
                      mu2 = double(), 
                      ll = double())
grid_mu1 = seq(20,60,len=N)
grid_mu2 = seq(20,60,len=N)
grid_params <- as.data.frame(matrix(0,N*N,2))
names(grid_params) <- c('mu1', 'mu2') 
grid_params[,1]<-  rep(grid_mu1,N)
grid_params[,2] <- rep(grid_mu2,each=N)
res_loglik <- mclapply(1:N^2,function(i){
  params.i = theta_true
  params.i$mu <- c(grid_params[i,1],grid_params[i,2])
  logLikelihood(params.i,y=ysim)},mc.cores = 4)
res_like = grid_params
res_like[,3]  = unlist(res_loglik)
names(res_like)[3] <- 'loglik'

wmax = which.max(res_like[,3])
grid_params[wmax,]

#res_all_ru
#---------------- RUNs  EM from various starting points on only MU 

theta0 = theta_true
params0 <- lapply(1:5,function(i){theta0})
params0[[1]]$mu = c(25,25)
params0[[2]]$mu = c(30,25)
params0[[3]]$mu = c(55,40)
params0[[4]]$mu = c(25,55)
params0[[5]]$mu = c(25,45)




res_all_run <-do.call(rbind,lapply(1:length(params0),function(run){
  res_EM_run  <- EM_2Mixture(params0[[run]],ysim,tol=10^{-5},estim_only_mu = TRUE)
  niter <- nrow(res_EM_run)
  res_EM_run <- cbind(res_EM_run,rep(run,niter))
  res_EM_run <- cbind(res_EM_run, 1:niter)
  names(res_EM_run)[7]='Init'
  names(res_EM_run)[8]='numIter'
  
  return(res_EM_run)}))



#------------ PLOT ------------------

initial_points <- purrr::map(params0, "mu") %>% 
  map_dfr(function(x){
    data.frame(mu1 = x[1],
               mu2 = x[2])
  }) %>% 
  mutate(Init = as.factor(1:5))

final_points <- res_all_run %>%
  group_by(Init) %>%
  filter(numIter == max(numIter)) %>%
  filter(Init %in% c(1, 2, 4))

# Ajouter les log-vraisemblances finales au dataframe des points finaux
final_points <- final_points %>%
  mutate(loglik_final = format(round(loglik, 2), nsmall = 2)) # Formater les valeurs


res_all_run$Init <- as.factor(res_all_run$Init)
res_all_run$numIter <- as.factor(res_all_run$numIter)
true_mu <- as.data.frame(matrix(theta_true$mu,nrow=1))
names(true_mu)  =c('mu1','mu2')
true_mu_df <- as.matrix(true_mu) %>% as.data.frame()
init_labels <-c(
  "1" = "Init 1 (25, 25)",
  "2" = "Init 2 (30, 25)",
  "3" = "Init 3 (55, 40)",
  "4" = "Init 4 (25, 55)",
  "5" = "Init 5 (25, 45)"
)

best_point <- arrange(final_points %>% ungroup(), loglik_final) %>% 
  slice(1)
em_traj_2d <- ggplot(res_like,aes(x=mu1,y=mu2)) + 
  geom_tile(aes(fill=loglik)) + 
  scale_fill_viridis() + 
  stat_contour(aes(z = loglik), color="white", size=0.25, breaks = c(-5000, -4000, -3000, -2000, -1500, -1300, -1200, -1100)) + 
  geom_point(data=res_all_run %>% dplyr::filter(as.numeric(numIter) > 1), 
             aes(color=Init),size=2) + 
  geom_line(data=res_all_run, aes(color =Init),linewidth=1) + 
  geom_point(shape = 4, data = initial_points, size = 3) + 
  scale_shape_manual(values=c(15, 20)) + 
  # geom_point(data=true_mu, aes(x=mu1, y=mu2,z=NULL),size=2,color='black') + 
  # theme(axis.text=element_text(size=30),
  #       axis.title=element_text(size=30),
  #       legend.text = element_text(size=30),
  #       legend.title = element_text(size=25)) +
  guides(shape="none") + 
  labs(x =expression(mu[1]),
       y = expression(mu[2]),
       fill = "Log-likelihood") + 
  geom_text(data = final_points, 
            aes(x = mu1, y = mu2, label = loglik_final), 
            color = "black", size = 3, vjust = -1) +
  geom_point(data = best_point, color = "darkgreen", size = 3) +
  geom_point(data = true_mu_df, color = "red", shape = 3, size = 3) +
  coord_equal(expand = FALSE) + 
  scale_colour_discrete(
    labels = init_labels, # Libellés personnalisés
    name = "Initialisation" # Titre de la légende
  ) + 
  theme(legend.position = "none", 
        axis.title = element_text(family = "LM Roman 10"))

if(exportFig)
  ggsave(filename = paste0(fig_path, "EM_convergence.pdf"),device = cairo_pdf,
         width = graph_width, height = graph_height, units = graph_unit)
# Transformée en pdf puis convertie en png ensuite






