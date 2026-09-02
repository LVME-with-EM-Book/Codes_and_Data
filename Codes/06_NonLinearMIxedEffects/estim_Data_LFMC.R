rm(list=ls())

require(ggplot2)
# Unified theme for figures
theme_set(theme_bw(base_family = "LM Roman 10"))
# Dimensions of output graphs (as asked by ggsave)
graph_unit <- "px"
graph_width <- 1500
graph_height <- 1500
exportFig <- FALSE

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scam)
source('SAEM-MCMC.R')
################################################
#-------------------les data --------------------
################################################"

#------------------------------- 
lfmc.data <- read.table("data_lfmc.txt", header = TRUE)
lfmc.data$plot <- as.factor(lfmc.data$plot)
lfmc.data$leaf.type = with(lfmc.data, factor(leaf.type, levels=c("Grass W", "Grass E", "M. spinosum", "S. bracteolactus")))
levels(lfmc.data$leaf.type)[levels(lfmc.data$leaf.type) == "S. bracteolactus"] <- "S.  filaginoides"
lfmc.data$plot_fact = with(lfmc.data, factor(plot), levels=c("4", "5", "6", "1", "2", "3"))
lfmc.data$plot <- as.numeric(as.character(lfmc.data$plot_fact))

str(lfmc.data)


myData <- lfmc.data[lfmc.data$leaf.type=='Grass E',]
nPlot <- 3


#------------ Plot des datas
cols <- setNames(
  scales::hue_pal()(6)[4:6],
  levels(myData$plot.fact)
)
g <- ggplot(data = myData, aes(x = time, y = lfmc)) + 
  geom_point(aes(color=plot_fact)) + geom_smooth(colour="gray", linewidth=0.7)+
  xlab("Time (days)") + #+  facet_wrap(~leaf.type,scales = "free") + 
  ylab("LFMC (%)") # The LFMC temporal dynamics is plotted by "leaf type" 
g +scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols)+theme(axis.text.x = element_text(size = 14),  # graduations axe X
                                         axis.text.y = element_text(size = 14))


################################################ 
#########   initialisation 
############################################ 

Z.init <- initZ(myData)
theta.init <-initTheta(Z.init,myData)






###################################################""
#------------------- ESTIM SAEM 
######################################################## 
paramsAlgo = list(nbIterMCMC=100,nbIterSAEM = 100,rho=c(1/100,1/10,1,2),keep = FALSE,print=FALSE)
paramsAlgo$gammaSAEM = c(rep(1,20),1/(1:paramsAlgo$nbIterSAEM)^(0.99))


#res_SAEM_MCMC <- SAEM_MCMC_NLME(myData,paramsAlgo,omegaDiag = TRUE)

#df_res_SAEM_MCMC <- as.data.frame(cbind(res_SAEM_MCMC$mu,res_SAEM_MCMC$sigma2))
#names(df_res_SAEM_MCMC) <- c('log(A)', 'log(w)', 'log(m)', 'log(s)',"sigma2")
#df_res_SAEM_MCMC$iterations  = 1:101
#save(df_res_SAEM_MCMC,file='res_SAEM_LFMC.Rdata')

load('res_SAEM_LFMC.Rdata')
df_res_SAEM_MCMC <- pivot_longer(
  df_res_SAEM_MCMC,
  cols = c('log(A)', 'log(w)', 'log(m)', 'log(s)',"sigma2"),
  names_to = "parameter",
  values_to = "value"
)

lab <- c(
  "log(A)" = expression(mu[log(A)]),
  "log(w)" = expression(mu[log(w)]),
  "log(m)" = expression(mu[log(m)]),
  "log(s)" = expression(mu[log(s)]),
  "sigma2" = expression("sigma^2")
)

ggiter <- ggplot(df_res_SAEM_MCMC, aes(x = iterations, y = value)) +
  geom_line() +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2,labeller = labeller(
    parameter = c(
      "log(A)" = "mu[log(A)]",
      "log(w)" = "mu[log(w)]",
      "log(m)" = "mu[log(m)]",
      "log(s)" = "mu[log(s)]",
      "sigma2" = "sigma^2"
    ),
    .default = label_parsed
  ))
ggiter <- ggiter + geom_vline(xintercept = 20,
                              colour = "grey80",
                              linewidth = 0.5)
ggiter <- ggiter +  theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          theme(axis.text.x = element_text(size = 14),  # graduations axe X
                                axis.text.y = element_text(size = 14),
                                axis.title.x = element_text(size = 14),
                                axis.title.y = element_text(size = 14),
                                strip.text = element_text(size = 14))
)
ggiter

if(exportFig)
  ggsave(paste0("LFMC_param_iterations.png"),
         plot = ggiter ,
         height = graph_height,
         width = graph_width,
         units = graph_unit)



################# posterior 
theta.estim = list(mu = res_SAEM_MCMC$mu[paramsAlgo$nbIterSAEM+1,])
theta.estim$sigma2 = res_SAEM_MCMC$sigma2[paramsAlgo$nbIterSAEM+1]
theta.estim$Omega = res_SAEM_MCMC$Omega[paramsAlgo$nbIterSAEM+1,,]

paramsAlgo2 <- paramsAlgo
paramsAlgo2$keep = TRUE
paramsAlgo2$nbIterMCMC<- 5000
paramsAlgo2$print=TRUE
#post_Ind_param_MCMC <- MCMC_NLME(res_SAEM_MCMC$Z,myData,theta.estim,paramsAlgo2)

#save(post_Ind_param_MCMC,file='res_postSample_LFMC.Rdata')
load('res_postSample_LFMC.Rdata')
seq_iter = seq(paramsAlgo2$nbIterMCMC/2, paramsAlgo2$nbIterMCMC,by=2)


gg <- ggplot()
for (i in 1:nPlot){  
  Data.i <- myData%>%filter(plot == i)
  seq_time <- seq(0,max(Data.i$time)+1,len=100)
  Z.post.i <- post_Ind_param_MCMC[seq_iter,i,]
  F.i <- sapply(1:nrow(Z.post.i),function(m){
    FourParamLogis(seq_time,exp(Z.post.i[m,]))}) 
  d1 = dim(F.i)[1]
  d2 = dim(F.i)[2]
  F.i = F.i #+ sqrt(theta.estim$sigma2)*matrix(rnorm(d1*d2,0,1),d1,d2)
  q.i <-  apply(F.i, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  Data_simpost.i=data.frame(time=seq_time )
  Data_simpost.i$plot_fact  = with(Data_simpost.i, factor(i, levels=c("4", "5", "6", "1", "2", "3")))
  
  
  Data_simpost.i$q05<- q.i[1,]
  Data_simpost.i$q95<- q.i[2,]
  
  gg <- gg +
    geom_ribbon(data = Data_simpost.i,
                aes(x = time, ymin = q05, ymax = q95,
                    fill = plot_fact, colour = plot_fact),
                alpha = 0.20,
                linewidth = 0.4
    ) + geom_point(
      data = Data.i,
      aes(x = time, y = lfmc, colour = plot_fact),
      size = 2
    ) 
}
gg <- gg +scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) + xlab("Time (days)") + 
  ylab("LFMC (%)")  + labs(
    colour = "Plot",
    fill = "Plot"
  ) + facet_wrap(~leaf.type)+theme(axis.text.x = element_text(size = 14),  # graduations axe X
                                   axis.text.y = element_text(size = 14),
                                   axis.title.x = element_text(size = 14),
                                   axis.title.y = element_text(size = 14),
                                   strip.text = element_text(size = 14))

gg
if(exportFig)
  ggsave(paste0("LFMC_post_ind_traj.png"),
         plot = gg ,
         height = graph_height,
         width = graph_width,
         units = graph_unit)


