require(ggplot2)
# Unified theme for figures
theme_set(theme_bw(base_family = "LM Roman 10"))
# Dimensions of output graphs (as asked by ggsave)
graph_unit <- "px"
graph_width <- 1500
graph_height <- 1500
exportFig <- FALSE
##########################"
#--------------- Packages 
library(ggplot2)
library(dplyr)

########################## Functions
source('SAEM-MCMC.R')
#######################
# ---- Data handling --------------------------------------------- 

# Database
lfmc.data <- read.table("data_lfmc.txt", header = TRUE)
str(lfmc.data)

# Data are ordered for temporal analysis
lfmc.data <- lfmc.data[with(lfmc.data,order(site, plot, time)), ] 

# "Plot" is coded as factor
lfmc.data$plot <- as.factor(lfmc.data$plot)
lfmc.data$leaf.type = with(lfmc.data, factor(leaf.type, levels=c("Grass W", "Grass E", "M. spinosum", "S. bracteolactus")))
levels(lfmc.data$leaf.type)[levels(lfmc.data$leaf.type) == "S. bracteolactus"] <- "S.  filaginoides"
lfmc.data$plot = with(lfmc.data, factor(plot, levels=c("4", "5", "6", "1", "2", "3")))
lfmc.data$plot = with(lfmc.data, factor(plot))

str(lfmc.data)

g <- ggplot(data = lfmc.data, aes(x = time, y = lfmc)) + 
  geom_point(aes(color=plot)) +
  xlab("Time (days)") +  facet_wrap(~leaf.type,scales = "free") + 
  ylab("LFMC (%)") # The LFMC temporal dynamics is plotted by "leaf type" 
g <- g +  labs(colour = "Plot") +theme(axis.text.x = element_text(size = 14),  # graduations axe X
                                       axis.text.y = element_text(size = 14),
                                       axis.title.x = element_text(size = 14),
                                       axis.title.y = element_text(size = 14),
                                       strip.text = element_text(size = 14))
g
if(exportFig)
  ggsave(paste0("LFMC_data_all.png"),
         plot = g ,
         height = graph_height,
         width = graph_width,
         units = graph_unit)


###################################### 
### Logistic

A = 10; 
m = 30
w = 85
s = 10 
phi = list(A= A,w=w,m=m,s=s)
phi2 <- phi; phi2$s = 5

t <- seq(0,80, length.out = 1000)
df2 <- tibble(Time = rep(t,2),
              LFMC = c(FourParamLogis(t,phi),FourParamLogis(t,phi2)))
df2$s <-as.factor(c(rep(phi$s,1000),rep(phi2$s,1000)))
gg <- ggplot(df2, aes(Time, LFMC,group=s,colour = s)) +   geom_line(size = 1.0) 
gg <- gg +  geom_vline(xintercept =phi$m, size = 0.5, linetype = "dashed",col = "#C4961A") 
gg <- gg +  geom_hline(yintercept =c(phi$A,phi$w), size = 0.5, linetype = "dashed",col = "#C4961A")
gg <- gg + annotate("text",y=-10,x=phi$m,label="m",col="#C4961A")+coord_cartesian(ylim=c(0,100),xlim=c(0,80),clip="off")
gg <- gg + annotate("text",y=phi$A,x=-5,label="A",col="#C4961A")
gg <- gg + annotate("text",y=phi$w,x=-5,label="w",col="#C4961A") 
gg <- gg +theme(axis.text.x = element_text(size = 14),  # graduations axe X
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                strip.text = element_text(size = 14))

if(exportFig)
  ggsave(paste0("LFMC_FourParamsLogis.png"),
         plot = gg ,
         height = graph_height,
         width = graph_width,
         units = graph_unit)

#################################"" 

lfmc.data$group <- with(lfmc.data, factor(plot:leaf.type)) # one col ("group") crossing both factors ("plot" x "leaf type") is generated
levels(lfmc.data$group)
lfmc.data$group=with(lfmc.data, factor(group, labels=c("GW plot 4","GW plot 5", "GW plot 6",
                                                       "GE plot 1","SM plot 1", "SS plot 1",
                                                       "GE plot 2","SM plot 2", "SS plot 2",
                                                       "GE plot 3","SM plot 3", "SS plot 3"))) # the factor "group" is relabeled according to the manuscript

lfmc.data.gd <- groupedData(lfmc ~ time | group, data = lfmc.data, order.groups=FALSE) # and a groupedData class object is created 
str(lfmc.data.gd)
lfmc.data.gd$plot = with(lfmc.data.gd, factor(plot, levels=c("4", "5", "6", "1", "2", "3"))) # "plot" is coded as factor in the groupedData object

# ---- Data exploration ------------------------------------------- 

par(mfrow=c(1,1), cex=0.85)
plot(lfmc.data$time, lfmc.data$lfmc, xlab="Time (days)", ylab="LFMC (%)", ylim=c(0,350), main="") # LFMC decreases over time
par(mfrow=c(1,2), cex=0.85)
plot(lfmc.data$leaf.type, lfmc.data$lfmc, xlab="Leaf type", ylab="LFMC (%)", ylim=c(0,350), col="grey", main="") # Shrubs are more heterogeneous
hist(lfmc.data$lfmc, main="", xlab="LFMC (%)", col="grey")


