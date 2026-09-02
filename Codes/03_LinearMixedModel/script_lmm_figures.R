library(tidyverse)

exportFig <- FALSE

donnees_NH4 <- read.table("../../Data/03_LinearMixedModel/donnees_NH4.txt",
                          sep = ";", header = TRUE)

graph_width <- 5
graph_height <- 5
graph_unit <- "in"

quadrat_levels <- donnees_NH4 %>% 
  group_by(Sol, Quadrat) %>% 
  summarise(mN = mean(logNH4)) %>% 
  ungroup(Quadrat) %>% 
  mutate(m = mean(mN)) %>% 
  ungroup() %>% 
  arrange(m, mN) %>% 
  pull(Quadrat)

nh4_plot <- ggplot(donnees_NH4 %>% 
                     mutate(Quadrat = factor(Quadrat, 
                                             levels = quadrat_levels,
                                             labels = paste0("Plot-", 1:9))),
                   aes(y =logNH4, x = Quadrat, fill = Sol)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.75, 0.25)) + 
  labs(x = "", fill = "Soil type",
       y = expression(NH[4]~"("~log[10]~")")) +
  scale_fill_viridis_d() 

if(exportFig)
  ggsave("../../Figures/borneo-NH4-plot.png",
         plot = nh4_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)

