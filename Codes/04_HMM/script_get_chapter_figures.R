rm(list = ls())

library(tidyverse)
library(sf)
library(geodata)
source("../utils_figure_parameters.R")
source("../utils_ggplot_theme.R")
exportFig <- FALSE

# zone_map --------------------------------------------------------------

zone_map <- bind_rows(gadm("Canada", path = "./", level = 2) %>% 
                        st_as_sf() %>% 
                        filter(str_detect(NAME_1, "Newfoundland"), 
                               as.numeric(CC_2) %in% (1:9)) %>% # Only newfoundland without labrador
                        dplyr::select(COUNTRY), 
                      gadm("Saint Pierre and Miquelon", 
                           path = "./", level = 0) %>% 
                        st_as_sf() %>% 
                        dplyr::select(COUNTRY))

if(!file.exists("../../Data/04_HMM/gannet_trajectory.csv")){
  source("script_create_clean_dataset.R")
}
donnees_brutes <- read.table("../../Data/04_HMM/gannet_trajectory.csv", 
                             header = TRUE, sep = ",",
                             colClasses = c(t = "POSIXct")) 

library(xtable)  
donnees_brutes %>% 
  mutate(t = as.character(t)) %>% 
  head() %>% 
  rename("log10(Step length)" = log10_step) %>% 
  xtable(digits = 3) %>% 
  print(include.rownames = FALSE)



donnees <-   st_as_sf(donnees_brutes,
                      coords = c("Longitude", "Latitude"), crs = st_crs(zone_map)) 

trajectory_plot <- ggplot(zone_map) + 
  geom_sf(fill = "darkgreen") +
  theme_bw() +
  geom_sf(data = donnees %>% 
            st_combine() %>% 
            st_cast("LINESTRING"), alpha = 0.5) +
  geom_sf(data = donnees, mapping = aes(color = t), size = 0.5) + 
  scale_color_gradient(low = "lightgray", high = "black") + 
  theme(legend.position = "none") +
  coord_sf(xlim = c(-57, -52),
           ylim = c(46, 48))

if(exportFig)
  ggsave("../../Figures/gannet-trajectory-plot.png",
         plot = trajectory_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)

step_length_plot <- ggplot(donnees) + 
  aes(x = log10_step) +
  geom_histogram(fill = "lightblue", color = "black", 
                 mapping = aes(y = after_stat(density)), bins = 100) +
  theme_bw() +
  labs(y = "Empirical density", x = expression(log[10]~"of step length"))

if(exportFig)
  ggsave("../../Figures/gannet-steplength-hist.png",
         plot = step_length_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)


all_post_densities <- read.table("results_gannets_post_densities.txt",
                                 sep = ";", header = TRUE)

selected_K <- c(8, 10)
post_densities_plot <- ggplot(all_post_densities %>% 
                                filter(n_states %in% selected_K) %>% 
                                mutate(n_states = factor(paste0("K=", n_states),
                                                         levels = paste0("K=", selected_K)))) + 
  aes(x = log10_step) + 
  geom_histogram(data = map_dfr(selected_K,
                                function(k){
                                  mutate(donnees, 
                                         n_states = factor(paste0("K=", k),
                                                           levels = paste0("K=", 
                                                                           selected_K)))
                                }),
                 mapping = aes(y = after_stat(density)), bins = 100,
                 color = "gray", fill = "lightgray") +
  geom_line(linewidth = 1,
            aes(y = density,
                color = factor(state, 
                               levels = 1:max(selected_K)))) +
  facet_wrap(~n_states, nrow = 2) +
  scale_color_viridis_d(option = "C") +
  labs(color = "State", y = "Empirical density", 
       x = expression(log[10]~"of step length")) +
  theme(legend.position = "inside",
        text = element_text(family = "LM Roman 10"),
        legend.position.inside = c(0.8, 0.28),
        legend.background = element_blank()) +
  guides(color = guide_legend(ncol = 2))

if(exportFig)
  ggsave("../../Figures/gannet-post-state-dist.png",
         plot = post_densities_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)

final_speed_plot <- all_post_densities %>% 
  filter(n_states == 6) %>% 
  ggplot() +
  aes(x = log10_step, y = density, color = factor(state, levels = 1:6)) +
  geom_line() +
  labs(color = "Speed regime", y = "Empirical density", 
       x = expression(log[10]~"of step length")) +
  theme(legend.position = "inside",
        text = element_text(family = "LM Roman 10"),
        legend.position.inside = c(0.8, 0.6),
        legend.background = element_blank()) +
  guides(color = guide_legend(ncol = 2))

if(exportFig)
  ggsave("../../Figures/gannet-post-state-dist-K6.png",
         plot = final_speed_plot,
         width = 2 * graph_width,
         height = graph_height,
         units = graph_unit)


model_criterions <- read.table("results_gannets_model_selection_criterions.txt",
                               sep = ";", header = TRUE)

model_criterions_plot <- model_criterions %>% 
  mutate(value = -0.5 * value) %>% 
  mutate(criterion = ifelse(criterion == "-2 Loglik.", "Log-lik.", criterion)) %>%  
  ggplot() +
  aes(x = K, y = value, color = criterion) + 
  geom_point(aes(shape = criterion)) + 
  geom_line() + 
  labs(x = "Number of hidden states", y = "Criterion value",
       color = "") +
  scale_color_viridis_d(option = "C") 

if(exportFig)
  ggsave("../../Figures/gannet-model-criterions.png",
         plot = model_criterions_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)

illustrative_model <- readRDS("results_gannets_fitted_model.rds")

track_lines <- donnees %>% 
  mutate(
    viter = illustrative_model$best_seq,
    geometry_lagged = lead(geometry)
  ) %>%
  slice(-n()) %>% 
  # drop the NA row created by lagging
  mutate(
    line = st_sfc(purrr::map2(
      .x = geometry, 
      .y = geometry_lagged, 
      .f = ~{st_union(c(.x, .y)) %>% st_cast("LINESTRING")}
    ), crs = st_crs(.)))

trajectory_plot_with_states <- ggplot(zone_map) + 
  geom_sf(fill = "darkgreen") +
  theme_bw() +
  geom_sf(data = track_lines, aes(geometry = line, color = factor(viter))) + 
  geom_sf(data = donnees %>% 
            mutate(viter = illustrative_model$best_seq %>% factor()), 
          mapping = aes(color = viter)) + 
  theme(legend.position = "inside", 
        legend.direction = "horizontal", legend.background = element_blank(),
        legend.position.inside = c(0.4, 0.9),
        text = element_text(family = "LM Roman 10")) +
  coord_sf(xlim = c(-56.3, -56),
           ylim = c(46.75, 46.95)) + 
  # scale_color_viridis_d(option = "C") +
  # scale_color_manual(values = ten_color_palette[1:6]) + 
  labs(color = "Speed regime") 

if(exportFig)
  ggsave("../../Figures/gannet-trajectory-plot-with-states.png",
         plot = trajectory_plot_with_states,
         width = graph_width,
         height = graph_height,
         units = graph_unit)

transition_matrix_plot <- ggcorrplot(orig_mat[, 6:1], 
                                     lab = TRUE, 
                                     type = "full", lab_col = "red") +
  theme(legend.position = "none") +
  scale_fill_viridis_c()

if(exportFig)
  ggsave("../../Figures/gannet-transition-mat.png",
         plot = transition_matrix_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)
orig_mat <- illustrative_model$P
colnames(orig_mat) <- rownames(orig_mat) <- paste0("S", 1:ncol(orig_mat))
colnames(my_mat)
rownames(my_mat) <- nrow(illustrative_model$P):1


corrplot::corrplot(illustrative_model$P, col.lim = c(0, 1))
ggcorrplot(illustrative_model$P) +
  scale_fill_gradientn(colors = viridis(256, option = 'D'))
