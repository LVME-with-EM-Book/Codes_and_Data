rm(list = ls())

library(tidyverse)
library(sf)
library(geodata)
source("../utils_figure_parameters.R")
source("../utils_ggplot_theme.R")
exportFig <- FALSE

if(!dir.exists("world_map")){
  dir.create("world_map")
}
# Get world land polygons (country boundaries)
world <- geodata::world(path = "world_map/")

# Convert to sf object if not already
world_sf <- st_as_sf(world)

my_data <- read.table("whale.txt", sep = ";", 
                      colClasses = c(date = "POSIXct"),
                      header = TRUE) %>% 
  rename(longitude = X, latitude = Y)
my_data_sf <- my_data  |> 
  dplyr::select(longitude, latitude) %>% 
  as.matrix() %>% 
  st_linestring() %>% 
  st_sfc() %>% 
  st_sf(geometry = ., crs = 4326)

# Plot the map
map_plot <- ggplot() +
  geom_sf(data = world_sf, fill = "darkgreen", color = "black") +
  geom_sf(data = my_data_sf) +
  coord_sf(xlim = c(-60, -25), ylim = c(35, 70), expand = FALSE) +
  labs(x = "", y = "")

Sys.setlocale("LC_TIME", "C")  # or "en_US.UTF-8" on some systems
coord_plot <- pivot_longer(my_data, cols = c("longitude", "latitude"),
               names_to = "Coordinate", values_to = "y", names_prefix = "location.", 
               names_transform = str_to_sentence) %>% 
  ggplot(aes(x = date, y = y)) +
  labs(x = "", y = "") +
  facet_wrap(Coordinate~., scales = "free", nrow = 2, strip.position = "left") +
  geom_line() +
  scale_x_datetime(breaks = seq.POSIXt(from = as.POSIXct("2009-05-01"), 
                                       to = as.POSIXct("2009-06-30"),
                                       length.out = 5),
                   date_labels = "%b %d")
trajectory_plot <- gridExtra::grid.arrange(map_plot, coord_plot, nrow = 1)

if(exportFig)
  ggsave("../../Figures/whale-trajectory-plot.png",
       plot = trajectory_plot,
       width = 2 * graph_width,
       height = graph_height,
       units = graph_unit)


# Representing results ----------------------------------------------------

load("res_Kalman_whale.RData")

## Ellipse on trajectory ---------------------------------------------------

original_mean <- attr(Y, "scaled:center")
selected_indexes <- 10:20
ellipses_conf <- map_dfr(selected_indexes,
                         function(t)
                           ellipse::ellipse(x = final_smooth$smooth_vars[1:2, 1:2, t],
                                            level = 0.5,
                                            centre = final_smooth$smooth_means[t,1:2],) |> 
                           as.data.frame() |> 
                           mutate(t = t)) 
ellipses <- ellipses_conf |> 
  mutate(x = (x + original_mean[1]) * 1000,
         y = (y + original_mean[2]) * 1000) |> 
  st_as_sf(coords = c("x", "y"), crs = 32624) |> 
  st_transform(crs = 4326) |> 
  st_coordinates() |> 
  as.data.frame() |> 
  mutate(t = ellipses_conf$t)
means <- final_smooth$smooth_means[selected_indexes, ] |> 
  as.data.frame() |> 
  rename(x = V1, y = V2) |> 
  mutate(x = (x + original_mean[1]) * 1000,
         y = (y + original_mean[2]) * 1000) |> 
  st_as_sf(coords = c("x", "y"), crs = 32624) |> 
  st_transform(crs = 4326) |> 
  st_coordinates() |> 
  as.data.frame()
original_data <- st_coordinates(my_data) |> 
  as.data.frame()
reconstruction_plot <- ggplot(ellipses, aes(X, Y)) +
  geom_polygon(aes(group = t), fill = "darkred", alpha = .5) +
  geom_path(data = means, color = "darkred", linetype = 2) +
  geom_point(data = means, color = "darkred") +
  geom_point(data = original_data[selected_indexes, ]) +
  geom_path(data = original_data[selected_indexes, ], linetype = 3) +
  labs(x = "Longitude", y = "Latitude")

if(exportFig)
ggsave("../../Figures/whale-reconstruction-plot.png",
       plot = reconstruction_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)

# Methodes de Davies (1980) pour calcul des quantiles
qXtX <- function(mu, Sigma, p = c(.025, 0.5, .975)) {
  eig <- eigen(Sigma, symmetric = TRUE)
  lambda <- eig$values
  m <- crossprod(eig$vectors, mu)
  # Si toutes les valeurs propres sont distinctes :
  h <- rep(1, length(lambda))
  delta <- m^2 / lambda
  # P(Q > q)
  surv_fun <- function(q) {
    CompQuadForm::davies(q, lambda, h, delta)$Qq
  }
  # Quantile par inversion de P(Q > q) = 1-p
  qfun <- function(p) {
    upper <- sum(lambda * (1 + delta)) * 2
    while (surv_fun(upper) > 1 - p)
      upper <- 2 * upper
    uniroot(
      function(q) surv_fun(q) - (1 - p),
      c(0, upper)
    )$root |> 
      sqrt()
  }
  sapply(p, qfun) |> 
    set_names(p)
}


# Figure for speed --------------------------------------------------------


all_quantiles <- sapply(1:nrow(Y), function(t)
  qXtX(mu = final_smooth$smooth_means[t, 3:4],
       Sigma = final_smooth$smooth_vars[3:4, 3:4, t])) |> 
  t() |> 
  as.data.frame() |> 
  rename_all(function(x) paste0("q", x)) |> 
  mutate(date = my_data$date)
speed_plot <- ggplot(all_quantiles) + 
  aes(x = date) + 
  geom_line(aes(y = q0.5)) + 
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), 
              alpha = 0.5, fill = "lightblue") +
  labs(y = "Estimated mean speed (km/h)",
       x = "Date")

if(exportFig)
ggsave("../../Figures/whale-speed-plot.png",
       plot = speed_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)
