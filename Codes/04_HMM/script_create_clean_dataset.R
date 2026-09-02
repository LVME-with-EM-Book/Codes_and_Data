library(tidyverse)
library(sf)
library(move2)
library(geodata)


# Si filtre 1h
# Prendre ID-9_T-189
# Si filtre 5h
# Prendre
donnees_completes <- read.table("../../Data/04_HMM/raw_gannet.csv", 
                                header = TRUE, row.names = 1, sep = ",") %>% 
  filter(!is.na(location.long), !is.na(location.lat)) %>% 
  mutate(timestamp = lubridate::as_datetime(timestamp)) %>% 
  arrange(individual.local.identifier,
          timestamp) %>%
  mutate(BIRD_ID = factor(individual.local.identifier) %>%
           as.numeric() %>%
           paste0("ID-", .)) %>%
  group_by(BIRD_ID) %>%
  mutate(dt = difftime(timestamp, lag(timestamp), units = "hours") %>% 
           as.numeric %>% replace_na(0),
         traj = paste0("T-", cumsum(dt > 5) + 1)) %>% 
  ungroup() %>% 
  mutate(ID = paste(BIRD_ID, traj, sep = "_")) %>% 
  dplyr::select(timestamp, location.long, location.lat, external.temperature, ID) %>% 
  rename(x = location.long, y = location.lat) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(crs = 6665) 

interpolate_values <- function(times, y_values) {
  # Convert times to numeric for interpolation
  numeric_times <- as.numeric(times)
  
  # Perform linear interpolation
  interpolated_y <- approx(x = numeric_times, 
                           y = y_values, 
                           xout = numeric_times, 
                           method = "linear", 
                           rule = 2)$y
  
  return(interpolated_y)
}


donnees <- filter(donnees_completes,
                  ID == "ID-9_T-63") %>% 
  move2::mt_as_move2(time_column = "timestamp", 
                     track_id_column = "ID") %>% 
  mt_interpolate(time = seq.POSIXt(.$timestamp[1], .$timestamp[nrow(.)], 
                                   by = "15 mins"), omit = FALSE) %>% 
  mutate(interpolated = is.na(external.temperature),
         temperature = interpolate_values(timestamp, 
                                       y = external.temperature)) %>% 
  filter(timestamp %in% seq.POSIXt(.$timestamp[1], .$timestamp[nrow(.)], 
                                   by = "15 mins")) %>% 
  dplyr::select(-external.temperature, -interpolated) %>% 
  mutate(x_m = st_coordinates(.)[, 1],
         y_m = st_coordinates(.)[, 2],
         long = st_transform(., crs = 4326) %>% 
           st_coordinates() %>% 
           {.[, 1]},
         lat = st_transform(., crs = 4326) %>% 
           st_coordinates() %>% 
           {.[, 2]}) %>% 
  st_drop_geometry() %>% 
  as.data.frame() %>% 
  moveHMM::prepData(coordNames = c("x_m", "y_m"), type = "UTM") %>% 
  replace_na(list(angle = 0)) %>% 
  mutate(Vx = step * cos(angle),
         Vy = step * sin(angle)) %>% 
  as.data.frame() %>% 
  na.omit() %>% 
  mutate(log10_step = log(step, base = 10)) %>% 
  rename(t = timestamp,
         Longitude = long,
         Latitude = lat)

rm(donnees_completes)
donnees %>% 
  dplyr::select(t, Longitude, Latitude, log10_step, temperature) %>% 
  write.table(file = "../../Data/04_HMM/gannet_trajectory.csv", sep = ",", 
              col.names = TRUE, row.names = FALSE)
