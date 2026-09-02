library(tidyverse)
# Data from this article
# Assessing Performance of Bayesian State-Space Models Fit to Argos Satellite Telemetry Locations Processed with Kalman Filtering
# Extracted from movebank
# https://www.movebank.org/cms/webapp?gwt_fragment=page%3Dstudies%2Cpath%3Dstudy72289508
whales <- read.csv("whales.csv")
whales_with_traj <- whales %>% 
  dplyr::select(location.long, location.lat, timestamp, tag.local.identifier) %>% 
  na.omit() %>% 
  unique() %>% 
  mutate(date = as_datetime(timestamp)) %>% 
  group_by(tag.local.identifier) %>% 
  arrange(date) %>% 
  mutate(time_lag = c(4, as.numeric(difftime(date, lag(date), units = "hours"))[-1])) %>%
  filter(time_lag > 2) %>% 
  mutate(traj = cumsum(time_lag > 6)) %>% 
  ungroup()
# Checking number of points per traj
group_by(whales_with_traj, traj, tag.local.identifier) %>%
  count() %>% ungroup() %>%
  arrange(desc(n))
my_tag <- "80711"
my_traj <- 0
my_data <- filter(whales_with_traj, tag.local.identifier == my_tag, traj == my_traj) %>% 
  rename(longitude = location.long,
         latitude = location.lat) %>% 
  dplyr::select(date, longitude, latitude)
write.table(my_data, file = "whale_trajectory.txt",
            sep = ";", col.names = TRUE, row.names = FALSE)



# add noise ---------------------------------------------------------------


set.seed(123)
whale_sf <- read.table("whale_trajectory.txt", sep = ";", 
                       colClasses = c(date = "POSIXct"),
                       header = TRUE) |> 
  st_as_sf(crs = 4326, coords = c("longitude", "latitude"))
whale_noise <- whale_sf |> 
  st_transform(crs = 32624) |> 
  st_coordinates() %>%
  {. + rnorm(nrow(.) * 2, sd = 10000)} |>
  as.data.frame() |> 
  st_as_sf(crs = 32624, coords = c("X", "Y")) |> 
  st_transform(crs = 4326) |> 
  st_coordinates() |> 
  as.data.frame() |> 
  mutate(date = whale_sf$date)
write.table(whale_noise, file = "whale.txt", sep = ";", row.names = FALSE,
            col.names = TRUE)
