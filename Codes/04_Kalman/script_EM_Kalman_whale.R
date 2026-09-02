rm(list = ls())
library(tidyverse)
library(sf)
source("utils_kalman_functions.R")


# Test on whale data ------------------------------------------------------

my_data <- read.table("../../Data/04_Kalman/whale.txt",
  sep = ";",
  colClasses = c(date = "POSIXct"),
  header = TRUE
) |>
  st_as_sf(crs = 4326, coords = c("X", "Y")) %>%
  mutate(
    long_m = st_coordinates(.)[, 1],
    lati_m = st_coordinates(.)[, 2]
  )

Y <- st_transform(my_data, crs = 32624) |>
  st_coordinates() %>%
  {
    . / 1000
  } |>
  scale(scale = FALSE)
tYY <- crossprod(Y)

Sigma0 <- diag(5, 2)
Omega0 <- diag(3, 4)
A_fixed <- cbind(diag(1, 2), matrix(0, 2, 2))
B_fixed <- diag(1, 4) + cbind(matrix(0, 4, 2), rbind(diag(4, 2), matrix(0, 2, 2)))
mu0 <- c(Y[1, ], apply(Y, 2, diff)[1, ] / 5)
init_pars <- list(
  A = A_fixed,
  B = B_fixed,
  Sigma = Sigma0,
  Omega = Omega0,
  mu0 = mu0,
  V0 = Omega0
)

n_steps <- 1000
Q_old <- rep(NA, n_steps)
Q_new <- rep(NA, n_steps)
ll <- rep(NA, n_steps + 1)
est_pars <- list(init_pars)
for (i in 1:n_steps) {
  filt_res <- get_filter_quants(pars_h = est_pars[[i]], Y = Y)
  ll[i] <- filt_res$loglik
  smooth_res <- get_smoother_quants(
    pars_h = est_pars[[i]],
    filt_pars = filt_res
  )
  Q_old[i] <- get_Q(
    pars = est_pars[[i]],
    smooth_pars = smooth_res,
    Y = Y
  )
  est_pars[[i + 1]] <- get_new_pars(smooth_res, Y,
    tYY = tYY,
    A = init_pars$A,
    B = init_pars$B
  )
  cat("\nIteration", i, "\n")
  Q_new[i] <- get_Q(
    pars = est_pars[[i + 1]],
    smooth_pars = smooth_res,
    Y = Y
  )
  if (i == n_steps) {
    final_filt <- get_filter_quants(est_pars[[n_steps + 1]], Y)
    ll[n_steps + 1] <- final_filt$loglik
    final_smooth <- get_smoother_quants(
      pars_h = est_pars[[n_steps + 1]],
      filt_pars = final_filt
    )
  }
}
Q_diff <- Q_new - Q_old
save(est_pars, ll, Q_diff, Y, final_filt, final_smooth, n_steps, my_data,
  file = "res_Kalman_whale.RData"
)
