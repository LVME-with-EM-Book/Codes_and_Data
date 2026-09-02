require(tidyverse)
require(abind)
require(mixtools)
require(parallel)
source("utils_gaussian_mixture_functions.R")

ys <- read.table("../../Data/03_GMM/Bohemia_vegetation.csv",
                 header = TRUE, row.names = 1, sep = ",") %>% 
  as.matrix()

get_best_estimate <- function(ks, ys, n_trys){
  purrr::map(ks, function(k){
    if(k == 1){
      all_res_k <- get_EM_estimation_mixture(ys = ys, 
                                theta0 = list(pi = 1,
                                              mu = rep(0, ncol(ys)),
                                              Sigma = array(diag(1, ncol(ys)),
                                                            dim = c(ncol(ys), ncol(ys), 1))), 
                                n_iter_max = 2000)
      return(c(all_res_k, k = k))
    }
    else{
      all_res_k <- mclapply(1:n_trys, function(i){
        set.seed(i)
        initial_theta <- get_theta0(k = k, ys = ys)
        all_res <- try(get_EM_estimation_mixture(ys = ys, 
                                                 theta0 = initial_theta, 
                                                 n_iter_max = 2000))
        if(inherits(all_res, "try-error")){
          return(NULL)
        }
        else{
          return(c(all_res, seed = i))
        }
      }, mc.cores = detectCores() - 2)
      lls <- map_dbl(all_res_k, function(x){
        if(is.null(x)){
          return(-Inf)
        }
        else{
          return(x$final_criterions$log_lik)
        }})
      return(c(all_res_k[[which.max(lls)]], k = k))
    }
  })
}


gmm_est_multiple_ks <- get_best_estimate(1:6, ys, 200)
gmm_est_k3 <- mclapply(1:200, function(i){
  set.seed(i)
  initial_theta <- get_theta0(k = 3, ys = ys)
  all_res <- try(get_EM_estimation_mixture(ys = ys, 
                                           theta0 = initial_theta, 
                                           n_iter_max = 2000))
  if(inherits(all_res, "try-error")){
    return(NULL)
  }
  else{
    return(c(all_res, seed = i))
  }
}, mc.cores = detectCores() - 2)
gmm_est_best <- map_dbl(gmm_est_k3, function(x){
  if(is.null(x))
    return(-Inf)
  else
    return(x$final_criterions$log_lik)
}) %>% 
  which.max() %>% 
  gmm_est_k3[[.]]

save(list = ls(pattern = "gmm_est"),
     file = "gmm_est.RData")
rm(list = ls(pattern = "get_"), simulate_gaussian_mixture, ys)
