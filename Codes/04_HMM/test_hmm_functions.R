library(tidyverse)
K <- 3
set.seed(123)
params_true <- list(init_distrib = rep(1 / K, K),
                   P =  lapply(1:K, function(i){
                     alphs <- rep(1, K)
                     alphs[i] <- 10
                     dirmult::rdirichlet(1, alphs)[1, ]
                   }) %>% 
                     do.call(what = rbind),
                   means = seq(from = -K, to = K, length.out = K) %>% floor(), 
                   sds = seq(0.1, 1, length.out = K) %>% round(2))
generate_hmm <- function(params, n_obs){
  xs <- rep(NA, n_obs)
  K <- length(params$means)
  xs[1] <- sample(1:K, prob = params$init_distrib, size = 1)
  for(t in 2:n_obs){
    xs[t] <- sample(1:K,
                 prob = params$P[xs[t - 1], ],
                 size = 1)
  }
  ys <- rnorm(n_obs, mean = params$means[xs],
              sd = params$sds[xs])
  data.frame(x = factor(xs), y = ys, t = 1:n_obs)
}

