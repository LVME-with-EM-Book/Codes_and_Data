rm(list = ls())
library(tidyverse)
library(depmixS4)
source("utils_hmm_inference_functions.R")

# Fit on gannets on data ---------------------------------------------------------

gannets <- read.table(file = "../../Data/04_HMM/gannet_trajectory.csv", 
                      sep = ",", 
                      header = TRUE)

set.seed(123) # Fore reproducibility

mods_wo_cov <- lapply(2:10, function(k){
  obs <- gannets$log10_step
  model_list <- parallel::mclapply(1:15, function(i){
    set.seed(i)
    depmix(log10_step ~ 1,
           data = gannets,
           nstates = k) %>%
      fit(emcontrol=em.control(maxit = 1000))
    },
    mc.cores = 15)
})

# mods_hclust <- lapply(2:10, function(k){
#   obs <- gannets$log10_step
#   model_list <- parallel::mclapply(1:1, function(i){
#     set.seed(i)
#     depmix(log10_step ~ 1,
#            data = gannets,
#            respstart = get_hclust_start(gannets$log10_step, K = k),
#            nstates = k) %>%
#       fit(emcontrol=em.control(maxit = 1000))
#   },
#   mc.cores = 1)
# })
# 
# mods_with_cov <- lapply(2:10, function(k){
#   model_list <- parallel::mclapply(1:15, function(i){
#     set.seed(i)
#     depmix(log10_step ~ 1,
#            nstates = k,
#            data = gannets,
#            transition = ~scale(temperature)) %>%
#       fit(emcontrol=em.control(maxit = 1000))
#   },
#   mc.cores = 15)
# })



all_crits_wo_cov <- lapply(mods_wo_cov, 
       function(m) sapply(m, function(x){
         c(ll = logLik(x), AIC = AIC(x), BIC = BIC(x), ICL_1 = BIC(x) + get_ICL_penalty(x),
           ICL_2 = get_viterbi_ICL(x, gannets$log10_step))
       }))

# all_crits_with_cov <- lapply(mods_with_cov, 
#                            function(m) sapply(m, function(x){
#                              c(ll = logLik(x), AIC = AIC(x), BIC = BIC(x), ICL_1 = BIC(x) + get_ICL_penalty(x),
#                                ICL_2 = get_viterbi_ICL(x, gannets$log10_step))
#                            }))
# 
# all_crits_hclust <- lapply(mods_hclust, 
#                              function(m) sapply(m, function(x){
#                                c(ll = logLik(x), AIC = AIC(x), BIC = BIC(x), ICL_1 = BIC(x) + get_ICL_penalty(x),
#                                  ICL_2 = get_viterbi_ICL(x, gannets$log10_step))
#                              }))

best_mods_wo_cov <- mapply(x = all_crits_wo_cov, y = mods_wo_cov, FUN = function(x, y){
  i <- which.max(x["ll", ])
  y[[i]]
})

best_crits_wo_cov <- map_dfr(best_mods_wo_cov,  function(x){
                             c(ll = -2 * logLik(x), AIC = AIC(x), 
                               BIC = BIC(x), 
                               ICL = BIC(x) + get_ICL_penalty(x)) %>% 
    enframe(name = "criterion") %>% 
    mutate(K = x@nstates,
           criterion = factor(criterion,
                              levels = c("ll", "AIC", "BIC", "ICL"),
                              labels = c("-2 Loglik.", "AIC", "BIC", "ICL")))
                           })

write.table(best_crits_wo_cov, file = "results_gannets_model_selection_criterions.txt",
            sep = ";", col.names = TRUE)

# all_post_densities_old <- map_dfr(best_mods_wo_cov, function(m){
#   new_params <- get_dS4_params(m)
#   my_ord <- order(old_params$means)
#   new_params$means <- old_params$means[my_ord]
#   new_params$sds <- old_params$sds[my_ord]
#   new_params$P <- old_params$P[my_ord, my_ord]
#   unnormed_statio_prob <- Re(eigen(t(new_params$P))$vectors[,1])
#   new_params$statio_prob <- unnormed_statio_prob / sum(unnormed_statio_prob)
#   new_params$init_distrib <- old_params$init_distrib[my_ord]
#   res <- pmap_dfr(new_params[c("means", "sds", "statio_prob")], function(means, sds, statio_prob){
#     x <- seq(0.75, 4.2, length.out = 501)
#     data.frame(log10_step = x, 
#                density = statio_prob * dnorm(x, means, sds))
#   }, .id = "state") %>% 
#     mutate(n_states = m@nstates)
# })
all_post_densities <- map_dfr(best_mods_wo_cov, function(m){
  new_params <- get_dS4_params(m, order_mean = TRUE)
  res <- pmap_dfr(new_params[c("means", "sds", "statio_prob")], function(means, sds, statio_prob){
    x <- seq(0.75, 4.2, length.out = 501)
    data.frame(log10_step = x, 
               density = statio_prob * dnorm(x, means, sds))
  }, .id = "state") %>% 
    mutate(n_states = m@nstates)
})


write.table(all_post_densities, file = "results_gannets_post_densities.txt",
            sep = ";", col.names = TRUE)

illustrative_model <- best_mods_wo_cov[[5]] %>% 
  get_dS4_params(order_mean = TRUE)
saveRDS(illustrative_model, file = "results_gannets_fitted_model.rds")


