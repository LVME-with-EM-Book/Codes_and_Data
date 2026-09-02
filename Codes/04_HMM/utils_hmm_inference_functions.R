compute_filt_prob <- function(params, ys){
  K <- length(params$init_distrib)
  densities <- sapply(ys, function(y){
    dnorm(y, params$means, params$sds)
  })
  n <- length(ys)
  filt_probs <- matrix(ncol = n, nrow = K)
  unnorm_probs <- params$init_distrib * densities[, 1]
  S <- sum(unnorm_probs)
  c0 <- log(S)
  ll <- c0
  filt_probs[,1] <- unnorm_probs / S
  for(t in 2:n){
    unnorm_probs <- (t(params$P) %*% filt_probs[,t-1]) * densities[, t]
    S <- sum(unnorm_probs)
    ct <- log(S)
    ll <- ll + ct
    filt_probs[, t] <- unnorm_probs / S
  }
  return(list(filt_probs = filt_probs, ll = ll))
}



get_logsumexp <- function(log_dens){
  max_log <- max(log_dens)
  max_log + log(sum(exp(log_dens - max_log)))
}



compute_alphas <- function(params, ys){
  K <- length(params$init_distrib)
  densities <- sapply(ys, function(y){
    dnorm(y, params$means, params$sds)
  })
  n <- length(ys)
  alphas <- matrix(ncol = n, nrow = K)
  alphas[, 1] <- params$init_distrib * densities[, 1]
  for(t in 2:n){
    alphas[, t] <- (t(params$P) %*% alphas[,t-1]) * densities[, t]
  }
  alphas
}

compute_betas <- function(params, ys){
  K <- length(params$init_distrib)
  densities <- sapply(ys, function(y){
    dnorm(y, params$means, params$sds)
  })
  n <- length(ys)
  betas <- matrix(ncol = n, nrow = K)
  betas[, n] <- 1
  for(t in 2:n){
    betas[, n - t + 1] <- params$P %*% (densities[, n - t + 2] *  betas[, n - t + 2])
  }
  betas
}

compute_xis_and_gammas <- function(params, ys){
  log_alphas <- log(compute_filt_prob(params, ys)$filt_probs)
  log_betas <- compute_log_betas(params, ys)
  log_alph_bet <- log_alphas + log_betas
  log_nrm_csts <- apply(log_alph_bet, 2, get_logsumexp)
  gammas <- exp(t(log_alph_bet) - log_nrm_csts)
  xis <- lapply(1:(n - 1),
      function(t){
        exp(logP + outer(log_alphas[, t],
                         log_densities[, t + 1] +  log_betas[, t + 1],
                         "+") - log_nrm_csts[t])
      }) %>% 
    do.call(what = function(...) abind::abind(..., along = 3))
  list(xis = xis, gammas = gammas)
}

compute_log_alphas <- function(params, ys){
  K <- length(params$init_distrib)
  log_densities <- sapply(ys, function(y){
    dnorm(y, params$means, params$sds, log = TRUE)
  })
  logP <- log(params$P)
  n <- length(ys)
  log_alphas <- log_betas <- matrix(ncol = n, nrow = K)
  log_alphas[, 1] <- log(params$init_distrib) + log_densities[, 1]
  for (t in 2:n) {
    ms_alphas <- apply(logP + log_alphas[, t - 1], 2, max)
    log_alphas[, t] <- log_densities[, t] +  ms_alphas + 
      log(colSums(exp(logP + log_alphas[, t - 1] - 
                        matrix(ms_alphas, nrow = K, ncol = K, byrow = TRUE))))
    
  }
  log_alphas
}

compute_log_betas <- function(params, ys){
  K <- length(params$init_distrib)
  log_densities <- sapply(ys, function(y){
    dnorm(y, params$means, params$sds, log = TRUE)
  })
  logP <- log(params$P)
  n <- length(ys)
  log_betas <- matrix(0, ncol = n, nrow = K)
  for(t in 2:n){
    ms_betas <- apply(logP, 1, function(x) max(x + log_densities[, n - t + 2] +  log_betas[, n - t + 2]))
    log_betas[, n - t + 1] <- ms_betas + log(colSums(exp(t(logP) + 
                                                           log_densities[, n - t + 2] +  
                                                           log_betas[, n - t + 2] -
                                                           matrix(ms_betas, nrow = K, ncol = K, byrow = TRUE))))
  }
  log_betas
}





compute_log_alphas_betas <- function(params, ys){
  K <- length(params$init_distrib)
  log_densities <- sapply(ys, function(y){
    dnorm(y, params$means, params$sds, log = TRUE)
  })
  logP <- log(params$P)
  n <- length(ys)
  log_alphas <- log_betas <- matrix(ncol = n, nrow = K)
  log_betas[, n] <- 0
  log_alphas[, 1] <- log(params$init_distrib) + log_densities[, 1]
  for(t in 2:n){
    ms_alphas <- apply(logP + log_alphas[, t - 1], 2, max)
    ms_betas <- apply(logP, 1, 
                      function(x) max(x + log_densities[, n - t + 2] +  
                                        log_betas[, n - t + 2]))
    log_alphas[, t] <- log_densities[, t] +  ms_alphas + 
      log(colSums(exp(logP + log_alphas[, t - 1] - 
                        matrix(ms_alphas, nrow = K, ncol = K, byrow = TRUE))))
    log_betas[, n - t + 1] <- ms_betas + log(colSums(exp(t(logP) + 
                                                           log_densities[, n - t + 2] +  
                                                           log_betas[, n - t + 2] -
                                                           matrix(ms_betas, nrow = K, 
                                                                  ncol = K, byrow = TRUE))))
  }
  log_lik <- max(log_alphas[, n]) + log(sum(exp(log_alphas[, n] - max(log_alphas[, n]))))
  
  list(log_alphas = log_alphas, 
       log_betas = log_betas,
       gammas = (log_alphas + log_betas) %>% 
         apply(2, function(x){
           exp(x - log_lik)
         }) %>% 
         t(),
       xis = map(1:(n - 1),
                 function(t){
                   exp(logP + outer(log_alphas[, t],
                                    log_densities[, t + 1] +  log_betas[, t + 1],
                                    "+") - log_lik)
                 }) %>% 
         do.call(what = function(...) abind::abind(..., along = 3)),
       log_lik = log_lik)
}

get_dS4_params <- function(mod_, order_mean = FALSE){
  K = length(mod_@prior@parameters$coefficients)
  ll <- logLik(mod_)
  vit_seq <- viterbi(mod_)[, "state"]
  old_params <- list(init_distrib =
         mod_@prior@parameters$coefficients %>% set_names(NULL),
       P = mod_@transition %>%
         map(function(x) x@parameters$coefficients) %>% do.call(what = rbind),
       means = getpars(mod_)[seq(from = K * (K + 1) + 1, by = 2, length.out = K)] %>%
         set_names(NULL),
       sds = getpars(mod_)[seq(from = K * (K + 1) + 2, by = 2, length.out = K)] %>%
         set_names(NULL),
       K = K,
       D = (2 * ll + AIC(ll)) * 0.5,
       best_seq = vit_seq)
  unnormed_statio_prob <- Re(eigen(t(old_params$P))$vectors[,1])
  old_params$statio_prob <- unnormed_statio_prob / sum(unnormed_statio_prob)
  if(order_mean){
    new_params <- old_params
    my_ord <- order(old_params$means)
    new_params$means <- old_params$means[my_ord]
    new_params$sds <- old_params$sds[my_ord]
    new_params$P <- old_params$P[my_ord, my_ord]
    new_params$statio_prob <- old_params$statio_prob[my_ord]
    new_params$init_distrib <- old_params$init_distrib[my_ord]
    new_params$best_seq <- factor(old_params$best_seq,
                                  levels = 1:K,
                                  labels = paste(rank(old_params$means))) %>% 
      as.character() %>% 
      as.numeric()
    return(new_params)
  }
  else{
    return(old_params)
  }
}

get_viterbi_ICL <- function(mod_, ys){
  n <- length(ys)
  vit_seq <- viterbi(mod_)[, "state"]
  params <- get_dS4_params(mod_)
  K <- length(params$init_distrib)
  log_init <- log(params$init_distrib[vit_seq[1]])
  log_obs_dens <- sum(dnorm(ys, params$means[vit_seq], params$sds[vit_seq], log = TRUE))
  log_trans <- sapply(2:n, function(i){
    log(params$P[vit_seq[i - 1], vit_seq[i]])
  }) %>% 
    sum()
  -2 * (log_init + log_obs_dens + log_trans) + params$D * log(n)
}

get_ICL_penalty <- function(m){
  frwd_bckwd_quants <- forwardbackward(m, return.all=TRUE)
  apply(frwd_bckwd_quants$gamma, 1, function(x){
    -sapply(x, function(z) ifelse(z == 0, 0, z * log(z)))
  }) %>%
    sum()
}

get_hclust_start <- function(obs, K){
  hclust_data <- dist(obs) %>% 
    hclust(method = "ward.D2")
  labels <- cutree(hclust_data, k = K)
  lst_guess <- split(obs, labels) %>% 
    lapply(function(x) c("(Intercept)" = mean(x), "sd" = sd(x))) %>% 
    set_names(NULL)
  ord <- order(map_vec(lst_guess, "(Intercept)"))
  unlist(lst_guess[ord])
}