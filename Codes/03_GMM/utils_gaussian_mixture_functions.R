require(mixtools)
# SIMULATION --------------------------------------------------------------

# Fonction principale -----------------------------------------------------

simulate_gaussian_mixture <- function(theta, n){
  
  K <- length(theta$pi) # Nombre de classes (se déduit des paramètres)
  
  # Tirage des variables latentes
  replicate(n,
            {
              z <- sample(x = 1:K, # Tirage dans le vecteur 1:K
                          size = 1, # Nombre de tirages
                          replace = TRUE, # Avec remplacement ou non (par défaut, non)
                          prob = theta$pi)
              rmvnorm(n = 1, mu = theta$mu[z,], 
                      sigma = theta$Sigma[,, z])[1,]
            }, simplify = TRUE) %>% 
    t()
}


# ESTIMATION --------------------------------------------------------------

get_mixture_log_likelihood <- function(theta, ys){
  weights <- theta$pi
  means <- theta$mu
  Sigmas <- theta$Sigma
  # sapply permet d'appliquer à chaque élément d'un vecteur (ici, ys)
  # un même fonction, ici la fonction qui à y associe
  # log(sum(weights * dnorm(y, means, standard_deviations))
  log_vraisemblance_individuelles <- apply(ys, 1, function(y){
    sapply(1:length(weights),
           function(k){
             weights[k] * mixtools::dmvnorm(y, means[k, ], Sigmas[,,k])
           }) %>% 
      sum() %>% 
      log()
  })
  return(sum(log_vraisemblance_individuelles))
}

# Algorithme EM -----------------------------------------------------------

get_theta0 <- function(k, ys, method = "pure_random"){
  if(method == "kmeans"){
    km_clust <- kmeans(ys, centers = k, nstart = 20)$cluster
    list(pi = as.numeric(table(km_clust) / nrow(ys)),
         mu = sapply(1:k, function(j){
           colMeans(subset(ys, km_clust == j))}) %>% 
           t(),
         Sigma =  lapply(1:k, function(j){
           subset(ys, km_clust == j) %>% 
             cov()  
         }) %>% 
           do.call(what = function(...) abind(..., along = 3))) %>% 
      return()
  }
  else if(method == "pure_random"){
    accepted <- FALSE
    while(!accepted){
      cand <- list(pi = rep(1/k, k),
                   mu = apply(ys, 2, function(z){
                     ymin <- min(z)
                     ymax = max(z)
                     runif(k, ymin, ymax)
                   }),
                   Sigma =  lapply(1:k, function(j){
                     cov(ys) / k  
                   }) %>% 
                     do.call(what = function(...) abind(..., along = 3)))
      accepted <- all(colSums(get_e_step_mixture(cand, ys)) > 3)
    }
    return(cand)
  }
}

get_e_step_mixture <- function(theta, ys){
  dim_y <- ncol(ys)
  n <- nrow(ys)
  weights <- theta$pi
  K <- length(weights)
  means <- matrix(theta$mu, nrow = K)
  Sigmas <- array(theta$Sigma, dim = c(dim_y, dim_y, K))
  # E step consists in computing matrix of class posterior probability 
  # for each class and each observations.
  # Therefore, it will return a matrix of size n times K
  # The sapply function returns a matrix of size K times, n
  post_weights <- apply(ys, 
                        1, # Pour chaque ys
                        function(y){ # On applique la fonction
                          # On clacule les K poids non normalisés
                          unnormed_log_w <- sapply(1:K,
                                                         function(k){
                                                           log(weights[k]) + mixtools::logdmvnorm(y, 
                                                                                                  means[k, ], 
                                                                                                  Sigmas[,,k])
                                                         })
                          max_uw <- max(unnormed_log_w)
                          log_nrm_cst <- max_uw + log(sum(exp(unnormed_log_w - max_uw)))
                          # Normalization
                          post_weight <- exp(unnormed_log_w - log_nrm_cst)
                        })
  # Le renvoi est sour la forme d'une matricce K*n
  return(t(post_weights)) # Donc on renvoit la transposée
}


get_q_theta_theta0 <- function(theta, theta0, ys){
  post_weights <- get_e_step_mixture(theta0, ys)
  K <- length(theta$pi)
  sum(sapply(1:K,
             function(k){
               sum(post_weights[, k] * (log(theta$pi[k]) + 
                                          mixtools::logdmvnorm(ys, 
                                                               theta$mu[k, ], 
                                                               theta$Sigma[,,k])))
             }))
}

get_m_step_mixture <- function(posterior_weights, ys){
  dim_y <- ncol(ys)
  n <- nrow(ys)
  K <- ncol(posterior_weights)
  Nks <- colSums(posterior_weights)
  # Actualisation des poids
  updated_pi <- Nks / n
  # Actualisation des mu
  updated_mu <- sapply(1:K, 
                       function(k){
                         colSums(posterior_weights[, k] * ys) / Nks[k]
                       }) %>% 
    t()
  # Matrice des ys centrées
  updated_Sigma <- lapply(1:K, 
                          function(k){
                            m <- updated_mu[k, ]
                            N <- Nks[k]
                            lapply(1:n, function(i){
                              posterior_weights[i, k] * (ys[i, ] - m) %*% t(ys[i, ] - m) / N
                            }) %>% 
                              Reduce(f = "+")
                          }) %>% 
    do.call(what = function(...) abind(..., along = 3))
  # Nouveau paramètres
  new_pars <- list(pi = updated_pi,
                   mu = updated_mu,
                   Sigma = updated_Sigma)
  return(new_pars)
}

get_EM_estimation_mixture <- function(ys, theta0, n_iter_max, all_params = FALSE){
  # Initialisation des paramètres et des quantités annexes
  K <- length(theta0$pi)
  dim_y <- ncol(ys)
  n <- nrow(ys)
  if(K == 1) {
    all_params <- FALSE
    new_theta <- list(pi = 1,
                      mu = matrix(colMeans(ys), nrow = 1),
                      Sigma = array(cov(ys)  * (n - 1) / n, 
                                    dim = c(dim_y, dim_y, 1)))
    log_likelihood = NULL
    new_weights = matrix(1, ncol = 1, nrow = n)
    get_mixture_log_likelihood(new_theta, ys)
    
  }
  else {
    params_list <- list(theta0)
    log_likelihood <- get_mixture_log_likelihood(theta0, ys) / n # Log vraisemblance
    # weights_list <- list() # Poids a posterior
    stop_criterion <- FALSE
    n_iter <- 0
    while(!stop_criterion){
      n_iter <- n_iter + 1
      # Etape E
      new_weights <- get_e_step_mixture(theta0, ys)
      # Etape M
      new_theta <- get_m_step_mixture(posterior_weights = new_weights, ys = ys)
      theta0 <- new_theta
      # Stockages des quantités annexes
      params_list <- c(params_list, list(new_theta))
      # weights_list <- c(weights_list, list(new_weigths))
      new_log_likelihood <- get_mixture_log_likelihood(new_theta, ys) / n
      log_likelihood <- c(log_likelihood, new_log_likelihood)
      # Critères d'arrêt
      likelihood_criterion <- (new_log_likelihood - log_likelihood[n_iter]) < 10^(-8)
      iteration_criterion <- (n_iter >= n_iter_max)
      stop_criterion <- (likelihood_criterion | iteration_criterion) 
    }
    # Dernières quantités à stocker
    new_weights <- get_e_step_mixture(theta0, ys)
    # weights_list <- c(weights_list, list(new_weigths))
    # Message sur la convergence
    if(likelihood_criterion){
      message(paste("Convergence was reached after", n_iter, "iterations."))
    }
    else{
      message(paste("Algorithm was stopped after", n_iter, "iterations."))
    }
  }
  result <- list( 
              posterior_weights = new_weights,
              log_likelihoods = log_likelihood,
              final_param = new_theta,
              final_criterions = data.frame(k = K,
                                            log_lik = get_mixture_log_likelihood(new_theta, ys),
                                            AIC = get_AIC(new_theta, ys),
                                            BIC = get_BIC(new_theta, ys),
                                            ICL = get_ICL(new_theta, ys)))
  if(all_params){
    result$parameter_sequence <- params_list
  }
  return(result)
}

get_AIC <- function(theta, ys){
  K <- length(theta$pi)
  dim_y <- ncol(ys)
  D_M <- (K - 1) + K * (dim_y + dim_y * 0.5 * (dim_y + 1))  
  n <- nrow(ys)
  -2 * get_mixture_log_likelihood(theta, ys) + 2 * D_M 
}

get_BIC <- function(theta, ys){
  K <- length(theta$pi)
  dim_y <- ncol(ys)
  D_M <- (K - 1) + K * (dim_y + dim_y * 0.5 * (dim_y + 1))  
  n <- nrow(ys)
  - 2 * get_mixture_log_likelihood(theta, ys) +  D_M * log(n) 
}

get_ICL <- function(theta, ys){
  K <- length(theta$pi)
  dim_y <- ncol(ys)
  D_M <- (K - 1) + K * (dim_y + dim_y * 0.5 * (dim_y + 1))  
  n <- nrow(ys)
  -2 * get_q_theta_theta0(theta, theta, ys) + D_M * log(n) 
}

