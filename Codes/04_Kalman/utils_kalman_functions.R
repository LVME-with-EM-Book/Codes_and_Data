get_logdnorm <- function(y, m, V){
  z <- y - m
  out <- -0.5*(2 * log(2 * pi) + determinant(V, logarithm = TRUE)$modulus +
                 sum(z * solve(V, z)))
  as.numeric(out)
}
get_filter_quants <- function(pars_h, Y){
  # Model infos
  dz <- ncol(pars_h$B)
  dy <- ncol(Y)
  n_obs <- nrow(Y)
  Iz <- diag(1, dz)
  # Creation of output objects
  filt_means <- matrix(nrow = n_obs, 
                       ncol = dz)
  filt_vars <- array(dim = c(dz, dz, n_obs))
  pred_vars <- array(dim = c(dz, dz, n_obs)) # W_t
  loglik <- 0
  # Parameter extraction
  A_h <- pars_h$A; B_h <- pars_h$B; 
  Sigma_h <- pars_h$Sigma; Omega_h <- pars_h$Omega
  V0 <- pars_h$V0
  mu0 <- pars_h$mu0
  # Starting the filter
  ## Initialization
  mean_Y0 <- as.numeric(A_h %*% mu0)
  var_Y0 <- Sigma_h + A_h %*% V0 %*% t(A_h)
  Kalman_gain <- V0 %*% t(solve(var_Y0, A_h))
  filt_means[1, ] <- mu0 + Kalman_gain %*% (Y[1, ] - mean_Y0)
  filt_vars[,, 1] <- (Iz - Kalman_gain %*% A_h) %*% V0
  # Initialization of likelihood
  all_ll <- rep(NA, nrow(Y))
  all_ll[1] <- get_logdnorm(Y[1, ],
                          mean_Y0,
                          var_Y0)
  for(t in 1:(n_obs-1)){
    # pred_vars is for prediction variance (Matrices W in the book)
    pred_vars[,, t+1] <- Omega_h + B_h %*% filt_vars[,, t] %*% t(B_h) 
    # Ensure strict symmetry
    pred_vars[,, t+1] <- (pred_vars[,, t+1] + t(pred_vars[,, t+1])) / 2
    var_Ytp1 <- Sigma_h + A_h %*% pred_vars[,, t+1] %*% t(A_h)
    mean_Ytp1 <- as.numeric(A_h %*% B_h %*% filt_means[t, ])
    Kalman_gain <- pred_vars[,, t+1] %*% t(solve(var_Ytp1, A_h))
    filt_means[t+1, ] <- B_h %*% filt_means[t, ] + 
      Kalman_gain %*% (Y[t+1, ] - mean_Ytp1)
    # filt_vars[,, t+1] <- (Iz - Kalman_gain %*% A_h) %*% pred_vars[,, t+1]
    I_KA <- Iz - Kalman_gain %*% A_h
    filt_vars[,, t+1] <-
      I_KA %*% pred_vars[,, t+1] %*% t(I_KA) +
      Kalman_gain %*% Sigma_h %*% t(Kalman_gain)
    filt_vars[,, t+1] <- # Ensure strict symmetry
      (filt_vars[,, t+1] + t(filt_vars[,, t+1])) / 2
    all_ll[t+1] <- get_logdnorm(Y[t+1, ], 
                            mean_Ytp1,
                            var_Ytp1)
  }
  loglik <- sum(all_ll)
  list(filt_means = filt_means,
       filt_vars = filt_vars,
       pred_vars = pred_vars,
       loglik = loglik)
}



get_smoother_quants <- function(pars_h, filt_pars){
  # Model infos
  dz <- ncol(pars_h$B)
  n_obs <- nrow(filt_pars$filt_means)
  Iz <- diag(1, dz)
  # Creation of output objects
  smooth_means <- matrix(nrow = n_obs, 
                         ncol = dz)
  smooth_vars <- array(dim = c(dz, dz, n_obs))
  smooth_Ss <- array(dim = c(dz, dz, n_obs))
  smooth_covs <- array(dim = c(dz, dz, n_obs - 1)) 
  smooth_Cs <- array(dim = c(dz, dz, n_obs - 1))
  # Parameter extraction
  filt_means <- filt_pars$filt_means
  filt_vars <- filt_pars$filt_vars
  pred_vars <- filt_pars$pred_vars # W_t
  A_h <- pars_h$A; B_h <- pars_h$B; 
  Sigma_h <- pars_h$Sigma; Omega_h <- pars_h$Omega
  V0 <- pars_h$V0
  mu0 <- pars_h$mu0
  # Starting the smoother
  ## Initialization with the filter
  smooth_means[n_obs, ] <- filt_means[n_obs, ]
  smooth_vars[,, n_obs] <- filt_vars[,, n_obs]
  smooth_Ss[,, n_obs] <- smooth_vars[,, n_obs] + 
    smooth_means[n_obs, ] %*% t(smooth_means[n_obs, ])
  smooth_Ss[,, n_obs] <- (smooth_Ss[,, n_obs] + t(smooth_Ss[,, n_obs])) / 2
  for(t in (n_obs - 1):1){
    Ktilde <- filt_vars[,, t] %*% t(solve(pred_vars[,, t + 1], B_h))
    smooth_means[t, ] <- filt_means[t, ] + 
      Ktilde %*% (smooth_means[t + 1, ] - B_h %*% filt_means[t, ])
    smooth_vars[,, t] <- filt_vars[,, t] + 
      Ktilde %*% (smooth_vars[,, t+1] - pred_vars[,, t+1]) %*% t(Ktilde)
    # Ensure symmetry
    smooth_vars[,, t] <- (smooth_vars[,, t] + t(smooth_vars[,, t])) / 2
    smooth_covs[,, t] <- Ktilde %*% smooth_vars[,, t+1]
    smooth_Ss[,, t] <- smooth_vars[,, t] + 
      smooth_means[t, ] %*% t(smooth_means[t, ])
    smooth_Ss[,, t] <- (smooth_Ss[,, t] + t(smooth_Ss[,, t])) / 2
    smooth_Cs[,, t] <- smooth_covs[,, t] + 
      smooth_means[t, ] %*% t(smooth_means[t + 1, ])
  }
  list(smooth_means = smooth_means,
       smooth_vars = smooth_vars,
       smooth_covs = smooth_covs,
       smooth_Cs = smooth_Cs,
       smooth_Ss = smooth_Ss)
}

get_new_pars <- function(smooth_pars, Y, tYY, A = NULL,
                         B = NULL){
  n_obs <- nrow(Y)
  sum_Y_mu <- crossprod(Y, smooth_pars$smooth_means)
  sum_mu_Y <- t(sum_Y_mu)
  sum_S <- apply(smooth_pars$smooth_Ss,
                 c(1, 2), sum)
  Sum_S_m_n <- sum_S - smooth_pars$smooth_Ss[,,n_obs]
  Sum_S_m_1 <- sum_S - smooth_pars$smooth_Ss[,,1]
  sum_C <- apply(smooth_pars$smooth_Cs,
                 c(1, 2), sum)
  if(is.null(A)){
    A <- solve(sum_S, sum_Y_mu)
  }
  Sigma <- (tYY + A %*% sum_S %*% t(A) 
            - sum_Y_mu %*% t(A) - A %*% sum_mu_Y) / n_obs
  if(is.null(B)){
    B <- t(solve(Sum_S_m_n, sum_C))
  }
  mu0 <- smooth_pars$smooth_means[1, ]
  V0 <- smooth_pars$smooth_vars[,,1]
  Omega <- (B %*%  Sum_S_m_n %*%  t(B) + Sum_S_m_1 
            - t(sum_C) %*% t(B) - B %*% sum_C) / (n_obs - 1)
  Omega <- (Omega + t(Omega)) / 2 # Ensure symetry
  list(A = A, B = B, Sigma = Sigma, Omega = Omega,
       mu0 = mu0, V0 = V0)
}

# Code get_Q

get_Q <- function(pars, smooth_pars, Y) {
  
  n_obs <- nrow(Y)
  T_max <- n_obs - 1
  
  A <- pars$A
  B <- pars$B
  
  Sigma <- pars$Sigma
  Omega <- pars$Omega
  
  mu0 <- pars$mu0
  V0  <- pars$V0
  
  m <- smooth_pars$smooth_means
  S <- smooth_pars$smooth_Ss
  C <- smooth_pars$smooth_Cs
  
  dz <- nrow(V0)
  dy <- nrow(Sigma)
  
  ## ---- Initial state ----------------------------------------
  
  m0 <- m[1, ]
  S0 <- S[, , 1]
  
  E0 <- S0 -
    tcrossprod(m0, mu0) -
    tcrossprod(mu0, m0) +
    tcrossprod(mu0)
  
  Q0 <- -0.5 * determinant(V0, logarithm = TRUE)$modulus -
    0.5 * sum(diag(solve(V0, E0)))
  
  
  ## ---- State transitions ------------------------------------
  
  Qtrans <- -0.5 * T_max *
    determinant(Omega, logarithm = TRUE)$modulus
  
  for (t in 1:T_max) {
    St   <- S[, , t]
    Stp1 <- S[, , t + 1]
    Ct   <- C[, , t]
    Et <-
      Stp1 -
      B %*% Ct -
      t(Ct) %*% t(B) +
      B %*% St %*% t(B)
    Qtrans <- Qtrans -
      0.5 * sum(diag(solve(Omega, Et)))
  }
  ## ---- Observations -----------------------------------------
  Qobs <- -0.5 * n_obs *
    determinant(Sigma, logarithm = TRUE)$modulus
  for (t in 1:n_obs) {
    yt <- Y[t, ]
    mt <- m[t, ]
    St <- S[, , t]
    Et <-
      tcrossprod(yt) -
      tcrossprod(yt, mt) %*% t(A) -
      A %*% tcrossprod(mt, yt) +
      A %*% St %*% t(A)
    
    Qobs <- Qobs -
      0.5 * sum(diag(solve(Sigma, Et)))
  }
  Q0 + Qtrans + Qobs
}


check_matrix <- function(M, name, tol = 1e-10) {
  
  sym_error <- max(abs(M - t(M)))
  
  eig <- eigen(
    (M + t(M)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  cat(
    sprintf(
      "%-20s symmetry error = %.3e, min eigenvalue = %.3e, kappa = %.3e\n",
      name,
      sym_error,
      min(eig),
      kappa(M)
    )
  )
  
  invisible(
    list(
      symmetry_error = sym_error,
      min_eigenvalue = min(eig),
      kappa = kappa(M)
    )
  )
}
