#############################################################
# EM for a two-component univariate Gaussian mixture
#############################################################
# The mixture model writes
#     y_i ~ p * N(mu1, var1) + (1 - p) * N(mu2, var2),
# where p is the mixing weight of component 1 and the unknown
# parameters are theta = (mu1, mu2, var1, var2, p).
# These two functions (logLikelihood, EM_2Mixture) perform the
# maximum-likelihood fit by the Expectation-Maximisation algorithm.
#############################################################

#############################################################
# Observed-data log-likelihood of the mixture
#############################################################
# For a candidate parameter set params (a list with components
#   - mu  : vector c(mu1, mu2)
#   - var : vector c(var1, var2)
#   - p   : mixing weight of component 1),
# returns the log-likelihood of the observed data y:
#
#   log L(theta; y) = sum_i log( p*phi(y_i;mu1,var1)
#                                + (1-p)*phi(y_i;mu2,var2) )
#
# i.e. the sum over observations of the log of the marginal
# density (with the latent memberships Z_i integrated out).
logLikelihood = function(params, y){
  # Unnormalised density contribution of component 1 (weight p)
  dy1 <- params$p*dnorm(y, params$mu[1], sqrt(params$var[1]))
  # Unnormalised density contribution of component 2 (weight 1-p)
  dy2 <- (1 - params$p)*dnorm(y, params$mu[2], sqrt(params$var[2]))
  # Marginal density of each observation
  dy <- dy1 + dy2
  return(sum(log(dy)))
}

#############################################################
# EM algorithm for the two-component Gaussian mixture
#############################################################
# Fits the mixture by alternating E and M steps until the L2
# distance between two successive means, delta, drops below tol.
#
# Arguments
#   params0         : initial parameters, a list(mu, var, p)
#   y               : vector of observed data
#   tol             : convergence tolerance on delta (default 1e-5)
#   estim_only_mu   : if TRUE, only the two means are updated
#                     (variances and p stay fixed at their
#                     initial values); if FALSE (default), all
#                     parameters are estimated.
#
# Returns a data.frame with one row per EM iteration containing
#   mu1, mu2, var1, var2, p, loglik
# where loglik is the observed-data log-likelihood.
EM_2Mixture = function(params0, y, tol = 10^{-5}, estim_only_mu = FALSE){

  #----------- initialisation
  # First row of the output: the initial parameters
  params_mat <- matrix(unlist(params0), ncol = 5, nrow = 1)
  # Initial log-likelihood
  ll_vec <- c(logLikelihood(params0, y))
  # Convergence criterion: start above the tolerance
  delta <- 10
  iter  <- 0
  params <- params0
  n <- length(y)

  #----------- EM iterations
  while (delta > tol){
    #print(iter)
    iter <- iter + 1

    #---------- E step : posterior membership weights
    # Unnormalised weight of component 1 for each observation
    rho1 <- params$p*dnorm(y, params$mu[1], sqrt(params$var[1]))
    # Unnormalised weight of component 2
    rho2 <- (1 - params$p)*dnorm(y, params$mu[2], sqrt(params$var[2]))
    # Normalised posterior probability that y_i comes from component 1
    tau <- rho1/(rho1 + rho2)
    # Guard against a degenerate fit (no weight allocated at all)
    if (sum(tau) == 0){ tau <- rep(0.001, n) }

    #---------- M step : responsibility-weighted updates
    # Effective sample size of each component
    N1 <- sum(tau)      # component 1
    N2 <- n - N1        # component 2
    # Update the means
    mu_new <- c(0, 0)
    mu_new[1] <- sum(tau*y)/N1
    mu_new[2] <- sum((1 - tau)*y)/N2

    # Update variances and mixing weight (only if requested)
    if (!estim_only_mu){
      var_new <- c(0, 0)
      var_new[1] <- 1/N1*sum(tau*(y - mu_new[1])^2)
      var_new[2] <- 1/N2*sum((1 - tau)*(y - mu_new[2])^2)
      p_new <- mean(tau)
    }

    #---------- convergence check (on the means only)
    delta <- sum((params$mu - mu_new)^2)

    #---------- commit the updates
    params$mu <- mu_new
    if (!estim_only_mu){
      params$var <- var_new
      params$p <- p_new
    }

    #---------- record the log-likelihood and the parameters
    ll_vec <- c(ll_vec, logLikelihood(params, y))
    params_mat <- rbind(params_mat, unlist(params))
  }

  #----------- format and return the trace of iterations
  params_mat <- as.data.frame(params_mat)
  params_mat[, 6] <- ll_vec
  names(params_mat) <- c('mu1', 'mu2', 'var1', 'var2', 'p', 'loglik')
  rownames(params_mat) <- 1:nrow(params_mat)

  return(params_mat)
}
