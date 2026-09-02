# EM written in the cas of a single random effect.
exportFig <- FALSE

library(tidyverse) # For the %>% and plots
library(xtable) # To output the tables
donnees_NH4 <- read.table("../../Data/03_LinearMixedModel/donnees_NH4.txt",
                          sep = ";", header = TRUE)


quadrat_levels <- donnees_NH4 %>% 
  group_by(Sol, Quadrat) %>% 
  summarise(mN = mean(logNH4)) %>% 
  ungroup(Quadrat) %>% 
  mutate(m = mean(mN)) %>% 
  ungroup() %>% 
  arrange(m, mN) %>% 
  pull(Quadrat)

NH4_data <- donnees_NH4 %>% 
  rename(Soil = Sol,
         Plot = Quadrat) %>%
  mutate(Plot = factor(Plot, 
                       levels = quadrat_levels,
                       labels = paste0("Plot-", 1:9))) %>% 
  arrange(Plot, Soil) 


library(tidyverse)

get_log_lik <- function(pars, y, x, njs){
  n <- nrow(x)
  J <- length(njs)
  log_det_Sigma <- (n - J) * log(pars$sigma2) +
    sum(log(pars$sigma2 + njs * pars$gamma2))
  cst <- -0.5 * n * log(2 * pi)
  blocks <- lapply(njs, function(nj){
    diag(1 / pars$sigma2, nj) -
      matrix(pars$gamma2 / (pars$sigma2 * (pars$sigma2 + nj * pars$gamma2)), 
             nrow = nj, ncol = nj)
  })
  # Here, the correctness relies on the sorting of ys
  Sigma_inv <- do.call(what = Matrix::bdiag, args = blocks) %>% 
    as.matrix()
  quadr_term <- sum((y - x %*% pars$beta) * (Sigma_inv %*% (y - x %*% pars$beta)))
  cst - 0.5 * log_det_Sigma - 0.5 * quadr_term
}


get_estimations <- function(pars0, n_iter, Y, X, U){
  n <- nrow(X)
  output_pars <- matrix(nrow = n_iter + 1,
                        ncol = length(pars0$beta) + 2,
                        dimnames = list(NULL,
                                        c(paste0("beta[",1:length(pars0$beta),"]"),
                                          "gamma^2",
                                          "sigma^2")))
  output_pars[1, ] <- unlist(pars0)
  UtU <- t(U) %*% U
  XtXm1 <- solve(t(X) %*% X)
  ll <- rep(NA, n_iter + 1)
  ll[1] <- get_log_lik(pars0, Y, X, colSums(U))
  for(i in 1:n_iter){
    # E step
    Omega <- solve(diag(1 / pars0$gamma2, ncol(U)) + 
                     UtU / pars0$sigma2)
    mu <- 1 / pars0$sigma2 * Omega %*% t(U) %*% (Y - X %*% pars0$beta)
    # M step
    pars0$beta <- as.numeric(XtXm1 %*% t(X) %*% (Y - U %*% mu))
    pars0$sigma2 <- mean((Y - X %*% pars0$beta)^2) -
      2 * mean((U %*% mu) * (Y - X %*% pars0$beta)) +
      sum(diag(UtU %*% (Omega + mu %*% t(mu)))) / n
    pars0$gamma2 <- sum(diag(Omega + mu %*% t(mu))) / ncol(U)
    ll[i + 1] <- get_log_lik(pars0, Y, X, colSums(U))
    output_pars[i + 1, ] <- unlist(pars0)
  }
  # Computing Fisher information estimate
  YmXB <- Y - X%*%pars0$beta
  var_score_beta <- t(X) %*% U %*% Omega %*% t(U) %*% X / (pars0$sigma2)^2
  Esp_JS <- -t(X) %*% X / pars0$sigma2
  V <- U %*% diag(pars0$gamma2, ncol(U)) %*% t(U) + 
    pars0$sigma2 * diag(1, n)
  # Formating results
  pars_df <- output_pars %>% 
    as.data.frame() %>% 
    mutate(Iteration = 0:n_iter) %>% 
    pivot_longer(cols = -c("Iteration"),
                 values_to = "Estimate",
                 names_to = "Parameter")
  list(ll = ll,
       best_par = pars0,
       mu_post = mu,
       Omega_post = Omega,
       pars_df = pars_df,
       Esp_JS = Esp_JS,
       var_score_beta = var_score_beta,
       cov_beta_hat =  -solve(var_score_beta + Esp_JS))
}


# Exemple d'utilisation et comparaison avec lmer --------------------------

library(lme4)

set.seed(123)
# dummy_x <- rnorm(nrow(NH4_data)) # Pour voir si ça marche en général
X <- model.matrix(lm(logNH4 ~ Soil, data = NH4_data)) 
# %>% cbind(x = dummy_x)
U <- model.matrix(lm(logNH4 ~ Plot - 1, data = NH4_data))
Y <- NH4_data$logNH4

subsample <- 1:nrow(NH4_data)


my_result <- get_estimations(pars0 = list(beta = rep(0, ncol(X)),
                                          sigma2 = 0.5,
                                          gamma2 = 0.1), 
                             n_iter =  500, 
                             Y = Y[subsample], 
                             X = X[subsample, ], 
                             U = U[subsample,])


# Comparison with lmer ----------------------------------------------------

result_lmer <- lmer(logNH4 ~ Soil + 
                      # x + 
                      (1|Plot), 
                    data = NH4_data %>% 
                      # mutate(x = dummy_x) %>% 
                      slice(subsample),
                    REML = FALSE)
## Fixed effects
beta_hat <- result_lmer@beta
beta_hat
my_result$best_par$beta

## Variances
print(VarCorr(result_lmer),comp="Variance")
c(my_result$best_par$gamma2, my_result$best_par$sigma2)

## Covariance of beta
cov_beta_hat <- my_result$cov_beta_hat
cov_beta_hat
vcov(result_lmer)

# Tables

## Variance beta_hat
printed_cov <- xtable(my_result$cov_beta_hat, digits = 3, align = "cccc")
print.xtable(printed_cov, tabular.environment="pmatrix", include.rownames = FALSE, include.colnames = FALSE)

## Confidence intervals

contrat_mat <- matrix(c(1, 0, 0,
                        1, 1, 0,
                        1, 0, 1), nrow = 3, byrow = TRUE)
alpha_lev <- 0.01
mean_variances <- matrixmy_result$cov_beta_hat
tibble(Soil = sort(unique(NH4_data$Soil)),
       mean_resp = as.numeric(contrat_mat %*% my_result$best_par$beta),
       IC_low = mean_resp - qnorm(1 - alpha_lev / 2) * sqrt(diag(contrat_mat %*% 
                                                                   cov_beta_hat %*%
                                                                   t(contrat_mat))) ,
       IC_sup = mean_resp + qnorm(1 - alpha_lev / 2) * sqrt(diag(contrat_mat %*% 
                                                                   cov_beta_hat %*%
                                                                   t(contrat_mat)))) %>% 
  rename("Expected log$_{10}$NH$_4$" = mean_resp) %>% 
  mutate("99% confidence interval" = paste0("[", round(IC_low, 3), 
                                            ", ", round(IC_sup, 3), "]"),
         .keep = "unused") %>% 
  xtable(digits = 3) %>% 
  print(include.rownames = FALSE, escape = FALSE)

# Figures -----------------------------------------------------------------

graph_width <- 5
graph_height <- 5
graph_unit <- "in"
library(furrr)

future::plan(strategy = "multisession", workers = 10)

set.seed(1234)
all_results <- future_map(1:5,
                          function(start_id){
                            start <- list(beta = runif(3, -1, 1),
                                          sigma2 = runif(1, 0.01, 1),
                                          gamma2 = runif(1, 0.01, 0.2))
                            result <- get_estimations(start, 
                                                      n_iter =  200, 
                                                      Y = Y, X = X, U = U)
                            result$pars_df$Replicate = start_id
                            return(result)
                          }, .options=furrr_options(seed = TRUE))

all_pars <- map_dfr(all_results, "pars_df")
beta_traj_plot <- ggplot(all_pars %>% 
                           filter(str_detect(Parameter, "beta"))) + 
  aes(x = Iteration, y = Estimate, group = Replicate, 
      color = factor(Replicate)) + 
  facet_wrap(~Parameter, labeller = label_parsed, nrow = 1, scales = "free") + 
  geom_path() +
  theme_bw() +
  scale_color_viridis_d() +
  theme(legend.position = "none")
if(exportFig)
  ggsave("../../Figures/borneo-beta-em-trajectory.png",
         plot = beta_traj_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)  


means <- all_results[[1]]$mu_post
vars <- diag(all_results[[1]]$Omega_post)

post_density_plot <- map_dfr(1:9, function(i){
  mu <- means[i]
  sd <- sqrt(vars[i])
  xs <- seq(-.4, .4, length.out = 301)
  data.frame(x = xs,
             dens = dnorm(xs, mu, sd),
             Plot = paste0("Plot-", i))
}) %>% 
  ggplot(aes(x = x, y = dens)) + 
  geom_line(aes(color = Plot)) +
  geom_ribbon(aes(xmin = x, xmax = x, ymin = 0, ymax = dens,
                  fill = Plot), alpha = .3) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.position = "none") +
  labs(x = "Z|Y", y = "Posterior density")
if(exportFig)
  ggsave("../../Figures/borneo-posterior-density.png",
         plot = post_density_plot,
         width = graph_width,
         height = graph_height,
         units = graph_unit)  

