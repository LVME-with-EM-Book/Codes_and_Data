EM_2Mixture= function(params0,y,tol=10^{-5},estim_only_mu=FALSE){
  
  #-------------------------- 
  params_mat <- matrix(unlist(params0),ncol=5,nrow = 1)
  ll_vec =c(logLikelihood(params0,y))
  #-------------------------
  delta = 10
  iter = 0 
  params =  params0 
  n <- length(y)
  
  while (delta > tol){
    print(iter)
    iter = iter + 1 
    #------------- E step 
    rho1 = params$p*dnorm(y,params$mu[1],sqrt(params$var[1]))
    rho2 = (1-params$p)*dnorm(y,params$mu[2],sqrt(params$var[2]))
    tau = rho1/(rho1 + rho2)
    if (sum(tau)==0) {tau=rep(0.001,n)}
    #----------------- M step
    N1 = sum(tau)
    N2 =  n-N1
    mu_new <- c(0,0)
    mu_new[1] <- sum(tau*y)/N1
    mu_new[2] <- sum((1-tau)*y)/N2
   
    if(!estim_only_mu){
      var_new <- c(0,0)
      var_new[1] <- 1/N1*sum(tau*(y-mu_new[1])^2)
      var_new[2] <- 1/N2*sum((1-tau)*(y-mu_new[2])^2)
      p_new <- mean(tau) 
    }
    #---------- delta
    delta = sum((params$mu-mu_new)^2)
    
    #----------------- LL evaluation
    params$mu <- mu_new
    if(!estim_only_mu){
      params$var <- var_new
      params$p <- p_new
    }
   
    ll_vec = c(ll_vec,logLikelihood(params,y))
    params_mat <- rbind(params_mat,unlist(params))
    
  }
  
  params_mat <- as.data.frame(params_mat)
  params_mat[,6] = ll_vec
  names(params_mat) = c('mu1','mu2','var1','var2','p','loglik')
  rownames(params_mat) = 1:nrow(params_mat)
  
  
  return(params_mat)
  
}


logLikelihood = function(params,y){
  dy1 <- params$p*dnorm(y,params$mu[1],sqrt(params$var[1])) 
  dy2 <- (1-params$p)*dnorm(y,params$mu[2],sqrt(params$var[2]))
  dy <- dy1 + dy2
  return(sum(log(dy)))
}
