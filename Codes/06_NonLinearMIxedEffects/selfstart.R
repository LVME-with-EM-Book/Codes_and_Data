

# 1)  SELFSTART  FUNCTION FOR DECLINING 4-PARAMETERS LOGISTC RESPONSE  ######################################

# Response function:  lfmc = (A - w) / (1 + exp((m - time)/s))) + w 

# Response variable: "lfmc" 

# Parameters: .- "A": upper asymptote
#             .- "w": lower asymptote 
#             .- "m": time when lfmc is midway between w and A 
#             .- "s": controls the slope 

# Predictor: "time"


# 1.1 ---- Initial values ---------------------------------------

dlfInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a dlf")
  }
  z <- xy[["y"]]
  x1 <- xy[["x"]]
  A <- max(xy[,"y"])
  w <- min(xy[,"y"])
  m <- NLSstClosestX(xy, (A + w)/2)
  s <- -m/2 # this is a crude estimate for "s"
  # this is an alternative method of getting a good guess for "s"
  # it uses the idea that this decreasing logistic is really an inverted logistic
  
  zz <- z - w # subtract W to have zero as the lowest value
  
  ## "SSlogis" is used to find out a guess for S, but be silent about it
  s2 <- try(coef(nls(zz ~ SSlogis(-x1, A, m, s)))[3], silent = TRUE)
  ## If the previous line returns an error we are silent about it and maintain the crude guess for "s"
  
  if(class(s2) != "try-error")  s <- -s2; #  print(c(-s2, s))
  value <- c(A, w, m, s)
  names(value) <- mCall[c("A","w","m","s")]
  value
  
}


# 1.2 Main function and partial derivatives (gradient) --------------------------------

dlf <- function(time, A, w, m, s){
  
  .expr1 <- (m - time)/s
  .expr2 <- 1 + exp(.expr1)
  .value <- (A - w) / .expr2 + w
  # derivative with respect to A
  .expi1 <- 1/.expr2 
  ## Derivative with respect to W
  .expi2 <- 1 - 1/.expr2 
  ## Derivative with respect to M
  .expr3 <- A - w
  .expr4 <- exp((m - time)/s)
  .expr5 <- 1 + .expr4
  .expi3 <- -(.expr3 * (.expr4 * (1/s))/.expr5^2)
  ## Derivative with respect to S
  .eexpr1 <- A - w
  .eexpr2 <- m - time
  .eexpr4 <- exp(.eexpr2/s)
  .eexpr5 <- 1 + .eexpr4
  .expi4 <- .eexpr1 * (.eexpr4 * (.eexpr2/s^2))/.eexpr5^2
  
  .actualArgs <- as.list(match.call()[c("A", "w", "m", "s")])
  
  #  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("A", "w", "m","s")))
    .grad[, "A"] <- .expi1
    .grad[, "w"] <- .expi2
    .grad[, "m"] <- .expi3 
    .grad[, "s"] <- .expi4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}


# 1.3 ---- selfStart function ---------------------------------------- 

SSdlf <- selfStart(dlf, initial = dlfInit, c("A","w","m","s")) # selfStart function

