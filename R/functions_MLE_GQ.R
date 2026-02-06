##########################################################
### functions for calculating MLE of parameters of trivariate (q=3) ADMM model
### using Gaussian quadrature approximation of log-likelihood 

### Autor: Tomas
### Version: August 2025


###################  function calculating the integrand exp(hd.ud) 

Integrand <- function(u.k, beta, sigma, yy.d, X.dk, phi){
  
  sigma.ud <- u.k * sigma
  
  eta.1 <- exp(X.dk[[1]] %*% beta[[1]] + sigma.ud[1])
  eta.2 <- exp(X.dk[[2]] %*% beta[[2]] + sigma.ud[2])
  
  sum.eta <- eta.1 + eta.2
  
  mu <- rep(0,3)
  
  mu[1] <- eta.1/(1+sum.eta)
  mu[2] <- eta.2/(1+sum.eta)
  mu[3] <- 1/(1+sum.eta)
  
  alpha <- phi * mu
  #sigma.v <- as.vector(Reduce(sigma, f=cbind))
  
  log.y <- log(yy.d)
  log.g <- log(gamma(alpha))
  
  hd <- (alpha-1) %*% log.y - sum(log.g) - sum(u.k^2)/2
  
  result <- exp(hd)    
  
  return (result)
  
}
  
### function evaluating the function Integrand for a matrix of values
### neede for the quadrature function

# Integrand.GQ <- function(u.k, beta, sigma, yy.d, X.dk, phi){
#   
#   SS <- dim(u.k)[1]
#   result <- rep(0,SS)
#   
#   for (ss in 1:SS){
#     
#     u.k.ss <- u.k[ss,] 
#     result[ss] <- Integrand(u.k.ss, beta, sigma, yy.d, X.dk, phi) 
#     
#   }
#   
#   return (result)
#   
# }
  
### another version without for cycle
Integrand.GQ <- function(u.k, beta, sigma, yy.d, X.dk, phi) {
  result <- apply(u.k, 1, function(u.k.ss) {
    Integrand(u.k.ss, beta, sigma, yy.d, X.dk, phi)
  })
  return(result)
}

  
# # H-cubature approximation of the log-likelihood  
# 
# log.lik.H <- function(param, XX.k, yy,  phi.d){
#   
#   beta.l <- list() 
#   beta.l[[1]] <- param[1:p1]
#   beta.l[[2]] <- param[(p1+1):p]
#   sigma.l <- param[(p+1):(p+2)]
#   
#   l.d <- rep(0,D)
#   
#   for (d in 1:D) {
#     
#     phi <- phi.d[d]
#     X.dk <- lapply(XX.k,function(x) x[d,])   # taking out lines d from all the elements of the list
#     # other option lapply(XX.k,'[',d,)
#     yy.d <- yy[d,]
#     
#     # calculating the inner integral 
#     pom <- hcubature(Integrand, rep(-Inf,2), upperLimit=rep(Inf,2), beta.l, sigma.l, yy.d, X.dk, phi)
#     int.val <- pom$integral
#     
#     # calculating the summand
#     l.d[d] <-  log(int.val)     # I did not include the constants log(Gamma(phi_d)) and -log(2pi) 
#     
#     }  # cycle d
#   
#   log.lik <- - sum(l.d)   # we put minus since we use minimization of log.lik in the sequel
#   
#   return(log.lik)     
# }


# Gaussian quadrature approximation of the log-likelihood  

log.lik.GQ <- function(param, XX.k, yy,  phi.d, nw){
  
  beta.l <- list() 
  beta.l[[1]] <- param[1:p1]
  beta.l[[2]] <- param[(p1+1):p]
  sigma.l <- param[(p+1):(p+2)]
  
  l.d <- rep(0,D)
  
  for (d in 1:D) {
    
    phi <- phi.d[d]
    X.dk <- lapply(XX.k,function(x) x[d,])   # taking out lines d from all the elements of the list
    # other option lapply(XX.k,'[',d,)
    yy.d <- yy[d,]
    
    # calculating the inner integral 
    
    int.val <- quadrature(Integrand.GQ, grid = nw, beta.l, sigma.l, yy.d, X.dk, phi)
    
    # calculating the summand
    l.d[d] <-  log(int.val)     # I did not include the constants log(Gamma(phi_d)) and -log(2pi) 
    
  }  # cycle d
  
  log.lik <- - sum(l.d)   # we put minus since we use minimization of log.lik in the sequel
  
  return(log.lik)     
}



###############################################################
# MLE Gaussian quadrature algorithm 

MLE_GQ <- function(XX.k, yy, phi.d, beta.s, sigma.s, nw){
  
  # beta  <- beta.s
  # sigma <- sigma.s
  param.0 <- c(beta.s[[1]], beta.s[[2]], sigma.s)   # seeds for beta, sigma
  
  
  # Turn warnings into errors so they can be trapped
  options(warn = 2)
  
  # optimization of the Gaussian quadrature approximation of the log-likelihood
  res.opt <- try( nloptr(param.0, eval_f = log.lik.GQ, lb = c(rep(-Inf, p), rep(0, q-1)), ub = c(rep(Inf, p+q-1)), opts = opts,
                         XX.k=XX.k, yy=yy, phi.d=phi.d, nw=nw), silent = TRUE)
  
  
  convergence <- TRUE
  
  if (inherits(res.opt, "try-error")) {     # checking if there was no error in nloptr
    convergence <- FALSE
  }
  
  if (res.opt$status != 1){             # checking convergence in nloptr
    convergence <- FALSE
  }
  
  #print(res.opt)
  param.1 <- res.opt$solution
  
  
  return(list(par.est = param.1, convergence = convergence))
  
}

