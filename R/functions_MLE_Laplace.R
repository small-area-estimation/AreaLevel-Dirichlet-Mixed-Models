##########################################################
### functions for calculating MLE of parameters of trivariate (q=3) ADMM model
### using Laplace approximation of log-likelihood 

### Autor: Tomas
### Version: September 2024


### function calculating the score vector and the Hessian matrix
### (defined on page 13 of the actual version of the paper)

ScoreHes <- function(u.k, beta, sigma, yy.d, X.dk, phi){
  
  Sd <- matrix(0,2,1)
  Hd <- matrix(0,2,2)

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
  
  
  # calculating first derivatives of mu_{dj}  (formula on page 11)
  d.mu <- rep(list(c(0,0)), 3)
  
  # first derivatives of mu_d1
  d.mu[[1]][1] <- sigma[1]*mu[1]*(1-mu[1])
  d.mu[[1]][2] <- -sigma[2]*mu[1]*mu[2]
  
  # first derivatives of mu_d2
  d.mu[[2]][1] <- -sigma[1]*mu[1]*mu[2]
  d.mu[[2]][2] <- sigma[2]*mu[2]*(1-mu[2])
  
  # first derivatives of mu_d3
  d.mu[[3]][1] <- -sigma[1]*mu[1]*mu[3]
  d.mu[[3]][2] <- -sigma[2]*mu[2]*mu[3]
  
  
  # calculating second derivatives of mu_{dj}  (formulas on page 12)
  dd.mu <- lapply(1:3, matrix, data= NA, nrow=2, ncol=2)
  
  # second derivatives of mu_d1
  dd.mu[[1]][1,1] <- sigma[1]^2*mu[1]*(1-mu[1])*(1-2*mu[1])
  dd.mu[[1]][2,2] <- -sigma[2]^2*mu[1]*mu[2]*(1-2*mu[2])
  dd.mu[[1]][1,2] <- dd.mu[[1]][2,1] <- -sigma[1]*sigma[2]*mu[1]*mu[2]*(1-2*mu[1])
  
  # second derivatives of mu_d2
  dd.mu[[2]][1,1] <- -sigma[1]^2*mu[1]*mu[2]*(1-2*mu[1])
  dd.mu[[2]][2,2] <- sigma[2]^2*mu[2]*(1-mu[2])*(1-2*mu[2])
  dd.mu[[2]][1,2] <- dd.mu[[2]][2,1] <- -sigma[1]*sigma[2]*mu[1]*mu[2]*(1-2*mu[2])
  
  # second derivatives of mu_d3
  dd.mu[[3]][1,1] <- -sigma[1]^2*mu[1]*mu[3]*(1-2*mu[1])
  dd.mu[[3]][2,2] <- -sigma[2]^2*mu[2]*mu[3]*(1-2*mu[2])
  dd.mu[[3]][1,2] <- dd.mu[[3]][2,1] <- 2*sigma[1]*sigma[2]*mu[1]*mu[2]*mu[3]
  
  
  log.y <- log(yy.d)
  psi <- digamma(alpha)
  psi.1 <- trigamma(alpha)
  
  
  for (k in 1:2){
    d.mu.k <- unlist(lapply(d.mu,'[',k))   # taking out elements k from all the elements of the list
    dd.mu.k <- unlist(lapply(dd.mu,'[', k, k))
    
    Sd[k,1] <- phi *( d.mu.k %*% log.y - d.mu.k %*% psi) - u.k[k]
    Hd[k,k] <- phi *(dd.mu.k %*% log.y - dd.mu.k %*% psi) - phi^2*(psi.1 %*% d.mu.k^2) - 1
    
  }
  
  dd.mu.12 <- unlist(lapply(dd.mu,'[', 1, 2))
  d.mu.1 <- unlist(lapply(d.mu,'[',1))
  d.mu.2 <- unlist(lapply(d.mu,'[',2))
  Hd[1,2] <- Hd[2,1] <- phi *(dd.mu.12 %*% log.y - dd.mu.12 %*% psi) - phi^2*(psi.1 %*% (d.mu.1*d.mu.2))
  
  
  #  calculating the function hd(ud)    (formula (5.5))
  log.g <- log(gamma(alpha))
  hd <- (alpha-1) %*% log.y - sum(log.g) - sum(u.k^2)/2
  
  result <- list(Sd=Sd, Hd=Hd, hd=hd)
  
  return (result)

}
  


# Newton-Raphson algorithm for random effects u  (formula (5.6))

NewtonRaphson_u <- function(ud.seed, beta, sigma, yy.d, X.dk, phi, epsilon, maxiter)
{
  
 udk.0 <- udk.1 <- ud.seed    
  
  num.iter <- 0
  convergence <- FALSE
  
  while (convergence == FALSE & num.iter<=maxiter)
  {
    # Turn warnings into errors so they can be trapped
    options(warn = 2)
    
    SH <- try(ScoreHes(udk.0, beta, sigma, yy.d, X.dk, phi), silent = TRUE)
    
    if (inherits(SH, "try-error")) {     # checking if there was no error in calculating the matrix Hd
      break
    }
    
    Sd <- SH$Sd
    Hd <- SH$Hd
    
    udk.1 <- try(udk.0 - solve(Hd) %*% Sd, silent = TRUE) 
    
    if (inherits(udk.1, "try-error")) {     # checking if there was no error in inverting the matrix Hd
      break
    }
    else{
    
    diff <- abs(udk.1-udk.0)
    
    if (any(diff>=epsilon)){
      udk.0 <- udk.1
      num.iter  <- num.iter+1
    }
    else{
      convergence <- TRUE
      num.iter  <- num.iter+1
    }
    
    }

    #cat("Iter for u ", num.iter,"\n")

  } # while
  
  result <- list(u=udk.1, num.iter, convergence)
  
  return (result)
}     



# Laplace approximation of the log-likelihood  (formula (5.7))

log.lik.L <- function(param, XX.k, yy,  phi.d, udk.NR){
  
  beta.l <- list() 
  beta.l[[1]] <- param[1:p1]
  beta.l[[2]] <- param[(p1+1):p]
  sigma.l <- param[(p+1):(p+2)]
  
  l.0d <- rep(0,D)
  
  for (d in 1:D) {
    
    phi <- phi.d[d]
    X.dk <- lapply(XX.k,function(x) x[d,])   # taking out lines d from all the lements of the list
    # other option lapply(XX.k,'[',d,)
    yy.d <- yy[d,]
    udk.0 <- udk.NR[[d]]
    
    SH <- ScoreHes(udk.0, beta.l, sigma.l, yy.d, X.dk, phi)
    
    hd <- SH$hd
    Hd <- SH$Hd
    
    l.0d[d] <-  hd - log(det(-Hd))/2     # I did not include the constant log Gamma(phi_d) and I use the function h_d (5.5) 
    
  }  # cycle d
  
  log.lik <- - sum(l.0d)   # we put minus since we use minimization of log.lik in the sequel
  
  return(log.lik)     
}


###############################################################
# MLE-Laplace algorithm  (page 14)

MLE_Laplace <- function(XX.k, yy, phi.d, beta.s, sigma.s, seeds.u, maxiter=1000, epsilon.u=0.0001, epsilon.par=0.00001){

  beta  <- beta.s
  sigma <- sigma.s
  param.0 <- param.1 <- c(beta[[1]], beta[[2]], sigma)   # seeds for beta, sigma
  
  udk.NR.0 <- udk.NR <- seeds.u       # seeds for u_dk
  
  
  num.iter <- 0
  convergence <- FALSE
  
  
  while (convergence == FALSE & num.iter<=maxiter)
  {
    
    for (d in 1:D) {
      
      phi <- phi.d[d]
      X.dk <- lapply(XX.k,function(x) x[d,])   # taking out lines d from all the lements of the list
      # other option lapply(XX.k,'[',d,)
      yy.d <- yy[d,]
      
      ud.seed <- udk.NR.0[[d]]
      
      res.u <- NewtonRaphson_u(ud.seed, beta, sigma, yy.d, X.dk, phi, epsilon.u, maxiter)
      
      if (res.u[[3]] == FALSE){
        return(list(par.est = param.1, udk.pred = udk.NR, num.iter = num.iter, convergence = convergence))   # to terminate the double loop
      }
      
      udk.NR[[d]] <- as.vector(res.u[[1]])
      
    }
    
    
    # Turn warnings into errors so they can be trapped
    options(warn = 2)
    
    # optimization of the Laplace approximation of the log-likelihood
    res.opt <- try( nloptr(param.0, eval_f = log.lik.L, lb = c(rep(-Inf, p), rep(0, q-1)), ub = c(rep(Inf, p+q-1)), opts = opts,
                      XX.k=XX.k, yy=yy, phi.d=phi.d, udk.NR=udk.NR), silent = TRUE)
    
    if (inherits(res.opt, "try-error")) {     # checking if there was no error in nloptr
      break
    }
    
    if (res.opt$status != 1){             # checking convergence in nloptr
      break
    }
    
    #print(res.opt)
    param.1 <- res.opt$solution
    
    #### here we check the diference between iterations
    
    u.NR.0 <- Reduce(udk.NR.0, f=rbind)
    u.NR <- Reduce(udk.NR, f=rbind)
    
    diff.u <- abs(u.NR - u.NR.0)
    diff.par <- abs(param.1 - param.0)
    
    if (any(diff.u>=epsilon.u) | any(diff.par>=epsilon.par)){
      udk.NR.0 <- udk.NR
      param.0 <- param.1
      num.iter  <- num.iter+1
    }
    else{
      convergence <- TRUE
      num.iter  <- num.iter+1
    }
    
    cat("Iter Laplace = ", num.iter,"\n")
    
  } # while
  
  return(list(par.est = param.1, udk.pred = udk.NR, num.iter = num.iter, convergence = convergence))
  
}








###################  function hd.ud for checking the NR algorithm for u.dk

hd.ud <- function(u.k, beta, sigma, yy.d, X.dk, phi){
  
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
  
  result <- -hd    #  for maximization we put -
  
  return (result)
  
}

