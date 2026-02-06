###############################################################################
###
###        Simulation experiment 1  Dirichlet regression ADMM           
###                                                            
###                            August 2025                         
###
### Author: Tomas Hobza


rm(list = ls())



source("functions_MLE_GQ.R")


library(lme4)
library(Matrix)
library(mvtnorm)
library(cubature)
library(DirichletReg)
library(nloptr)



opts <- list( "algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=1000)


# package for Gaussian qudrature 
library(mvQuad)

# create grid
nw <- createNIGrid(dim=2, type="GHe", level=25)

#nw1 <- createNIGrid(dim=2, type="GHN", level=25)
#plot(nw)



set.seed(17032005)

D <- 50    # number of areas, used values {50, 75, 100, 150, 200, 300}
q <- 3        # number of categories

beta.true <- rbind(c(-1/2,1),c(1,-1/3))    # vectors beta_k in rows
beta.k <- split(beta.true, row(beta.true))

p1 <- 2    #  number of parameters beta1, beta2
p2 <- 2
p <- p1+p2

sigma.true <- as.matrix(c(1/3, 1/2))            #  variance parameters sigma_1, sigma_2
sigma.k <- split(sigma.true, row(sigma.true))

ind.d <- seq(1,D)
phi.d <- exp(2 + 1/ind.d)                 # parameters phi_d  (assumed to be known in application)


### generating the regressors x_dk

mu.x <- c(1,1)
sigma.x <- c(1,3/2)

XX.k <- list()

for (k in 1:(q-1)){
  
  unif <- as.matrix(runif(D))
  x.dk1 <- as.matrix(rep(1,D))
  x.dk2 <- mu.x[k] + sqrt(sigma.x[k]) * unif
  
  XX.k[[k]] <- cbind(x.dk1, x.dk2)
  
}


### SIMULATION STARTS

ITER <- 1000   # number of iterations

param.est <- list() 

###############################################
# start of the simulation loop (iter)

iter <- 1
no.convergence <- 0

while (iter <= ITER) {
  
  # Simulation of random effects ud (normal)  
  udk <-  split(rmvnorm(D, rep(0, (q-1))), 1:(q-1))
  
  
  # calculating mu_dk (step 2.2 of the simulation)
  
  sigma.ud <- Map(udk, f="*", sigma.k)
  
  eta.k <- Map(Map(XX.k, f="%*%", beta.k), sigma.ud, f="+")
  
  mu.k <- list()
  
  eta <- exp(Reduce(eta.k, f=cbind))
  
  sum.eta <- rowSums(eta)
  
  mu.k[[q]] <- 1/(1+sum.eta)
  
  for (k in 1:(q-1)){
    mu.k[[k]] <- eta[,k]/(1+sum.eta)
  }
  
  mu <- Reduce(mu.k, f=cbind )    # matrix of mu_dk
  
  ## generating y_d
  
  alpha <- phi.d * mu
  
  yy <- rdirichlet(D, alpha)
  
  
  ### fitting of the model with MLE - Gaussian quadrature
  
  # seeds for parameters 
  beta.s <- beta.k
  sigma.s <-  unlist(sigma.k)
  
  #seeds.u <- split(matrix(0, D, 2), 1:D)   # seeds for u_dk
  
  res.MLE <- MLE_GQ(XX.k, yy, phi.d, beta.s, sigma.s, nw)
  
  param.est[[iter]] <- res.MLE$par.est
  
  if (res.MLE$convergence == FALSE){
    iter <- iter - 1
    no.convergence <- no.convergence + 1
  }
  
  cat("Iteration overall:  ", iter, "\n")
  
  iter <- iter + 1
  
  }  ### cycle iter


#### Calculating output (Step 3 of the simulation 1)

param.true <- c(beta.k[[1]], beta.k[[2]], sigma.true)

Bias <- (Reduce( lapply(param.est, param.true, FUN="-"), f="+"))/ITER

Rmse <- ((Reduce(lapply(lapply(param.est, param.true, FUN="-"), FUN="^",2), f="+"))/ITER)^0.5

RBias <- Bias/abs(param.true)*100
RRmse <- Rmse/abs(param.true)*100


#### Saving results 

results.param <- data.frame(Bias, Rmse, RBias, RRmse)

rownames(results.param) <- c("beta_11", "beta_12", "beta_21", "beta_22", "sigma1", "sigma2")

# write.table(results.param, file=paste0("Results_GQ/Sim1_D=", D, ".txt"), sep="\t")

#save(list=ls(), file=paste("Results_GQ/Sim1_D=", D, ".Rdata", sep=""))









