###############################################################################
###
###        Simulation experiment 3  Dirichlet regression ADMM           
###                                                            
###                            July 2024                         
###
### Author: Maria Dolores Esteban, Agustin Perez and Tomas Hobza




rm(list = ls())


source("functions_MLE_Laplace.R")
source("functions_predictors.R")
source("boot_Dirichlet.R")


library(lme4)
library(Matrix)
library(mvtnorm)
library(cubature)
library(DirichletReg)
library(nloptr)



opts <- list( "algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=10000)


set.seed(17032005)

D <- 100    # number of areas
q <- 3      # number of categories

beta.true <- rbind(c(-1/2,1),c(1,-1/3))    # vectors beta_k in rows
beta.k <- split(beta.true, row(beta.true))

p1 <- 2    #  number of parameters beta1, beta2
p2 <- 2
p <- p1+p2

sigma.true <- as.matrix(c(1/3, 1/2))    #  variance parameters sigma_1, sigma_2
sigma.k <- split(sigma.true, row(sigma.true))

ind.d <- seq(1,D)
phi.d <- exp(2 + 1/ind.d)      # parameters phi_d  (assumed to be known in application)


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



# Reading theoretical MSEs for the plug-in predictors (need to be calculated from Simulation 2)

#load("Results/REd_D=100.Rdata")

mse.teo.mud.in <- REd$REd.mud.in^2
mse.teo.Rd.in <- REd$REd.Rd.in^2



### SIMULATION STARTS

BB <- 50
ITER <- 2   # number of iterations

param.est <- list() 
mud.true <- mud.in <- list()
Rd.true <- Rd.in  <- list() 
mse.mud.in <- mse.Rd.in <- list()




###############################################
# start of the simulation loop (iter)

iter <- 1
no.convergence <- 0
no.convergence.b <- rep(0,BB)



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
  
  mud.true[[iter]] <- mu
  Rd.true[[iter]] <- mu.k[[2]]/(mu.k[[1]]+mu.k[[2]])
  
  ## generating y_d
  
  alpha <- phi.d * mu
  
  yy <- rdirichlet(D, alpha)
  
  
  
  ### fitting of the model with MLE Laplace
  
  # seeds for parameters 
  beta.s <- beta.k
  sigma.s <-  unlist(sigma.k)
  
  seeds.u <- split(matrix(0, D, 2), 1:D)   # seeds for u_dk
  
  res.MLE <- MLE_Laplace(XX.k, yy, phi.d, beta.s, sigma.s, seeds.u, maxiter=30, epsilon.u=0.0001, epsilon.par=0.00001)
  
  
  if (res.MLE$convergence == FALSE){
    iter <- iter - 1
    no.convergence <- no.convergence + 1
  }
  else {
    
    param.est.MLE <- res.MLE$par.est
    
    
    ### Calculating Plug-in
    ud.mode <- res.MLE$udk.pred
    ud.mod <- split(Reduce(ud.mode, f=cbind),1:(q-1))
    
    plug <- plug.in(param.est.MLE, ud.mod, XX.k, yy,  phi.d)
    
    mud.in[[iter]] <- Reduce(plug$mud.in, f=cbind)
    Rd.in[[iter]] <- plug$Rd.in
    
    Boot.DIRICH <-  Boot.Dirichlet(XX.k, param.est.MLE,  MAXITER = 100, precision = 10^-4)

    mse.mud.in[[iter]] <- Boot.DIRICH[[2]]
    mse.Rd.in[[iter]] <- Boot.DIRICH[[1]]
    no.convergence.b[iter] <- Boot.DIRICH[[3]]
    
  }
  
  cat("Iteration overall:  ", iter, "\n")
  
  iter <- iter + 1
  
}  ### cycle iter




#### Bd
Bd.boot.Rd.in <- (Reduce(lapply(mse.Rd.in, FUN="-", mse.teo.Rd.in), f="+"))/ITER
Bd.boot.mud.in <- (Reduce(lapply(mse.mud.in, FUN="-", mse.teo.mud.in), f="+"))/ITER

##REd

REd.boot.Rd.in <- (Reduce(lapply(lapply(mse.Rd.in, FUN="-", mse.teo.Rd.in),FUN="^",2), f="+")/ITER)^0.5
REd.boot.mud.in <- (Reduce(lapply(lapply(mse.mud.in, FUN="-", mse.teo.mud.in), FUN="^",2), f="+")/ITER)^0.5
  
#### AB

AB.boot.R.in <-  mean(abs(Bd.boot.Rd.in))
AB.boot.mu.in <- colMeans(abs(Bd.boot.mud.in))

RE.boot.Rd.in <-  mean(REd.boot.Rd.in)
RE.boot.mud.in <- colMeans(REd.boot.mud.in)


#### RREd

RREd.boot.Rd.in <- REd.boot.Rd.in/mse.teo.Rd.in*100
RREd.boot.mud.in <- REd.boot.mud.in/mse.teo.mud.in*100

#### RBd

RBd.boot.Rd.in <-  Bd.boot.Rd.in/mse.teo.Rd.in*100
RBd.boot.mud.in <- Bd.boot.mud.in/mse.teo.mud.in*100


#### RRE

RRE.boot.Rd.in <-  mean(RREd.boot.Rd.in)
RRE.boot.mud.in <- colMeans(RREd.boot.mud.in)


#### RAB

RAB.boot.Rd.in <-  mean(abs(RBd.boot.Rd.in))
RAB.boot.mud.in <- colMeans(abs(RBd.boot.mud.in))


#### Saving results 

AB <- c(AB.boot.mu.in[1], AB.boot.mu.in[2], AB.boot.mu.in[3], AB.boot.R.in)
RAB <- c(RAB.boot.mud.in[1], RAB.boot.mud.in[2],  RAB.boot.mud.in[3], RAB.boot.Rd.in)

RE <- c(RE.boot.mud.in[1], RE.boot.mud.in[2], RE.boot.mud.in[3], RE.boot.Rd.in)
RRE <- c(RRE.boot.mud.in[1], RRE.boot.mud.in[2], RRE.boot.mud.in[3], RRE.boot.Rd.in)

results <- data.frame(AB, RE, RAB, RRE)

rownames(results) <- c("mud_1_in",  "mud_2_in",  "mud_3_in", "Rd.in")









#write.table(results,file=paste0("Results/Sim3_D=", D, " B=", BB,".txt"), sep="\t")


#save(list=ls(), file=paste("Results/Sim3_D=", D, " B=", BB, ".Rdata", sep=""))


