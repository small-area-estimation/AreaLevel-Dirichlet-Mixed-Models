###############################################################################
###
###        Simulation experiment 2  Dirichlet regression ADMM           
###                                                            
###                            June 2024                         
###
### Author: Tomas Hobza




rm(list = ls())


source("functions_MLE_Laplace.R")
source("functions_predictors.R")


library(lme4)
library(Matrix)
library(mvtnorm)
library(cubature)
library(DirichletReg)
library(nloptr)



opts <- list( "algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=10000)


set.seed(17032005)

D <- 100     # number of areas
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


### SIMULATION STARTS

ITER <- 1000   # number of iterations

param.est <- list() 
mud.true <- mud.in <- mud.ebp <- list()
Rd.true <- Rd.in <- Rd.ebp <- Rd.in.ebp <- list() 


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
  
  param.est.MLE <- res.MLE$par.est
  
 # param.est[[iter]] <- res.MLE$par.est
  
  if (res.MLE$convergence == FALSE){
    iter <- iter - 1
    no.convergence <- no.convergence + 1
  }
  else {
    
    ### calculating plug-in
    ud.mode <- res.MLE$udk.pred
    ud.mod <- split(Reduce(ud.mode, f=cbind),1:(q-1))
    
    plug <- plug.in(param.est.MLE, ud.mod, XX.k, yy,  phi.d)
    
    mud.in[[iter]] <- Reduce(plug$mud.in, f=cbind)
    Rd.in[[iter]] <- plug$Rd.in
    
    ### calculating EBP
    
    S <-500
    
    eb <- ebp(param.est.MLE, XX.k, yy,  phi.d, S)
    
    mud.e <- eb$mud.ebp
    
    mud.ebp[[iter]] <- Reduce(mud.e, f=cbind)
    Rd.ebp[[iter]] <- eb$Rd.ebp
    
    Rd.in.ebp[[iter]] <- mud.e[[2]]/(mud.e[[1]]+mud.e[[2]])
    
  }
  
  cat("Iteration overall:  ", iter, "\n")
  
  iter <- iter + 1
  
  }  ### cycle iter





#### Calculating output (Step 3 of the simulation 2)

#### B_d

Bd.Rd.in <- (Reduce( Map(Rd.in, f="-", Rd.true), f="+"))/ITER
Bd.Rd.ebp <- (Reduce( Map(Rd.ebp, f="-", Rd.true), f="+"))/ITER
Bd.Rd.in.ebp <- (Reduce( Map(Rd.in.ebp, f="-", Rd.true), f="+"))/ITER

Bd.mud.in <- (Reduce( Map(mud.in, f="-", mud.true), f="+"))/ITER
Bd.mud.ebp <- (Reduce( Map(mud.ebp, f="-", mud.true), f="+"))/ITER

#### RE_d

REd.Rd.in <- ((Reduce(lapply(Map(Rd.in, f="-", Rd.true),FUN="^",2),f="+"))/ITER)^0.5
REd.Rd.ebp <- ((Reduce(lapply(Map(Rd.ebp, f="-", Rd.true),FUN="^",2),f="+"))/ITER)^0.5
REd.Rd.in.ebp <- ((Reduce(lapply(Map(Rd.in.ebp, f="-", Rd.true),FUN="^",2),f="+"))/ITER)^0.5

REd.mud.in <- ((Reduce(lapply(Map(mud.in, f="-", mud.true),FUN="^",2),f="+"))/ITER)^0.5
REd.mud.ebp <- ((Reduce(lapply(Map(mud.ebp, f="-", mud.true),FUN="^",2),f="+"))/ITER)^0.5

#### eta_d

etad.Rd <- (Reduce(Rd.true, f="+"))/ITER
etad.mud <- (Reduce(mud.true, f="+"))/ITER

#### RE

RE.Rd.in <-  mean(REd.Rd.in)
RE.Rd.ebp <- mean(REd.Rd.ebp)
RE.Rd.in.ebp <- mean(REd.Rd.in.ebp)

RE.mud.in <- colMeans(REd.mud.in)
RE.mud.ebp <- colMeans(REd.mud.ebp)

#### AB

AB.Rd.in <-  mean(abs(Bd.Rd.in))
AB.Rd.ebp <- mean(abs(Bd.Rd.ebp))
AB.Rd.in.ebp <- mean(abs(Bd.Rd.in.ebp))

AB.mud.in <- colMeans(abs(Bd.mud.in))
AB.mud.ebp <- colMeans(abs(Bd.mud.ebp))

#### RRE_d

RREd.Rd.in <- REd.Rd.in/abs(etad.Rd)*100
RREd.Rd.ebp <- REd.Rd.ebp/abs(etad.Rd)*100
RREd.Rd.in.ebp <- REd.Rd.in.ebp/abs(etad.Rd)*100

RREd.mud.in <- REd.mud.in/abs(etad.mud)*100
RREd.mud.ebp <- REd.mud.ebp/abs(etad.mud)*100

#### RB_d

RBd.Rd.in <- Bd.Rd.in/abs(etad.Rd)*100
RBd.Rd.ebp <- Bd.Rd.ebp/abs(etad.Rd)*100
RBd.Rd.in.ebp <- Bd.Rd.in.ebp/abs(etad.Rd)*100

RBd.mud.in <- Bd.mud.in/abs(etad.mud)*100
RBd.mud.ebp <- Bd.mud.ebp/abs(etad.mud)*100

#### RRE

RRE.Rd.in <-  mean(RREd.Rd.in)
RRE.Rd.ebp <- mean(RREd.Rd.ebp)
RRE.Rd.in.ebp <- mean(RREd.Rd.in.ebp)

RRE.mud.in <- colMeans(RREd.mud.in)
RRE.mud.ebp <- colMeans(RREd.mud.ebp)

#### RAB

RAB.Rd.in <-  mean(abs(RBd.Rd.in))
RAB.Rd.ebp <- mean(abs(RBd.Rd.ebp))
RAB.Rd.in.ebp <- mean(abs(RBd.Rd.in.ebp))

RAB.mud.in <- colMeans(abs(RBd.mud.in))
RAB.mud.ebp <- colMeans(abs(RBd.mud.ebp))


#### Saving results 

AB <- c(AB.mud.ebp[1], AB.mud.in[1], AB.mud.ebp[2], AB.mud.in[2], AB.mud.ebp[3], AB.mud.in[3], AB.Rd.ebp, AB.Rd.in.ebp, AB.Rd.in)
RAB <- c(RAB.mud.ebp[1], RAB.mud.in[1], RAB.mud.ebp[2],  RAB.mud.in[2], RAB.mud.ebp[3], RAB.mud.in[3], RAB.Rd.ebp, RAB.Rd.in.ebp, RAB.Rd.in)

RE <- c(RE.mud.ebp[1], RE.mud.in[1], RE.mud.ebp[2], RE.mud.in[2], RE.mud.ebp[3], RE.mud.in[3], RE.Rd.ebp, RE.Rd.in.ebp, RE.Rd.in)
RRE <- c(RRE.mud.ebp[1], RRE.mud.in[1], RRE.mud.ebp[2], RRE.mud.in[2], RRE.mud.ebp[3], RRE.mud.in[3], RRE.Rd.ebp, RRE.Rd.in.ebp, RRE.Rd.in)

results <- data.frame(AB, RE, RAB, RRE)

rownames(results) <- c("mud_1_ebp", "mud_1_in", "mud_2_ebp", "mud_2_in", "mud_3_ebp", "mud_3_in","Rd.ebp", "Rd.in.ebp", "Rd.in")


#write.table(results, file=paste0("Results/Sim2_D=", D, ".txt"), sep="\t")


#save(list=ls(), file=paste("Results/Sim2_D=", D, ".Rdata", sep=""))


REd <- list(REd.Rd.in=REd.Rd.in, REd.Rd.ebp=REd.Rd.ebp, REd.Rd.in.ebp=REd.Rd.in.ebp, REd.mud.in=REd.mud.in, REd.mud.ebp=REd.mud.ebp)


#save(REd, file=paste("Results/REd_D=", D, ".Rdata", sep=""))




