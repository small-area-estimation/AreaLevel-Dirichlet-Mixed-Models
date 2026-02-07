###############################################################################
###
###        Simulation experiment 2  Dirichlet regression ADMM, BFH-alr and MMM           
###                                                            
###                            January 2026                         
###
### Autor: Tomas and Esteban




rm(list = ls())


source("functions_MLE_Laplace.R")
source("functions_predictors.R")
source("functions.bfh.R")
source("transf.logratio.R")

library(lme4)
library(Matrix)
library(mvtnorm)
library(cubature)
library(DirichletReg)
library(sae)
library(mme)

library(nloptr)


opts <- list( "algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=10000)

set.seed(17032005)

D <- 50     # number of areas
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

ITER <- 10^3   # number of iterations

param.est <- list() 
mud.true <- mud.in <- mud.ebp <- list()
mud.in.bfh <- pred.MMM <- Rd.in.MMM <- Rd.in.bfh <- list()
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
  
  ## generating sigma^2_edk (BFH-alr)
  Ved <- vector("list",D)
  
  for(d in 1:D){
    Sigma_d <- matrix(NA, nrow = 2, ncol = 2)
    
    off_diag <- -mu[,1][d]*mu[,2][d]/(1 + phi.d[d]) 
    Sigma_d[upper.tri(Sigma_d)] <- off_diag
    Sigma_d[lower.tri(Sigma_d)] <- off_diag
    diag(Sigma_d) <- c(mu[,1][d]*(1-mu[,1][d]), mu[,2][d]*(1-mu[,2][d]))/(1+phi.d[d]) 

    Ved[[d]] <- Sigma_d
  }
  

  
  ## generating y_d
  
  alpha <- phi.d * mu
  
  yy <- rdirichlet(D, alpha)
  
  ## generating y_d (BFH-alr)
  y.alr <- transf.alr(df = data.frame(yy), V = Ved)
  y1.alr <- y.alr$alr.y[,1]
  y2.alr <- y.alr$alr.y[,2]
  
  covZ <- y.alr$cov.mean.alr
  Ved.alr <- y.alr$Ved.alr
  
  ## Fitting FH model for obtaining seeds for variance components
  vardiry1.alr <- covZ[,1]; vardiry2.alr <- covZ[,4]
  

  fit.FH.sae1 <- eblupFH(y1.alr ~ XX.k[[1]][,2], vardiry1.alr) # FH model for y1.alr
  fit.FH.sae2 <- eblupFH(y2.alr ~ XX.k[[2]][,2], vardiry2.alr) # FH model for y2.alr
  
  
  # Seeds for variance component parameters
  thetas.0 <- c(fit.FH.sae1$fit$refvar, fit.FH.sae2$fit$refvar, 0)
  Vud <- UveU(thetas.0)
  
  y.alr <- lapply(1:D,
                  function(d) matrix(c(y1.alr[d], y2.alr[d])))
  X <- lapply(1:D,
              function(d) as.matrix(t(bdiag(as.numeric(c(1,XX.k[[1]][d,2])),
                                            as.numeric(c(1,XX.k[[2]][d,2]))))))
  
  fit.BFH <- try(REML.BFH(X, y.alr, D, Ved.alr, Vud), TRUE)
  
 #### fitting the MMM model
  
  nd <- ceiling(1 + phi.d)
  ndk <- matrix(0, nrow = D, ncol = q)
  
  for(d in 1:D){
    raw_counts <- yy[d, ] * nd[d]
    int_counts <- round(raw_counts)
    diff <- nd[d] - sum(int_counts)
    if(diff != 0){
      int_counts[q] <- int_counts[q] + diff
    }
    
    ndk[d, ] <- int_counts
  }
  
  pp <- c(1,1) #vector with the number of auxiliary variables in each category
  time <- rep(1,D)
  
  mult.data.eff <- data.frame(Area = 1:D, Time = time, Sample = nd, 
                              Population = rep(10000, D), 
                              eff.Y1 = ndk[,1], eff.Y2 = ndk[,2], eff.Y3 = ndk[,3],
                              X1 = XX.k[[1]][,2], X2 = XX.k[[2]][,2]
  )
  
  mod <- 1 #Type of model
  
  datar <- data.mme(mult.data.eff, q, pp, mod)
  
  mult.eff.fit <- modelfit1(pp, Xk = datar$Xk, X = datar$X, Z = datar$Z, 
                            initial = datar$initial, y = datar$y[,1:(q-1)],
                            M = datar$n, MM = datar$N)
  
  
  ### fitting of the ADMM with MLE Laplace
  
  # seeds for parameters 
  beta.s <- beta.k
  sigma.s <-  unlist(sigma.k)
  
  seeds.u <- split(matrix(0, D, 2), 1:D)   # seeds for u_dk
  
  res.MLE <- MLE_Laplace(XX.k, yy, phi.d, beta.s, sigma.s, seeds.u, maxiter=30, epsilon.u=0.0001, epsilon.par=0.00001)
  
  param.est.MLE <- res.MLE$par.est
  
  # param.est[[iter]] <- res.MLE$par.est
  
  if (res.MLE$convergence == FALSE || inherits(fit.BFH, "try-error") || inherits(mult.eff.fit, "try-error")){
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
    
    ### calculation EBLUPs and plug-in rates (BFH-alr)
    beta.u.hat <- BETA.U.BFH(X, y.alr, D, Ved.alr, fit.BFH[[1]])
    beta.hat <- beta.u.hat[[1]] # Estimates of regression parameters
    theta.hat <- fit.BFH[[1]] # Estimates of variance components
    u <- beta.u.hat[[2]] # Predictions of random effects
    
    eblup.bfh <- t(mapply(FUN="+", lapply(X, FUN="%*%", beta.hat), u))
    
    eblup.inv.alr <- transf.inv.alr(eblup.bfh)
    
    eblup.inv.alr.1 <- eblup.inv.alr[,1] # EBLUPs of employed people
    eblup.inv.alr.2 <- eblup.inv.alr[,2] # EBLUPs of unemployed people
    eblup.inv.alr.3 <- eblup.inv.alr[,3] # EBLUPs of inactive people
    
    mud.in.bfh[[iter]] <- Reduce(eblup.inv.alr, f=cbind)
    Rd.in.bfh[[iter]] <- eblup.inv.alr.2/(eblup.inv.alr.1 + eblup.inv.alr.2)
    
    ### calculating probabilites and plug-in rates (MMM)
    pred.MMM[[iter]] <- mult.eff.fit$Estimated.probabilities
    Rd.in.MMM[[iter]] <- mult.eff.fit$Estimated.probabilities[,2]/(mult.eff.fit$Estimated.probabilities[,1] + mult.eff.fit$Estimated.probabilities[,2])
    
  }
  
  cat("Iteration overall:  ", iter, "\n")
  
  iter <- iter + 1
  
}  ### cycle iter






#### Calculating output (Step 3 of the simulation 2)

#### B_d

ADMM.Bd.Rd.in <- (Reduce( Map(Rd.in, f="-", Rd.true), f="+"))/ITER
BFHalr.Bd.Rd.in <- (Reduce( Map(Rd.in.bfh, f="-", Rd.true), f="+"))/ITER
MMM.Bd.Rd.in <- (Reduce( Map(Rd.in.MMM, f="-", Rd.true), f="+"))/ITER

ADMM.Bd.mud.in <- (Reduce( Map(mud.in, f="-", mud.true), f="+"))/ITER
BFHalr.Bd.mud.in <- (Reduce( Map(mud.in.bfh, f="-", mud.true), f="+"))/ITER
MMM.Bd.mud.in <- (Reduce( Map(pred.MMM, f="-", mud.true), f="+"))/ITER

#### RE_d

ADMM.REd.Rd.in <- ((Reduce(lapply(Map(Rd.in, f="-", Rd.true),FUN="^",2),f="+"))/ITER)^0.5
BFHalr.REd.Rd.in <- ((Reduce(lapply(Map(Rd.in.bfh, f="-", Rd.true),FUN="^",2),f="+"))/ITER)^0.5
MMM.REd.Rd.in <- ((Reduce(lapply(Map(Rd.in.MMM, f="-", Rd.true),FUN="^",2),f="+"))/ITER)^0.5

ADMM.REd.mud.in <- ((Reduce(lapply(Map(mud.in, f="-", mud.true),FUN="^",2),f="+"))/ITER)^0.5
BFHalr.REd.mud.in <- ((Reduce(lapply(Map(mud.in.bfh, f="-", mud.true),FUN="^",2),f="+"))/ITER)^0.5
MMM.REd.mud.in <- ((Reduce(lapply(Map(pred.MMM, f="-", mud.true),FUN="^",2),f="+"))/ITER)^0.5

#### eta_d

etad.Rd <- (Reduce(Rd.true, f="+"))/ITER
etad.mud <- (Reduce(mud.true, f="+"))/ITER

#### RE

ADMM.RE.Rd.in <-  mean(ADMM.REd.Rd.in)
BFHalr.RE.Rd.in <-  mean(BFHalr.REd.Rd.in)
MMM.RE.Rd.in <-  mean(MMM.REd.Rd.in)

ADMM.RE.mud.in <- colMeans(ADMM.REd.mud.in)
BFHalr.RE.mud.in <- colMeans(BFHalr.REd.mud.in)
MMM.RE.mud.in <- colMeans(MMM.REd.mud.in)

#### AB

ADMM.AB.Rd.in <-  mean(abs(ADMM.Bd.Rd.in))
BFHalr.AB.Rd.in <-  mean(abs(BFHalr.Bd.Rd.in))
MMM.AB.Rd.in <-  mean(abs(MMM.Bd.Rd.in))


ADMM.AB.mud.in <- colMeans(abs(ADMM.Bd.mud.in))
BFHalr.AB.mud.in <- colMeans(abs(BFHalr.Bd.mud.in))
MMM.AB.mud.in <- colMeans(abs(MMM.Bd.mud.in))

#### RRE_d

ADMM.RREd.Rd.in <- ADMM.REd.Rd.in/abs(etad.Rd)*100
BFHalr.RREd.Rd.in <- BFHalr.REd.Rd.in/abs(etad.Rd)*100
MMM.RREd.Rd.in <- MMM.REd.Rd.in/abs(etad.Rd)*100


ADMM.RREd.mud.in <- ADMM.REd.mud.in/abs(etad.mud)*100
BFHalr.RREd.mud.in <- BFHalr.REd.mud.in/abs(etad.mud)*100
MMM.RREd.mud.in <- MMM.REd.mud.in/abs(etad.mud)*100

#### RB_d

ADMM.RBd.Rd.in <- ADMM.Bd.Rd.in/abs(etad.Rd)*100
BFHalr.RBd.Rd.in <- BFHalr.Bd.Rd.in/abs(etad.Rd)*100
MMM.RBd.Rd.in <- MMM.Bd.Rd.in/abs(etad.Rd)*100


ADMM.RBd.mud.in <- ADMM.Bd.mud.in/abs(etad.mud)*100
BFHalr.RBd.mud.in <- BFHalr.Bd.mud.in/abs(etad.mud)*100
MMM.RBd.mud.in <- MMM.Bd.mud.in/abs(etad.mud)*100

#### RRE

ADMM.RRE.Rd.in <-  mean(ADMM.RREd.Rd.in)
BFHalr.RRE.Rd.in <-  mean(BFHalr.RREd.Rd.in)
MMM.RRE.Rd.in <-  mean(MMM.RREd.Rd.in)


ADMM.RRE.mud.in <- colMeans(ADMM.RREd.mud.in)
BFHalr.RRE.mud.in <- colMeans(BFHalr.RREd.mud.in)
MMM.RRE.mud.in <- colMeans(MMM.RREd.mud.in)


#### RAB

ADMM.RAB.Rd.in <-  mean(abs(ADMM.RBd.Rd.in))
BFHalr.RAB.Rd.in <-  mean(abs(BFHalr.RBd.Rd.in))
MMM.RAB.Rd.in <-  mean(abs(MMM.RBd.Rd.in))


ADMM.RAB.mud.in <- colMeans(abs(ADMM.RBd.mud.in))
BFHalr.RAB.mud.in <- colMeans(abs(BFHalr.RBd.mud.in))
MMM.RAB.mud.in <- colMeans(abs(MMM.RBd.mud.in))


#### Saving results 
AB <- c(ADMM.AB.mud.in[1], BFHalr.AB.mud.in[1], MMM.AB.mud.in[1],
               ADMM.AB.mud.in[2], BFHalr.AB.mud.in[2], MMM.AB.mud.in[2],
               ADMM.AB.mud.in[3], BFHalr.AB.mud.in[3], MMM.AB.mud.in[3],
               ADMM.AB.Rd.in, BFHalr.AB.Rd.in, MMM.AB.Rd.in)

RAB <- c(ADMM.RAB.mud.in[1], BFHalr.RAB.mud.in[1], MMM.RAB.mud.in[1],
                ADMM.RAB.mud.in[2], BFHalr.RAB.mud.in[2], MMM.RAB.mud.in[2],
                ADMM.RAB.mud.in[3], BFHalr.RAB.mud.in[3], MMM.RAB.mud.in[3],
                ADMM.RAB.Rd.in, BFHalr.RAB.Rd.in, MMM.RAB.Rd.in)

RE <- c(ADMM.RE.mud.in[1], BFHalr.RE.mud.in[1], MMM.RE.mud.in[1],
               ADMM.RE.mud.in[2], BFHalr.RE.mud.in[2], MMM.RE.mud.in[2],
               ADMM.RE.mud.in[3], BFHalr.RE.mud.in[3], MMM.RE.mud.in[3],
               ADMM.RE.Rd.in, BFHalr.RE.Rd.in, MMM.RE.Rd.in)

RRE <- c(ADMM.RRE.mud.in[1], BFHalr.RRE.mud.in[1], MMM.RRE.mud.in[1],
                ADMM.RRE.mud.in[2], BFHalr.RRE.mud.in[2], MMM.RRE.mud.in[2],
                ADMM.RRE.mud.in[3], BFHalr.RRE.mud.in[3], MMM.RRE.mud.in[3],
                ADMM.RRE.Rd.in, BFHalr.RRE.Rd.in, MMM.RRE.Rd.in)
                                                                     


results <- data.frame(AB, RE, RAB, RRE)

rownames(results) <- c("mud_1_ADMM", "mud_1_BFH-alr", "mud_1_MMM",
                       "mud_2_ADMM", "mud_2_BFH-alr", "mud_2_MMM",
                       "mud_3_ADMM", "mud_3_BFH-alr", "mud_3_MMM",
                       "Rd_ADMM", "Rd_BFH-alr", "Rd_MMM")


#write.table(results, file=paste0("Results/Sim2_Models_D=", D, ".txt"), sep="\t")


#save(list=ls(), file=paste("Results/Sim2_Models_D=", D, ".Rdata", sep=""))


#REd <- list(REd.Rd.in=REd.Rd.in, REd.Rd.ebp=REd.Rd.ebp, REd.Rd.in.ebp=REd.Rd.in.ebp, REd.mud.in=REd.mud.in, REd.mud.ebp=REd.mud.ebp)


#save(REd, file=paste("Results/REd_D=", D, ".Rdata", sep=""))

