##########################################################
### functions for calculating predictors under the ADMM model
### 
### Autor: Tomas
### Version: June 2024


### Plug-in predictors 
#   (formulas (6.3), (6.4) or step 2.6 of Simulation 2)

plug.in <- function(param.e, ud.m, XX.k, yy,  phi.d ){
 
  beta.e <- list() 
  beta.e[[1]] <- param.e[1:p1]
  beta.e[[2]] <- param.e[(p1+1):p]
  sigma.e <- param.e[(p+1):(p+2)]
  sigma.e <- as.list(sigma.e)
  
  # calculating mu_dk
  
  sigma.ud <- Map(ud.m, f="*", sigma.e)
  
  eta.k <- Map(Map(XX.k, f="%*%", beta.e), sigma.ud, f="+")
  
  mu.in <- list()
  
  eta <- exp(Reduce(eta.k, f=cbind))
  
  sum.eta <- rowSums(eta)
  
  mu.in[[q]] <- 1/(1+sum.eta)
  
  for (k in 1:(q-1)){
    mu.in[[k]] <- eta[,k]/(1+sum.eta)
  }
  
  Rd.in <- mu.in[[2]]/(mu.in[[1]]+mu.in[[2]])
  
  return(list(mud.in = mu.in, Rd.in = Rd.in))
  
}

##############################################
# EBP predictors 
# formula (6.2)  or step 2.5 of Simulation 2)

ebp <- function(param.e, XX.k, yy,  phi.d, S){
 
  # generating random effects
  udk.s <- list()
  
  for (s in 1:S) {
    aux <- split(rmvnorm(D, rep(0, (q-1))), 1:(q-1))  
    udk.s[[s]] <- aux 
    udk.s[[s+S]] <- lapply(aux,"*",-1)
  }
  
  mu.s <- Rd.s <- Dd.s <- list()
  
  for (s in 1:(2*S)) {
    
    ud.s <- udk.s[[s]]
    
    # calculating quantities of step 2.5 ii. using the already implemented 
    # function plugin
    plug <- plug.in(param.est.MLE, ud.s, XX.k, yy,  phi.d)    
    
    mud.s <- Reduce(plug$mud.in, f=cbind)
    Rd.s[[s]] <- plug$Rd.in
    
    mu.s[[s]] <- mud.s
    
    alpha.s <- phi.d * mud.s
    
    log.y <- log(yy)
    log.g <- log(gamma(alpha.s))
    
    Dd.s[[s]] <- exp(rowSums((alpha.s-1)*log.y - log.g))
    
  }
  
  Dd <- Reduce(`+`, Dd.s)     # sum of D_d^s over s
  
  Rd.Dd <- Map(Rd.s, f="*", Dd.s)
  AdR <- Reduce(`+`, Rd.Dd)
  
  mud.Dd <- Map(mu.s, f="*", Dd.s)
  Adk <- Reduce(`+`, mud.Dd)
  
  mud.ebp <- Adk / Dd
  Rd.ebp <- AdR / Dd
  
  muk.ebp <- split(mud.ebp, rep(1:q, each = D))
  
  return(list(mud.ebp = muk.ebp, Rd.ebp = Rd.ebp))
  
}