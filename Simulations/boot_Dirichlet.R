Boot.Dirichlet <- function(XX.k, param.e,  MAXITER = 100, precision = 10^-4){
  
  
beta.b <- sigma.b <- list()
beta.b[[1]] <- param.e[1:p1]
beta.b[[2]] <- param.e[(p1+1):p]
sigma.b <- param.e[(p+1):(p+2)]
sigma.k.b <- as.list(sigma.b)


mud.true.b <- Rd.true.b <- mud.in.b <- Rd.in.b <- list()  
iter.b <- 1
no.convergence <- 0

while (iter.b <= BB) {

udk.b <-  split(rmvnorm(D, rep(0, (q-1))), 1:(q-1))

sigma.ud.b <- Map(udk.b, f="*", sigma.k.b)

eta.k.b <- Map(Map(XX.k, f="%*%", beta.b), sigma.ud.b, f="+")

mu.k.b <- list()

eta.b <- exp(Reduce(eta.k.b, f=cbind))

sum.eta.b <- rowSums(eta.b)

mu.k.b[[q]] <- 1/(1+sum.eta.b)

for (k in 1:(q-1)){
  mu.k.b[[k]] <- eta.b[,k]/(1+sum.eta.b)
}

mu.b <- Reduce(mu.k.b, f=cbind )    # matrix of mu_dk

mud.true.b[[iter.b]] <- mu.b
Rd.true.b[[iter.b]] <- mu.k.b[[2]]/(mu.k.b[[1]] + mu.k.b[[2]])

## generating y_d

alpha.b <- phi.d * mu.b

yy.b <- rdirichlet(D, alpha.b)

seeds.u <- split(matrix(0, D, 2), 1:D)   # seeds for u_dk

res.MLE.b <- MLE_Laplace(XX.k, yy.b, phi.d, beta.b, sigma.b, seeds.u, maxiter=30, epsilon.u=0.0001, epsilon.par=0.00001)

if (res.MLE.b$convergence == FALSE){
  iter.b <- iter.b - 1
  no.convergence <- no.convergence + 1
}
else{
  
  param.est.MLE.b <- res.MLE.b$par.est
  
  ### calculating plug-in
  ud.mode.b <- res.MLE.b$udk.pred
  ud.mod.b <- split(Reduce(ud.mode.b, f=cbind),1:(q-1))
  
  plug.b <- plug.in(param.est.MLE.b, ud.mod.b, XX.k, yy.b,  phi.d)
  
  mud.in.b[[iter.b]] <- Reduce(plug.b$mud.in, f=cbind)
  Rd.in.b[[iter.b]] <- plug.b$Rd.in
 }

cat("Iteration overall_Boot:  ", "iter:", iter, " iter.b:", iter.b, "\n")

iter.b <- iter.b + 1

}  ### cycle iter


REd.Rd.in.b <- ((Reduce(lapply(Map(Rd.in.b, f="-", Rd.true.b),FUN="^",2),f="+"))/BB)

REd.mud.in.b <- ((Reduce(lapply(Map(mud.in.b, f="-", mud.true.b),FUN="^",2),f="+"))/BB)


results <- list(REd.Rd.in.b, REd.mud.in.b, no.convergence)

return (results)
}