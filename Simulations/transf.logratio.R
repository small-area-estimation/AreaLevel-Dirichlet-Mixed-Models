########################################################################
########################################################################
######################## LOG-RATIO TRANSFORMATIONS ####################
########################################################################
#######################################################################



#### Additive (ALR)

transf.alr <- function(df, V){
  
  ### df: data frame 
  ### V : VeT (sampling error matrix)
  
  R <- length(df)
  alr.y <- log(apply(df[,-R], 2, FUN = "/", df[,R])) 
  
  #### Variance-Covariance matrix
  col1 <- t(matrix(rep(1,R-1), nrow = R-1))
  
  H0 <- R*(diag(1,nrow = (R-1))+t(col1)%*%col1)
  
  VeTcomp <- lapply(V,function(x) H0 %*% x %*% t(H0))
  cov.meanZ <- mapply(function(i,j) unlist(lapply(VeTcomp, function(x)  x[i,j])), rep(1:(R-1), each = R-1), rep(1:(R-1), times = R-1))
  colnames(cov.meanZ) <- paste0("cov.mean",rep(1:(R-1), each = R-1), rep(1:(R-1), times = R-1))

  return(list(alr.y = alr.y, Ved.alr = VeTcomp, cov.mean.alr = cov.meanZ))
  
  
}

transf.inv.alr <- function(df){
  
  ### df: data frame 
  
  inv.alr.y <- apply(exp(df), 2, function(x) x /(1 + apply(exp(df),1,sum)))
  
  X3 <- 1 - apply(inv.alr.y,1,sum) ## Adding last sector
  df <- data.frame(inv.alr.y,X3)
  
  return(df)
  
}




#### Isometric (ILR)

transf.ilr <- function(df, V){
  
  ### df: data frame 
  ### V : VeT (sampling error matrix)
  
  R <- length(df)
  gz <- (apply(df,1,prod))^(1/R) ## Funcion g(z)
  
  clr.z <- log(apply(df,2,FUN = "/", gz)) ### Transformamos a clr+(z+)
  
  phi <- matrix(0, ncol = R, nrow = R-1) ## Matriz phi
  
  for(i in 1:(R-1)){
    for(j in 1:R)
      if(i+j <= R)
        phi[i,j] <- 1/(sqrt((R-i)*(R-i+1)))
    else if (i+j == R+1) 
      phi[i,j] <- -sqrt(R-i)/(sqrt(R-i+1))
    else 
      phi[i,j] <- 0
  }
  
  ilr.y <- t(apply(clr.z,1,function(x) phi%*%x)) ### Transformamos a y_i (dimension vector R-1)
  colnames(ilr.y) <- paste0("mean.Sector",1:(R-1))
  
  H0 <- R*phi[,1:(R-1)]
  
  
  VeTcomp <- lapply(V,function(x) H0%*%x%*%t(H0))
  cov.mean.ilr.y <- mapply(function(i,j) unlist(lapply(VeTcomp, function(x)  x[i,j])), rep(1:(R-1), each = R-1), rep(1:(R-1), times = R-1))
  colnames(cov.mean.ilr.y) <- paste0("cov.meanSector",rep(1:(R-1), each = R-1), rep(1:(R-1), times = R-1))
  cov.mean.ilr.y <- cov.mean.ilr.y[,c(1,8,15,22,29,36,2:6,9:12,16:18,23:24,30)]
  
  return(list(ilr_y = ilr.y, VeT.irl = VeTcomp, cov.mean.ilr = cov.mean.ilr.y, phi = phi))
  
}



transf.inv.ilr <- function(df, phi){
  
  phi.t.y <- t(apply(df,1,function(x) t(phi)%*%x))
  cy <- apply(exp(phi.t.y),1,sum)  ### Funcion c(y) = sum(exp(y_i))
  transf <- apply(exp(phi.t.y),2,FUN = "/", cy) ## Deshacemos el cambio (clr^-1)
  
  colnames(transf) <- paste0("X",1:7)
  
  return(transf)
  
}





#### Centered (CLR)

transf.clr <- function(df, V){
  
  ### df: data frame 
  ### V : VeT (sampling error matrix)
  
  R <- length(df)
  
  gz <- (apply(df,1,prod))^(1/R)
  
  clr.y <- log(apply(df[,-R],2,FUN = "/", gz))
  
  
  #### Variance-Covariance matrix
  
  H0 <- R*diag(1,nrow = (R-1))
  
  VeTcomp <- lapply(V,function(x) H0%*%x%*%t(H0))
  
  cov.meanZ <- mapply(function(i,j) unlist(lapply(VeTcomp, function(x)  x[i,j])), rep(1:(R-1), each = R-1), rep(1:(R-1), times = R-1))
  colnames(cov.meanZ) <- paste0("cov.meanSector",rep(1:(R-1), each = R-1), rep(1:(R-1), times = R-1))
  cov.meanZ <- cov.meanZ[,c(1,8,15,22,29,36,2:6,9:12,16:18,23:24,30)]
  
  return(list(clr.y = clr.y, VeT.crl = VeTcomp, cov.mean.crl = cov.meanZ))
  
}



transf.inv.clr <- function(df){
  
  dflast <- -apply(df,1,sum) ### y_q = -y_1-...-y_(q-1)
  df <- cbind(df,dflast) ### Unimos en el mismo data frame
  cy <- apply(exp(df),1,sum)  ### Funcion c(y) = sum(exp(y_i))
  df <- apply(exp(df),2,FUN = "/", cy) ## Deshacemos el cambio
  
  colnames(df) <- paste0("X",1:7)
  
  return(df)
  
}
