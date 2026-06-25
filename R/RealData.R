#######################################################################################
###
###                           Example of Application to real data
###                     Area Level Dirichlet Regression 
###                             
###                           
###                      Real data of EPA 2022
###                           
###                             June 2026                           
###
### AUTHOR: Esteban Cabello
### Date of creation: 24 September 2024
### Edited by: Esteban Cabello
### Edited on: June 2026

source("functions_MLE_Laplace.R")
source("functions_predictors.R")
source("boot_Dirichlet.R")


library(openxlsx)
library(nloptr)
library(DirichletReg)
library(dplyr)
library(xtable)

##### Reading data
data <- read.xlsx("real.data.xlsx")
## We only consider the 16-34 and 35-64 age groups (Age1 and Age2)
## The retired group is not of interest (Age3)
data <- data %>% filter(Age <=2)
data$Prov <- as.numeric(data$Prov)


##### Fitting ADMM (Area-level Dirichlet Mixed Model)

q <- 3 
D <- nrow(data)

Intercept <- rep(1,D)

p1 <- 3    #  number of parameters beta1, beta2
p2 <- 3
p <- p1+p2
XX.k <- list(XX.1 = cbind(Intercept, data$edu1.MEAN, data$work1.MEAN), ## High levels of education (edu3) and full-time workers (work1)
             XX.2 = cbind(Intercept, data$edu3.MEAN, data$situ3.MEAN)) ## Low levels of education (edu1) and public sector salaried employees (situ3)


### Dirichlet Regression for fixed effects (necessary to beta.s, sigma.s and seeds.u)


data$Z <- DR_data(data[,10:12], base = 3, trafo = TRUE) 
yy <- data$Z

phi.d <- (1/(q-1)*(data$Z[,1]*(1-data$Z[,1])/data$varY1+(data$Z[,2]*(1-data$Z[,2])/data$varY2))-1)

## Marginal Dirichlet Regressions
reg.d1 <- DirichletReg::DirichReg(formula = Z ~ edu1.MEAN +  work1.MEAN, data, model = "alternative", base = 3)
summary(reg.d1)
reg.d2 <- DirichletReg::DirichReg(formula = Z ~ edu3.MEAN +  situ3.MEAN, data, model = "alternative", base = 3)
summary(reg.d2)


beta.k <- rbind(c(reg.d1$coefficients[1:3]), c(reg.d2$coefficients[4:6]))
beta.s <- split(beta.k, row(beta.k))

sigma.s <- c(0.25,0.25) ### Vector of seeds for the variance components
seeds.u <- split(rep(0, D*(q-1)), 1:D)


opts <- list("algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=50000)

dirichletM <- MLE_Laplace(XX.k, yy, phi.d, beta.s, sigma.s, seeds.u, maxiter=1000, epsilon.u=0.001, epsilon.par=0.001)
dirichletM$par.est
dirichletM$convergence



ud.mode <- dirichletM$udk.pred
ud.mod <- split(Reduce(ud.mode, f=cbind),1:(q-1))



####### Predictors of proportions and rates

pred.in <- plug.in(dirichletM$par.est, ud.mod, XX.k, yy,  phi.d)
mud.in <- Reduce(pred.in$mud.in, f=cbind)

Rd.in  <- pred.in$Rd.in

Rd.y <- yy[,2]/(yy[,2] + yy[,1]) ## Direct estimator of rates

### Standardized Residuals (QQ-plots and Boxplots)

sresiduals <- (yy - mud.in)/sqrt(mud.in*(1-mud.in)/(phi.d+1))


boxplot(sresiduals[,1],sresiduals[,2],sresiduals[,3], 
        main = "sresiduals by Category", cex.main = 1.8, axes = F)
axis(2, col="black", cex.axis = 1.4)
axis(1, at = 1:3, line = 0, col="black", labels = rep("",3))
mtext(text = c("Employed", "Unemployed", "Inactive"), 
      side = 1, at = 1:3, cex = 1.6, line = 2)
abline(h = 0, lwd = 2, lty = 2, col = "red")




library(ggplot2)
library(qqplotr)

sresiduals <- data.frame(apply(sresiduals, 2, as.numeric))

ggplot(sresiduals , aes(sample = Y1)) +
  stat_qq_band(distribution = "norm", bandType = "pointwise",
               fill = NA, color = "red", linewidth = 0.7) +
  # y = x (reference) line in red
  stat_qq_line(distribution = "norm", color = "red", linewidth = 0.9) +
  # hollow points (no fill)
  stat_qq_point(distribution = "norm",
                shape = 1,            # hollow circle
                size = 1.8,
                stroke = 0.9,         # outline width
                color = "black") + 
  labs(x = "Theoretical quantiles (Normal)", y = "Sample quantiles",
       title = "QQ plot sresiduals Employed") +
  theme_classic(base_size = 14)  + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 18))


ggplot(sresiduals , aes(sample = Y2)) +
  stat_qq_band(distribution = "norm", bandType = "pointwise",
               fill = NA, color = "red", linewidth = 0.7) +
  # y = x (reference) line in red
  stat_qq_line(distribution = "norm", color = "red", linewidth = 0.9) +
  # hollow points (no fill)
  stat_qq_point(distribution = "norm",
                shape = 1,            # hollow circle
                size = 1.8,
                stroke = 0.9,         # outline width
                color = "black") + 
  labs(x = "Theoretical quantiles (Normal)", y = "Sample quantiles",
       title = "QQ plot sresiduals Unemployed") +
  theme_classic(base_size = 14)  + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 18))

ggplot(sresiduals , aes(sample = Y3)) +
  stat_qq_band(distribution = "norm", bandType = "pointwise",
               fill = NA, color = "red", linewidth = 0.7) +
  # y = x (reference) line in red
  stat_qq_line(distribution = "norm", color = "red", linewidth = 0.9) +
  # hollow points (no fill)
  stat_qq_point(distribution = "norm",
                shape = 1,            # hollow circle
                size = 1.8,
                stroke = 0.9,         # outline width
                color = "black") + 
  labs(x = "Theoretical quantiles (Normal)", y = "Sample quantiles",
       title = "QQ plot sresiduals Inactive") +
  theme_classic(base_size = 14)  + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 18))





df.pred <- data.frame(cbind(Prov = data$Prov, Age = data$Age,Sex = data$Sex, yy, mud.in, Rd.y, Rd.in, nd = data$nd))
colnames(df.pred) <- c("Prov","Age","Sex","DirectY1","DirectY2", "DirectY3", "Plug.inY1", "Plug.inY2", "Plug.inY3", 
                       "RateDirect", "Rate.in", "nd")






##########################################
############### Bootstrap ################
##########################################

cond <- T
BB <- 5

if(cond == T){ # 
  
  iter <- 1
  library(mvtnorm)
  
  set.seed(03201)
  boot.b <- boot.Dirichlet(XX.k, phi.d, param.e = dirichletM$par.est,  D, MAXITER = 100, precision = 10^-3, BB)
  # saveRDS(boot.b, "Boot.dirichlet.rds")
}




mse.Rd.yy <- boot.b[[1]]
mse.Rd.in <- boot.b[[2]]
mse.mu.in.b <- boot.b[[3]]
mse.yy.b <- boot.b[[4]]



df.mse <- data.frame(cbind(Prov = data$Prov, Age = data$Age,Sex = data$Sex, 
                           mse.yy.b, mse.mu.in.b, nd = data$nd))
colnames(df.mse) <- c("Prov","Age","Sex",
                      "DirectY1","DirectY2", "DirectY3", 
                      "Plug.inY1", "Plug.inY2", "Plug.inY3",
                      "nd")
