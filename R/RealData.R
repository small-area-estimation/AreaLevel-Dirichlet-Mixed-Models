#######################################################################################
###
###                           Application for paper
###                     Area Level Dirichlet Regression 
###                             
###                           
###                      Real data of EPA 2022
###                           
###                             October 2024                           
###
### AUTHOR: Esteban Cabello
### Date of creation: 24 September 2024
### Edited by: Esteban Cabello
### Edited on: 16 September 2025

source("functions_MLE_LaplaceT.R")
source("functions_predictors.R")


library(openxlsx)
library(nloptr)
library(DirichletReg)
library(dplyr)
library(xtable)

##### Lectura de datos
data <- read.xlsx("data.Dirichlet.ProvAge2Sex.xlsx")
## Unicamente nos quedamos con los grupos de poblacion 16-34,35-64. 
## Grupo de jubilados no interesa.
data <- data %>% filter(Age <=2)
data$Prov <- as.numeric(data$Prov)


var.rate <- data$Y2^2/((data$Y1+data$Y2)^4)*data$varY1 + 
            data$Y1^2/((data$Y1+data$Y2)^4)*data$varY2 - 
            2*(data$Y1*data$Y2)/((data$Y1+data$Y2)^4)*data$covY1Y2 

quantile.nd <- round(quantile(data$nd, probs = seq(0,1,0.1)),0);quantile.nd
fd <- round(quantile(100*data$nd/data$Nd, probs = seq(0,1,0.1)),3);fd
hist(data$nd, main = expression(n[d] ~ "distribution"), xlab = "")

CV.y1 <- 100*sqrt(data$varY1)/data$Y1
CV.y2 <- 100*sqrt(data$varY2)/data$Y2
CV.y3 <- 100*sqrt(data$varY3)/data$Y3
CV.rate <- 100*sqrt(var.rate)/(data$Y2/(data$Y1+data$Y2))

### CVs basados en el diseño. Desempleados tiene CVs altos pero debidos a la distorsión de la media (muy cercana a 0).
boxplot(CV.y1, CV.y2, CV.y3, CV.rate, main = "Design-based CVs", names = c("Employed", "Unemployed", "Inactive", "Rate"))
abline(h = 20,lty = 2, col = "red")

quantile.CV.y1 <- round(quantile(CV.y1, probs = seq(0,1,0.1)),0);quantile.CV.y1
quantile.CV.y2 <- round(quantile(CV.y2, probs = seq(0,1,0.1)),0);quantile.CV.y2
quantile.CV.y3 <- round(quantile(CV.y3, probs = seq(0,1,0.1)),0);quantile.CV.y3
quantile.CV.rate <- round(quantile(CV.rate, probs = seq(0,1,0.1)),0);quantile.CV.rate

boxplot(data$Y1,data$Y2,data$Y3, main = "Domain target variables", names = c("Employed", "Unemployed", "Inactive"))


boxplot(data$covY1Y2, data$covY1Y3, data$covY2Y3, main = "Design-based Covariances", 
        axes = F, cex.main = 1.8)
axis(2, col="grey", cex.axis = 1.6)  
axis(1, col="grey", cex.axis = 1.6, at = 1:3, labels = rep("",3)) 
mtext(text = c(expression(Y[d1]-Y[d2]), expression(Y[d1]-Y[d3]), expression(Y[d2]-Y[d3])), 
      side = 1, at = 1:3, cex = 1.8, line = 2)
abline(h = 0,lwd = 4, lty = 3, col = "red")

sd.edu3 <- sqrt(data$edu3.var)
sd.work1 <- sqrt(data$work1.var)
sd.edu1 <- sqrt(data$edu1.var)
sd.situ3 <- sqrt(data$situ3.var)

quantile.sd.edu3 <- round(quantile(sd.edu3, probs = seq(0,1,0.25)),4);quantile.sd.edu3
quantile.sd.work1 <- round(quantile(sd.work1, probs = seq(0,1,0.25)),4);quantile.sd.work1
quantile.sd.edu1 <- round(quantile(sd.edu1, probs = seq(0,1,0.25)),4);quantile.sd.edu1
quantile.sd.situ3 <- round(quantile(sd.situ3, probs = seq(0,1,0.25)),4);quantile.sd.situ3

sd.y1 <- sqrt(data$varY1)
sd.y2 <- sqrt(data$varY2)

quantile.sd.y1 <- round(quantile(sd.y1, probs = seq(0,1,0.25)),4);quantile.sd.y1
quantile.sd.y2 <- round(quantile(sd.y2, probs = seq(0,1,0.25)),4);quantile.sd.y2

xtable(rbind(quantile.sd.edu3, quantile.sd.work1, quantile.sd.edu1, quantile.sd.situ3,
             quantile.sd.y1, quantile.sd.y2), digits = 4)

qd1.1 <- data$varY1/data$edu3.var
qd1.2 <- data$varY1/data$work1.var
qd2.1 <- data$varY2/data$edu1.var
qd2.2 <- data$varY2/data$situ3.var

xtable(rbind(quantile(qd1.1), quantile(qd1.2), quantile(qd2.1), quantile(qd2.2)), digits = 4)





##### ADMM (Area-level Dirichlet Mixed Model)

q <- 3 
D <- nrow(data)

Intercept <- rep(1,D)

p1 <- 3    #  number of parameters beta1, beta2
p2 <- 3
p <- p1+p2
XX.k <- list(XX.1 = cbind(Intercept, data$edu3.MEAN, data$work1.MEAN), 
             XX.2 = cbind(Intercept, data$edu1.MEAN, data$situ3.MEAN))


### Dirichlet Reg for fixed effects (necessary to beta.s, sigma.s and seeds.u)


data$Z <- DR_data(data[,10:12], base = 3, trafo = TRUE) 
yy <- data$Z



phi.d <- (1/(q-1)*(data$Z[,1]*(1-data$Z[,1])/data$varY1+(data$Z[,2]*(1-data$Z[,2])/data$varY2))-1)
hist(phi.d, main = expression(phi[d]~"distribution"), xlab = "")

reg.d1 <- DirichletReg::DirichReg(formula = Z ~ edu3.MEAN +  work1.MEAN, data, model = "alternative", base = 3)
summary(reg.d1)
reg.d2 <- DirichletReg::DirichReg(formula = Z ~ edu1.MEAN +  situ3.MEAN, data, model = "alternative", base = 3)
summary(reg.d2)


beta.k <- rbind(c(reg.d1$coefficients[1:3]), c(reg.d2$coefficients[4:6]))
beta.s <- split(beta.k, row(beta.k))

sigma.s <- c(0.25,0.25) ### Vector of seeds
seeds.u <- split(rep(0, D*(q-1)), 1:D)


opts <- list("algorithm" = "NLOPT_LN_NEWUOA_BOUND", maxeval=50000)

dirichletM <- MLE_Laplace(XX.k, yy, phi.d, beta.s, sigma.s, seeds.u, maxiter=1000, epsilon.u=0.001, epsilon.par=0.001)
dirichletM$par.est
dirichletM$convergence



ud.mode <- dirichletM$udk.pred
ud.mod <- split(Reduce(ud.mode, f=cbind),1:(q-1))



####### Predicciones

pred.in <- plug.in(dirichletM$par.est, ud.mod, XX.k, yy,  phi.d)
mud.in <- Reduce(pred.in$mud.in, f=cbind)
head(mud.in)
Rd.in  <- pred.in$Rd.in
Rd.y <- yy[,2]/(yy[,2] + yy[,1])

### Residuos (Normalidad, QQ-plots y Boxplots)

sresiduals <- (yy - mud.in)/sqrt(mud.in*(1-mud.in)/(phi.d+1))
boxplot(sresiduals, main = "Boxplot of sresiduals", frame = F, cex.main = 1.8, cex.axis = 1.6, ylim = c(-2.1,2.1))
abline(h = 0, col = "red", lwd = 2, lty = 2)
abline(h = -2, col = "black", lwd = 2, lty = 2)
abline(h = 2, col = "black", lwd = 2, lty = 2)


hist(sresiduals[,1], freq = F,ylim = c(0,0.7), breaks = 10,cex.main = 1.8, cex.axis = 1.6, cex.lab = 1.6, main = "sresiduals - Employed", xlab = "")
lines(density(sresiduals[,1]), col = "red")
hist(sresiduals[,2], freq = F, ylim = c(0,0.8), breaks = 10, cex.main = 1.8, cex.axis = 1.6, cex.lab = 1.6, main = "sresiduals - Unemployed", xlab = "")
lines(density(sresiduals[,2]), col = "red")
hist(sresiduals[,3], freq = F, ylim = c(0,0.6), breaks = 10,cex.main = 1.8, cex.axis = 1.6, cex.lab = 1.6, main = "sresiduals - Inactive", xlab = "")
lines(density(sresiduals[,3]), col = "red")

boxplot(sresiduals[,1],sresiduals[,2],sresiduals[,3], 
        main = "sresiduals by Category", cex.main = 1.8, axes = F)
axis(2, col="black", cex.axis = 1.4)
axis(1, at = 1:3, line = 0, col="black", labels = rep("",3))
mtext(text = c("Employed", "Unemployed", "Inactive"), 
      side = 1, at = 1:3, cex = 1.6, line = 2)
abline(h = 0, lwd = 2, lty = 2, col = "red")

boxplot(sresiduals[c(T,F,F,F),1], sresiduals[c(F,T,F,F),1], sresiduals[c(F,F,T,F),1],sresiduals[c(F,F,F,T),1],
        sresiduals[c(T,F,F,F),1], sresiduals[c(F,T,F,F),2], sresiduals[c(F,F,T,F),2],sresiduals[c(F,F,F,T),2],
        sresiduals[c(T,F,F,F),1], sresiduals[c(F,T,F,F),3], sresiduals[c(F,F,T,F),3],sresiduals[c(F,F,F,T),3],
        main = "sresiduals by Category-Sex-Age", 
        cex.main = 1.8, axes = F, ylim = c(-2,2.2)
)
axis(2, col="black", cex.axis = 1.4)
axis(1, at = 1:12, line = 0, col="black", labels = rep("",12))
mtext(text = c(expression("E"[M1]), expression("E"[W1]),expression("E"[M2]), expression("E"[W2]),
               expression("U"[M1]), expression("U"[W1]),expression("U"[M2]), expression("U"[W2]),
               expression("I"[M1]), expression("I"[W1]),expression("I"[M2]), expression("I"[W2])),
      side = 1, at = 1:12, cex = 1.6, line = 2.5)
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


shapiro.test(sresiduals$Y1)
shapiro.test(sresiduals$Y2)
shapiro.test(sresiduals$Y3)

nortest::lillie.test(sresiduals$Y1)
nortest::lillie.test(sresiduals$Y2)
nortest::lillie.test(sresiduals$Y3)




df.pred <- data.frame(cbind(Prov = data$Prov, Age = data$Age,Sex = data$Sex, yy, mud.in, Rd.y, Rd.in, nd = data$nd))
colnames(df.pred) <- c("Prov","Age","Sex","DirectY1","DirectY2", "DirectY3", "Plug.inY1", "Plug.inY2", "Plug.inY3", 
                       "RateDirect", "Rate.in", "nd")


df.pred.men.age1 <- df.pred %>% filter(Age == 1 & Sex == 1) %>% arrange(nd)
df.pred.men.age2 <- df.pred %>% filter(Age == 2 & Sex == 1) %>% arrange(nd)
df.pred.women.age1 <- df.pred %>% filter(Age == 1 & Sex == 2) %>% arrange(nd)
df.pred.women.age2 <- df.pred %>% filter(Age == 2 & Sex == 2) %>% arrange(nd)


t.test(df.pred.men.age1$Rate.in, df.pred.women.age1$Rate.in, paired = T, conf.level = 0.95)
t.test(df.pred.men.age2$Rate.in, df.pred.women.age2$Rate.in, paired = T, conf.level = 0.95)

######
#### Plug-in vs Direct estimates (Men - Age1)
#####

u <- max(df.pred.men.age1$DirectY1, df.pred.men.age1$Plug.inY1)
l <- min(df.pred.men.age1$DirectY1, df.pred.men.age1$Plug.inY1)
nsel <- df.pred.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age1$DirectY1,type="b",col="red",ylim=c(round(l,1),u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age1$Plug.inY1, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Employed Men 16-34 years old")), cex.main=1.9)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, at = seq(0.3,0.7,0.1), col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.74, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)




u <- max(df.pred.men.age1$DirectY2, df.pred.men.age1$Plug.inY2)
l <- min(df.pred.men.age1$DirectY2, df.pred.men.age1$Plug.inY2)
nsel <- df.pred.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age1$DirectY2,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age1$Plug.inY2, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployed Men 16-34 years old")), cex.main=1.9)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.38, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.pred.men.age1$DirectY3, df.pred.men.age1$Plug.inY3)
l <- min(df.pred.men.age1$DirectY3, df.pred.men.age1$Plug.inY3)
nsel <- df.pred.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age1$DirectY3,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age1$Plug.inY3, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Inactive Men 16-34 years old")), cex.main=1.9)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.57, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



#### Rates of unemployment (Men - Age1)

u <- max(df.pred.men.age1$RateDirect, df.pred.men.age1$Rate.in)
l <- min(df.pred.men.age1$RateDirect, df.pred.men.age1$Rate.in)
nsel <- df.pred.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age1$RateDirect,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age1$Rate.in, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployment Rate Men 16-34 years old")), cex.main=1.9)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.55, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



#######
#### Plug-in vs Direct estimates (Men - Age2)
######

u <- max(df.pred.men.age2$DirectY1, df.pred.men.age2$Plug.inY1)
l <- min(df.pred.men.age2$DirectY1, df.pred.men.age2$Plug.inY1)
nsel <- df.pred.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age2$DirectY1,type="b",col="red",ylim=c(0.6,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age2$Plug.inY1, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Employed Men 35-64 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, at = seq(0.6,0.9,0.1), col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.94, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



u <- max(df.pred.men.age2$DirectY2, df.pred.men.age2$Plug.inY2)
l <- min(df.pred.men.age2$DirectY2, df.pred.men.age2$Plug.inY2)
nsel <- df.pred.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age2$DirectY2,type="b",col="red",ylim=c(0,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age2$Plug.inY2, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployed Men 35-64 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, at = seq(0,0.2,0.05), col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.22, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.pred.men.age2$DirectY3, df.pred.men.age2$Plug.inY3)
l <- min(df.pred.men.age2$DirectY3, df.pred.men.age2$Plug.inY3)
nsel <- df.pred.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age2$DirectY3,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age2$Plug.inY3, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Inactive Men 35-64 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.28, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


#######
#### Rates of unemployment (Men - Age2)
######
u <- max(df.pred.men.age2$RateDirect, df.pred.men.age2$Rate.in)
l <- min(df.pred.men.age2$RateDirect, df.pred.men.age2$Rate.in)
nsel <- df.pred.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.men.age2$RateDirect,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.men.age2$Rate.in, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployment Rate Men 35-64 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.24, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



#### Plug-in vs Direct estimates (Women - Age1)

u <- max(df.pred.women.age1$DirectY1, df.pred.women.age1$Plug.inY1)
l <- min(df.pred.women.age1$DirectY1, df.pred.women.age1$Plug.inY1)
nsel <- df.pred.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age1$DirectY1,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age1$Plug.inY1, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Employed Women 16-34 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.64, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



u <- max(df.pred.women.age1$DirectY2, df.pred.women.age1$Plug.inY2)
l <- min(df.pred.women.age1$DirectY2, df.pred.women.age1$Plug.inY2)
nsel <- df.pred.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age1$DirectY2,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age1$Plug.inY2, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployed Women 16-34 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.35, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.pred.women.age1$DirectY3, df.pred.women.age1$Plug.inY3)
l <- min(df.pred.women.age1$DirectY3, df.pred.women.age1$Plug.inY3)
nsel <- df.pred.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age1$DirectY3,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age1$Plug.inY3, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Inactive Women 16-34 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.55, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



#### Rates of unemployment (Women - Age1)

u <- max(df.pred.women.age1$RateDirect, df.pred.women.age1$Rate.in)
l <- min(df.pred.women.age1$RateDirect, df.pred.women.age1$Rate.in)
nsel <- df.pred.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age1$RateDirect,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age1$Rate.in, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployment  Rate Women 16-34 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.52, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



#### Plug-in vs Direct estimates (Women - Age2)

u <- max(df.pred.women.age2$DirectY1, df.pred.women.age2$Plug.inY1)
l <- min(df.pred.women.age2$DirectY1, df.pred.women.age2$Plug.inY1)
nsel <- df.pred.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age2$DirectY1,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age2$Plug.inY1, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Employed Women 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.88, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



u <- max(df.pred.women.age2$DirectY2, df.pred.women.age2$Plug.inY2)
l <- min(df.pred.women.age2$DirectY2, df.pred.women.age2$Plug.inY2)
nsel <- df.pred.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age2$DirectY2,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age2$Plug.inY2, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployed Women 35-64 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.28, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.pred.women.age2$DirectY3, df.pred.women.age2$Plug.inY3)
l <- min(df.pred.women.age2$DirectY3, df.pred.women.age2$Plug.inY3)
nsel <- df.pred.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age2$DirectY3,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age2$Plug.inY3, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Inactive Women 35-64 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.42, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



#### Rates of unemployment (Women - Age2)

u <- max(df.pred.women.age2$RateDirect, df.pred.women.age2$Rate.in)
l <- min(df.pred.women.age2$RateDirect, df.pred.women.age2$Rate.in)
nsel <- df.pred.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.pred.women.age2$RateDirect,type="b",col="red",ylim=c(l,u+0.05),xlim = c(0,55),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.pred.women.age2$Rate.in, type="b", col = "blue", pch=16)
# abline(h=0, lty = 2, lwd = 1, col = "gray70")
title(main=expression(paste("Unemployment Rate Women 35-64 years old")), cex.main=1.8)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(40, 0.37, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


##########################################
############### Bootstrap ################
##########################################

cond <- F
BB <- 500

if(cond == T){ # Computed in Castle-black
  
  source("boot_Dirichlet_app.R")
  iter <- 1
  library(mvtnorm)
  
  set.seed(03201)
  boot.b <- Boot.Dirichlet(XX.k, param.e = dirichletM$par.est,  MAXITER = 100, precision = 10^-3)
  saveRDS(boot.b, "Boot.dirichlet.age2.rds")
}


boot.b <- readRDS("Boot.dirichlet.age2.rds")


mse.Rd.yy <- boot.b[[1]]
mse.Rd.in <- boot.b[[2]]
mse.mu.in.b <- boot.b[[3]]
mse.yy.b <- boot.b[[4]]



df.mse <- data.frame(cbind(Prov = data$Prov, Age = data$Age,Sex = data$Sex, 
                           mse.yy.b, mse.mu.in.b,nd = data$nd))
colnames(df.mse) <- c("Prov","Age","Sex",
                      "DirectY1","DirectY2", "DirectY3", 
                      "Plug.inY1", "Plug.inY2", "Plug.inY3",
                      "nd")

df.mse.men.age1 <- df.mse %>% filter(Age == 1 & Sex == 1) %>% arrange(nd)
df.mse.men.age2 <- df.mse %>% filter(Age == 2 & Sex == 1) %>% arrange(nd)
df.mse.women.age1 <- df.mse %>% filter(Age == 1 & Sex == 2) %>% arrange(nd)
df.mse.women.age2 <- df.mse %>% filter(Age == 2 & Sex == 2) %>% arrange(nd)


#### Mean reduction MSE by Category-Age-Sex

### Men 16-34
mean(100*(df.mse.men.age1$DirectY1 - df.mse.men.age1$Plug.inY1)/df.mse.men.age1$DirectY1)
mean(100*(df.mse.men.age1$DirectY2 - df.mse.men.age1$Plug.inY2)/df.mse.men.age1$DirectY2)
mean(100*(df.mse.men.age1$DirectY3 - df.mse.men.age1$Plug.inY3)/df.mse.men.age1$DirectY3)

### Men 35-64
mean(100*(df.mse.men.age2$DirectY1 - df.mse.men.age2$Plug.inY1)/df.mse.men.age2$DirectY1)
mean(100*(df.mse.men.age2$DirectY2 - df.mse.men.age2$Plug.inY2)/df.mse.men.age2$DirectY2)
mean(100*(df.mse.men.age2$DirectY3 - df.mse.men.age2$Plug.inY3)/df.mse.men.age2$DirectY3)

### Women 16-34
mean(100*(df.mse.women.age1$DirectY1 - df.mse.women.age1$Plug.inY1)/df.mse.women.age1$DirectY1)
mean(100*(df.mse.women.age1$DirectY2 - df.mse.women.age1$Plug.inY2)/df.mse.women.age1$DirectY2)
mean(100*(df.mse.women.age1$DirectY3 - df.mse.women.age1$Plug.inY3)/df.mse.women.age1$DirectY3)

### Women 35-64
mean(100*(df.mse.women.age2$DirectY1 - df.mse.women.age2$Plug.inY1)/df.mse.women.age2$DirectY1)
mean(100*(df.mse.women.age2$DirectY2 - df.mse.women.age2$Plug.inY2)/df.mse.women.age2$DirectY2)
mean(100*(df.mse.women.age2$DirectY3 - df.mse.women.age2$Plug.inY3)/df.mse.women.age2$DirectY3)


###########
### mse* Men Age 1
#########

u <- max(df.mse.men.age1$DirectY1, df.mse.men.age1$Plug.inY1)
l <- min(df.mse.men.age1$DirectY1, df.mse.men.age1$Plug.inY1)
nsel <- df.mse.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.men.age1$DirectY1, type="b",col="red", ylim=c(l,u), xlim = c(0,52), xlab="", ylab="", xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.men.age1$Plug.inY1, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Employed Men 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.men.age1$DirectY2, df.mse.men.age1$Plug.inY2)
l <- min(df.mse.men.age1$DirectY2, df.mse.men.age1$Plug.inY2)
nsel <- df.mse.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.men.age1$DirectY2,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.men.age1$Plug.inY2, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployed Men 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.men.age1$DirectY3, df.mse.men.age1$Plug.inY3)
l <- min(df.mse.men.age1$DirectY3, df.mse.men.age1$Plug.inY3)
nsel <- df.mse.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.men.age1$DirectY3,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.men.age1$Plug.inY3, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Inactive Men 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2, col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



###########
### mse* Men Age 2
#########

u <- max(df.mse.men.age2$DirectY1, df.mse.men.age2$Plug.inY1)
l <- min(df.mse.men.age2$DirectY1, df.mse.men.age2$Plug.inY1)
nsel <- df.mse.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.men.age2$DirectY1, type="b",col="red", ylim=c(l,u), xlim = c(0,52), xlab="", ylab="", xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.men.age2$Plug.inY1, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Employed Men 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.men.age2$DirectY2, df.mse.men.age2$Plug.inY2)
l <- min(df.mse.men.age2$DirectY2, df.mse.men.age2$Plug.inY2)
nsel <- df.mse.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.men.age2$DirectY2,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.men.age2$Plug.inY2, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployed Men 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.men.age2$DirectY3, df.mse.men.age2$Plug.inY3)
l <- min(df.mse.men.age2$DirectY3, df.mse.men.age2$Plug.inY3)
nsel <- df.mse.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.men.age2$DirectY3,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.men.age2$Plug.inY3, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Inactive Men 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


###########
### mse* Women Age 1
#########

u <- max(df.mse.women.age1$DirectY1, df.mse.women.age1$Plug.inY1)
l <- min(df.mse.women.age1$DirectY1, df.mse.women.age1$Plug.inY1)
nsel <- df.mse.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.women.age1$DirectY1,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.women.age1$Plug.inY1, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Employed Women 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.women.age1$DirectY2, df.mse.women.age1$Plug.inY2)
l <- min(df.mse.women.age1$DirectY2, df.mse.women.age1$Plug.inY2)
nsel <- df.mse.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.women.age1$DirectY2,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.women.age1$Plug.inY2, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployed Women 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.women.age1$DirectY3, df.mse.women.age1$Plug.inY3)
l <- min(df.mse.women.age1$DirectY3, df.mse.women.age1$Plug.inY3)
nsel <- df.mse.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.women.age1$DirectY3,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.women.age1$Plug.inY3, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Inactive Women 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


###########
### mse* Women Age 2
#########

u <- max(df.mse.women.age2$DirectY1, df.mse.women.age2$Plug.inY1)
l <- min(df.mse.women.age2$DirectY1, df.mse.women.age2$Plug.inY1)
nsel <- df.mse.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.women.age2$DirectY1,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.women.age2$Plug.inY1, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Employed Women 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.women.age2$DirectY2, df.mse.women.age2$Plug.inY2)
l <- min(df.mse.women.age2$DirectY2, df.mse.women.age2$Plug.inY2)
nsel <- df.mse.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.women.age2$DirectY2,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.women.age2$Plug.inY2, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployed Women 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.women.age2$DirectY3, df.mse.women.age2$Plug.inY3)
l <- min(df.mse.women.age2$DirectY3, df.mse.women.age2$Plug.inY3)
nsel <- df.mse.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.women.age2$DirectY3,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.women.age2$Plug.inY3, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Inactive Women 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)




##########
### mse* Rd 
#########

df.mse.Rd <- data.frame(cbind(Prov = data$Prov, Age = data$Age,Sex = data$Sex, mse.Rd.yy, mse.Rd.in, nd = data$nd))
colnames(df.mse.Rd) <- c("Prov","Age","Sex","Rd.yy", "Rd.in", "nd")

df.mse.Rd.men.age1 <- df.mse.Rd %>% filter(Age == 1 & Sex == 1) %>% arrange(nd)
df.mse.Rd.men.age2 <- df.mse.Rd %>% filter(Age == 2 & Sex == 1) %>% arrange(nd)
df.mse.Rd.women.age1 <- df.mse.Rd %>% filter(Age == 1 & Sex == 2) %>% arrange(nd)
df.mse.Rd.women.age2 <- df.mse.Rd %>% filter(Age == 2 & Sex == 2) %>% arrange(nd)


#### Mean reduction MSE by Category-Age-Sex

### Men 16-34
mean(100*(df.mse.Rd.men.age1$Rd.yy - df.mse.Rd.men.age1$Rd.in)/df.mse.Rd.men.age1$Rd.yy)

### Men 35-64
mean(100*(df.mse.Rd.men.age2$Rd.yy - df.mse.Rd.men.age2$Rd.in)/df.mse.Rd.men.age2$Rd.yy)

### Women 16-34
mean(100*(df.mse.Rd.women.age1$Rd.yy - df.mse.Rd.women.age1$Rd.in)/df.mse.Rd.women.age1$Rd.yy)

### Women 35-64
mean(100*(df.mse.Rd.women.age2$Rd.yy - df.mse.Rd.women.age2$Rd.in)/df.mse.Rd.women.age2$Rd.yy)


u <- max(df.mse.Rd.men.age1$Rd.yy, df.mse.Rd.men.age1$Rd.in)
l <- min(df.mse.Rd.men.age1$Rd.yy, df.mse.Rd.men.age1$Rd.in)
nsel <- df.mse.Rd.men.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.Rd.men.age1$Rd.yy,type="b",col="red",ylim=c(0,u+0.001),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.Rd.men.age1$Rd.in, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployment Rate Men 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2, line=4)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.4) #Axis OX
axis(2,col = "gray70", cex.axis=1.5)
text <- c("Direct", "Plug-in")
legend(35, u+0.001, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)

u <- max(df.mse.Rd.men.age2$Rd.yy, df.mse.Rd.men.age2$Rd.in)
l <- min(df.mse.Rd.men.age2$Rd.yy, df.mse.Rd.men.age2$Rd.in)
nsel <- df.mse.Rd.men.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.Rd.men.age2$Rd.yy,type="b",col="red",ylim=c(0,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.Rd.men.age2$Rd.in, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployment Rate Men 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



u <- max(df.mse.Rd.women.age1$Rd.yy, df.mse.Rd.women.age1$Rd.in)
l <- min(df.mse.Rd.women.age1$Rd.yy, df.mse.Rd.women.age1$Rd.in)
nsel <- df.mse.Rd.women.age1$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.Rd.women.age1$Rd.yy,type="b",col="red",ylim=c(0,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.Rd.women.age1$Rd.in, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployment Rate Women 16-34 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)


u <- max(df.mse.Rd.women.age2$Rd.yy, df.mse.Rd.women.age2$Rd.in)
l <- min(df.mse.Rd.women.age2$Rd.yy, df.mse.Rd.women.age2$Rd.in)
nsel <- df.mse.Rd.women.age2$nd[round(seq(1, 52, length=13),0)]

plot(df.mse.Rd.women.age2$Rd.yy,type="b",col="red",ylim=c(l,u),xlim = c(0,52),xlab="",ylab="",xaxt = "n", pch=15, cex.axis=1.5, axes = F)
lines(df.mse.Rd.women.age2$Rd.in, type="b", col = "blue", pch=16)
title(main=expression(paste("mse* - Unemployment Rate Women 35-64 years old")), cex.main=2)
title(xlab = "Sample size", cex.lab=2.5, line=3.5)
axis(side=1, at=round(seq(1,52,length=13),0), labels=nsel,line=0, col="gray70",cex.axis=1.6) #Axis OX
axis(2, col = "gray70", cex.axis=1.6)
text <- c("Direct", "Plug-in")
legend(35, u, text, col = c("red","blue"), pch = c(15, 16), bty="n", lwd = 2, cex=1.5)



##### CI and p-values Bootstrap for model parameters


param.b <- boot.b[[7]]  


alpha2 <- 0.05/2

percent.theta <- apply(param.b, 1, function(x) quantile(x, probs = c(alpha2,1-alpha2)))
sd.theta.hat <- apply(param.b, 1, sd)

CI.low.b <- dirichletM$par.est - sd.theta.hat*1.96
CI.upp.b <- dirichletM$par.est + sd.theta.hat*1.96
CI.b <- rbind(CI.low.b, CI.upp.b)
colnames(CI.b) <- c("Intercept", "edu1", "nac1", "Intercept", "edu1", "situ3", "sigma1", "sigma2")



percentb.theta <- apply(sqrt(BB)*(param.b-dirichletM$par.est), 1, function(x) quantile(x, probs = c(alpha2,1-alpha2)))

CI.low.basic <- dirichletM$par.est - percentb.theta[2,]/sqrt(BB)
CI.upp.basic <- dirichletM$par.est - percentb.theta[1,]/sqrt(BB)
CI.basic <- rbind(CI.low.basic, CI.upp.basic)
colnames(CI.basic) <- c("Intercept", "edu1", "nac1", "Intercept", "edu1", "situ3", "sigma1", "sigma2")

xtable(CI.b, digits = 4)
xtable(CI.basic, digits = 4)





#####################################################################################################
################### MAPEO PROPORCION DESEMPLEADOS Y TASAS DESEMPLEO (2022.T4) #######################
#####################################################################################################

library(mapSpain)
library(ggplot2)


###############
###############
##### Men 16-34
###############
###############

prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.in.men.age1 <- (df.pred.men.age1 %>% arrange(Prov))$Rate.in
quantile(rate.in.men.age1)


breaks.m <- cut(rate.in.men.age1, breaks=c(0.05,0.10,0.15,0.20,0.25,0.3,Inf))
breaks.pred.m <- data.frame(unique(prov$Prov), breaks.m)
names(breaks.pred.m) <- c('ine.prov.name', 'breaks.m')

prov <- merge(prov, breaks.pred.m, by='ine.prov.name')

legend_labs <- c("0.05-0.10", "0.10-0.15", "0.15-0.20", 
                 "0.20-0.25", "0.25-0.30", "> 0.30")
legend_cols <- hcl.colors(7, "Oslo")   

ggplot(prov) +
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.m), color = "grey70") +
  labs(title = "Plug-in Unemployment Rate Men 16-34 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.m),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )




prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.dir.men.age1 <- (df.pred.men.age1 %>% arrange(Prov))$RateDirect
quantile(rate.dir.men.age1)


breaks.m <- cut(rate.dir.men.age1, breaks=c(0.05,0.10,0.15,0.20,0.25,0.3,Inf))
breaks.pred.m <- data.frame(unique(prov$Prov), breaks.m)
names(breaks.pred.m) <- c('ine.prov.name', 'breaks.m')

prov <- merge(prov, breaks.pred.m, by='ine.prov.name')

legend_labs <- c("0.05-0.10", "0.10-0.15", "0.15-0.20", 
                 "0.20-0.25", "0.25-0.30", "> 0.30")
legend_cols <- hcl.colors(7, "Oslo")   

ggplot(prov) +
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.m), color = "grey70") +
  labs(title = "Direct Unemployment Rate Men 16-34 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.m),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )


###############
###############
##### Women 16-34
###############
###############


prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.dir.women.age1 <- (df.pred.women.age1 %>% arrange(Prov))$RateDirect
quantile(rate.dir.women.age1)


breaks.w <- cut(rate.dir.women.age1, breaks=c(0.05,0.10,0.15,0.20,0.25,0.3,Inf))
breaks.pred.w <- data.frame(unique(prov$Prov), breaks.w)
names(breaks.pred.w) <- c('ine.prov.name', 'breaks.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

legend_labs <- c("0.05-0.10", "0.10-0.15", "0.15-0.20", 
                 "0.20-0.25", "0.25-0.30", "> 0.30")
legend_cols <- hcl.colors(7, "Oslo")   

ggplot(prov) +
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.w), color = "grey70") +
  labs(title = "Direct Unemployment Rate Women 16-34 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.w),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )


prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.in.women.age1 <- (df.pred.women.age1 %>% arrange(Prov))$Rate.in
quantile(rate.in.women.age1)


breaks.w <- cut(rate.in.women.age1, breaks=c(0.05,0.10,0.15,0.20,0.25,0.3,Inf))
breaks.pred.w <- data.frame(unique(prov$Prov), breaks.w)
names(breaks.pred.w) <- c('ine.prov.name', 'breaks.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

legend_labs <- c("0.05-0.10", "0.10-0.15", "0.15-0.20", 
                 "0.20-0.25", "0.25-0.30", "> 0.30")
legend_cols <- hcl.colors(7, "Oslo")   

ggplot(prov) +
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.w), color = "grey70") +
  labs(title = "Plug-in Unemployment Rate Women 16-34 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.w),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )




###############
###############
##### Men 35-64
###############
###############

prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.in.men.age2 <- (df.pred.men.age2 %>% arrange(Prov))$Rate.in
quantile(rate.in.men.age2)


breaks.m <- cut(rate.in.men.age2, breaks=seq(0.0, 0.2, by= 0.05))
breaks.pred.m <- data.frame(unique(prov$Prov), breaks.m)
names(breaks.pred.m) <- c('ine.prov.name', 'breaks.m')

prov <- merge(prov, breaks.pred.m, by='ine.prov.name')

legend_labs <- c("<0.05", "0.05-0.10", "0.10-0.15", "0.15-0.20")
legend_cols <- hcl.colors(7, "Oslo")   

ggplot(prov) +    
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.m), color = "grey70") +
  labs(title = "Plug-in Unemployment Rate Men 35-64 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.m),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )


prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.dir.men.age2 <- (df.pred.men.age2 %>% arrange(Prov))$RateDirect
quantile(rate.dir.men.age2)


breaks.m <- cut(rate.dir.men.age2, breaks=seq(0.0, 0.2, by= 0.05))
breaks.pred.m <- data.frame(unique(prov$Prov), breaks.m)
names(breaks.pred.m) <- c('ine.prov.name', 'breaks.m')

prov <- merge(prov, breaks.pred.m, by='ine.prov.name')

legend_labs <- c("<0.05", "0.05-0.10", "0.10-0.15", "0.15-0.20")
legend_cols <- hcl.colors(7, "Oslo")   

ggplot(prov) +    
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.m), color = "grey70") +
  labs(title = "Direct Unemployment Rate Men 35-64 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.m),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )


###############
###############
##### Women 35-64
###############
###############


prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.in.women.age2 <- (df.pred.women.age2 %>% arrange(Prov))$Rate.in
quantile(rate.in.women.age2)


breaks.w <- cut(rate.in.women.age2, breaks=seq(0.05, 0.35, by= 0.05))
breaks.pred.w <- data.frame(unique(prov$Prov), breaks.w)
names(breaks.pred.w) <- c('ine.prov.name', 'breaks.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

legend_labs <- c("0.05-0.10", "0.10-0.15", "0.15-0.20", "0.20-0.25",
                 "0.25-0.30", ">0.30")
legend_cols <- hcl.colors(7, "Oslo")[2:7]   


ggplot(prov) +    
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.w), color = "grey70") +
  labs(title = "Plug-in Unemployment Rate Women 35-64 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.w),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )



prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

rate.dir.women.age2 <- (df.pred.women.age2 %>% arrange(Prov))$RateDirect
quantile(rate.dir.women.age2)


breaks.w <- cut(rate.dir.women.age2, breaks=seq(0.05, 0.35, by= 0.05))
breaks.pred.w <- data.frame(unique(prov$Prov), breaks.w)
names(breaks.pred.w) <- c('ine.prov.name', 'breaks.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

legend_labs <- c("0.05-0.10", "0.10-0.15", "0.15-0.20", "0.20-0.25",
                 "0.25-0.30", ">0.30")
legend_cols <- hcl.colors(7, "Oslo")[2:7]   

ggplot(prov) +    
  geom_sf(data = esp_get_country(),    color = "grey80", fill = NA) +
  geom_sf(data = esp_get_can_box(),    color = "grey70", fill = NA) +
  geom_sf(data = esp_get_can_provinces(), color = "grey70", fill = NA) +
  geom_sf(aes(fill = breaks.w), color = "grey70") +
  labs(title = "Direct Unemployment Rate Women 35-64 years old") +
  scale_fill_manual(
    name   = "",                          # legend title
    values = legend_cols,                 # colors
    breaks = levels(prov$breaks.w),       # order of categories in the legend
    labels = legend_labs,                 # your labels
    drop   = FALSE
  ) +
  theme_classic() +
  theme(
    text        = element_text(size = 18),
    legend.position = c(0.9, 0.23),
    plot.title  = element_text(hjust = 0.5)
  )




###############################################################################
################### MAPEO CV TASAS DESEMPLEO (2022.T4) #######################
##############################################################################


###################
##### Plug-in Men (Age1)
################### 


CVrate.men.age1 <- data.frame(Prov = df.pred.men.age1$Prov, 
                              CV.in = 100*sqrt(df.mse.Rd.men.age1$Rd.in)/(df.pred.men.age1$Rate.in),
                              CV.dir = 100*sqrt(df.mse.Rd.men.age1$Rd.yy)/df.pred.men.age1$Rate.in)

CVrate.men.age2 <- data.frame(Prov = df.pred.men.age2$Prov, 
                              CV.in = 100*sqrt(df.mse.Rd.men.age2$Rd.in)/(df.pred.men.age2$Rate.in),
                              CV.dir = 100*sqrt(df.mse.Rd.men.age2$Rd.yy)/df.pred.men.age2$Rate.in)

CVrate.women.age1 <- data.frame(Prov = df.pred.women.age1$Prov, 
                                CV.in = 100*sqrt(df.mse.Rd.women.age1$Rd.in)/(df.pred.women.age1$Rate.in),
                                CV.dir = 100*sqrt(df.mse.Rd.women.age1$Rd.yy)/df.pred.women.age1$Rate.in)

CVrate.women.age2 <- data.frame(Prov = df.pred.women.age2$Prov, 
                                CV.in = 100*sqrt(df.mse.Rd.women.age2$Rd.in)/(df.pred.women.age2$Rate.in),
                                CV.dir = 100*sqrt(df.mse.Rd.women.age2$Rd.yy)/df.pred.women.age2$Rate.in)


prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

CVrate.men.age1 <- merge(CVrate.men.age1, 
                         data.frame(Prov = as.numeric(prov$cpro), ine.prov.name = prov$ine.prov.name), by = "Prov")
CVrate.men.age2 <- merge(CVrate.men.age2, 
                         data.frame(Prov = as.numeric(prov$cpro), ine.prov.name = prov$ine.prov.name), by = "Prov")
CVrate.women.age1 <- merge(CVrate.women.age1, 
                           data.frame(Prov = as.numeric(prov$cpro), ine.prov.name = prov$ine.prov.name), by = "Prov")
CVrate.women.age2 <- merge(CVrate.women.age2, 
                           data.frame(Prov = as.numeric(prov$cpro), ine.prov.name = prov$ine.prov.name), by = "Prov")

prov$cpro <- as.numeric(prov$cpro)

CVrate.in.men.age1 <- (CVrate.men.age1 %>% arrange(Prov))$CV.in
quantile(CVrate.in.men.age1)

breaks.in.m <- cut(CVrate.in.men.age1, breaks=seq(10, 30, by= 5))
breaks.pred.in.m <- data.frame(unique(prov$Prov), breaks.in.m)
names(breaks.pred.in.m) <- c('ine.prov.name', 'breaks.in.m')

prov <- merge(prov, breaks.pred.in.m, by='ine.prov.name')


ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.in.m), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{"in"}, " Unemployment Rate Men 16-34 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[7:4]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(0.9,.23),
                          plot.title = element_text(hjust = 0.5) )



CVrate.in.men.age2 <- (CVrate.men.age2 %>% arrange(Prov))$CV.in
quantile(CVrate.in.men.age2)

breaks.in.m <- cut(CVrate.in.men.age2, breaks=seq(5, 30, by= 5))
breaks.pred.in.m <- data.frame(unique(prov$Prov), breaks.in.m)
names(breaks.pred.in.m) <- c('ine.prov.name', 'breaks.in.m')

prov <- merge(prov, breaks.pred.in.m, by='ine.prov.name')


ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.in.m), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{"in"}, " Unemployment Rate Men 35-64 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[8:4]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(0.9,.23), 
                          plot.title = element_text(hjust = 0.5) )

## Sorted by CV.in

CVrate.men.age1 %>% arrange(CV.in)
CVrate.men.age2 %>% arrange(CV.in)

CVrate.women.age1 %>% arrange(CV.in)
CVrate.women.age2 %>% arrange(CV.in)



###################
##### Direct Men 
###################


prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

CVrate.dir.men.age1 <- (CVrate.men.age1 %>% arrange(Prov))$CV.dir
quantile(CVrate.dir.men.age1)


breaks.dir.m <- cut(CVrate.dir.men.age1, breaks=seq(10, 40, by= 5))
breaks.pred.dir.m <- data.frame(unique(prov$Prov), breaks.dir.m)
names(breaks.pred.dir.m) <- c('ine.prov.name', 'breaks.dir.m')

prov <- merge(prov, breaks.pred.dir.m, by='ine.prov.name')

ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.dir.m), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{dir}, " Unemployment Rate Men 16-34 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[7:2]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(0.9,.23), 
                          plot.title = element_text(hjust = 0.5) )




CVrate.dir.men.age2 <- (CVrate.men.age2 %>% arrange(Prov))$CV.dir
quantile(CVrate.dir.men.age2)


breaks.dir.m <- cut(CVrate.dir.men.age2, breaks=seq(5, 45, by= 5))
breaks.pred.dir.m <- data.frame(unique(prov$Prov), breaks.dir.m)
names(breaks.pred.dir.m) <- c('ine.prov.name', 'breaks.dir.m')

prov <- merge(prov, breaks.pred.dir.m, by='ine.prov.name')

ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.dir.m), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{dir}, " Unemployment Rate Men 35-64 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[8:1]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(0.9,.23), 
                          plot.title = element_text(hjust = 0.35) )



###################
##################
##### Plug-in Women
###################
##################

prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

CVrate.in.women.age1 <- (CVrate.women.age1 %>% arrange(Prov))$CV.in
quantile(CVrate.in.women.age1)


breaks.in.w <- cut(CVrate.in.women.age1, breaks=seq(10, 25, by= 5))
breaks.pred.in.w <- data.frame(unique(prov$Prov), breaks.w)
names(breaks.pred.in.w) <- c('ine.prov.name', 'breaks.in.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.in.w), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{"in"}, " Unemployment Rate Women 16-34 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[7:5]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(.9,.23), 
                          plot.title = element_text(hjust = 0.5) )


CVrate.in.women.age2 <- (CVrate.women.age2 %>% arrange(Prov))$CV.in
quantile(CVrate.in.women.age2)


breaks.in.w <- cut(CVrate.in.women.age2, breaks=seq(5, 30, by= 5))
breaks.pred.in.w <- data.frame(unique(prov$Prov), breaks.w)
names(breaks.pred.in.w) <- c('ine.prov.name', 'breaks.in.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.in.w), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{"in"}, " Unemployment Rate Women 35-64 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[8:4]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(.9,.23), 
                          plot.title = element_text(hjust = 0.5) )

###################
##################
##### Direct Women
###################
##################

prov <- esp_get_prov() %>% 
  mutate(Provincia = ine.prov.name) %>% arrange(cpro)

prov$cpro <- as.numeric(prov$cpro)

CVrate.dir.women.age1 <- (CVrate.women.age1 %>% arrange(Prov))$CV.dir
quantile(CVrate.dir.women.age1)


breaks.dir.w <- cut(CVrate.dir.women.age1, breaks=seq(10, 40, by= 5))
breaks.pred.dir.w <- data.frame(unique(prov$Prov), breaks.dir.w)
names(breaks.pred.dir.w) <- c('ine.prov.name', 'breaks.dir.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.dir.w), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{dir}, " Unemployment Rate Women 16-34 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[7:2]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(.9,.23), 
                          plot.title = element_text(hjust = 0.5) )


CVrate.dir.women.age2 <- (CVrate.women.age2 %>% arrange(Prov))$CV.dir
quantile(CVrate.dir.women.age2)


breaks.dir.w <- cut(CVrate.dir.women.age2, breaks=seq(5, 40, by= 5))
breaks.pred.dir.w <- data.frame(unique(prov$Prov), breaks.dir.w)
names(breaks.pred.dir.w) <- c('ine.prov.name', 'breaks.dir.w')

prov <- merge(prov, breaks.pred.w, by='ine.prov.name')

ggplot(prov) +
  geom_sf(data = esp_get_country(), color = "grey80") +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill =breaks.dir.w), color = "grey70") +
  labs(title = expression(paste(widehat(CV)^{dir}, " Unemployment Rate Women 35-64 years old"))) +
  scale_fill_discrete('', type = hcl.colors(8, "YlOrRd")[8:2]) +
  theme_classic() + theme(text = element_text(size=18), 
                          legend.position=c(.9,.23), 
                          plot.title = element_text(hjust = 0.5) )


CVrate.men.age1 %>% filter(CV.dir >=30) %>% arrange(CV.dir)
CVrate.men.age2 %>% filter(CV.dir >=30) %>% arrange(CV.dir)

CVrate.women.age1 %>% filter(CV.dir >=30) %>% arrange(CV.dir)
CVrate.women.age2 %>% filter(CV.dir >=30) %>% arrange(CV.dir)




##### Tables

nd.men1 <- data %>% filter(Age == 1 & Sex == 1) %>% select(nd)
nd.men2 <- data %>% filter(Age == 2 & Sex == 1) %>% select(nd)
nd.women1 <- data %>% filter(Age == 1 & Sex == 2) %>% select(nd)
nd.women2 <- data %>% filter(Age == 2 & Sex == 2) %>% select(nd)

table.pred.CV.men1 <- data.frame(CVrate.men.age1,
                                 CV.design = CV.rate[c(T,F,F,F)], 
                                 nd = nd.men1,
                                 Plug.in = rate.in.men.age1,
                                 Direct = Rd.y[c(T,F,F,F)]) %>% select(ine.prov.name, nd, Direct, Plug.in, CV.design, CV.dir, CV.in) %>%
                      
                      arrange(nd) %>% slice(seq(1, n(), by = 5))


names.prov <- table.pred.CV.men1$ine.prov.name

table.pred.CV.men2 <- data.frame(CVrate.men.age2,
                                 CV.design = CV.rate[c(F,F,T,F)], 
                                 nd = nd.men2,
                                 Plug.in = rate.in.men.age2,
                                 Direct = Rd.y[c(F,F,T,F)]) %>% select(ine.prov.name, nd, Direct, Plug.in, CV.design, CV.dir, CV.in) %>%
  
  filter(ine.prov.name %in% names.prov) %>% arrange(nd)

table.pred.CV.women1 <- data.frame(CVrate.women.age1,
                                 CV.design = CV.rate[c(F,T,F,F)], 
                                 nd = nd.women1,
                                 Plug.in = rate.in.women.age1,
                                 Direct = Rd.y[c(F,T,F,F)]) %>% select(ine.prov.name, nd, Direct, Plug.in, CV.design, CV.dir, CV.in) %>%
  
  filter(ine.prov.name %in% names.prov) %>% arrange(nd)


table.pred.CV.women2 <- data.frame(CVrate.women.age2,
                                 CV.design = CV.rate[c(F,F,F,T)], 
                                 nd = nd.women2,
                                 Plug.in = rate.in.women.age2,
                                 Direct = Rd.y[c(F,F,F,T)]) %>% select(ine.prov.name, nd, Direct, Plug.in, CV.design, CV.dir, CV.in) %>%
  
  filter(ine.prov.name %in% names.prov) %>% arrange(nd)


print(xtable(table.pred.CV.men1, digits = 3), include.rownames = F)

print(xtable(cbind(table.pred.CV.men1, table.pred.CV.women1)[,-8], digits = 2), include.rownames = F)
print(xtable(table.pred.CV.men2, digits = 3), include.rownames = F)
print(xtable(table.pred.CV.women2, digits = 3), include.rownames = F)



par(mfrow = c(1,3))
boxplot(df.mse.Rd.men.age1$Rd.in, df.mse.Rd.men.age2$Rd.in, 
        df.mse.Rd.women.age1$Rd.in, df.mse.Rd.women.age2$Rd.in, main = "mse*", names = c("M1", "M2", "W1", "W2"))
boxplot(df.mse.Rd.men.age1$nd, df.mse.Rd.men.age2$nd, 
        df.mse.Rd.women.age1$nd,df.mse.Rd.women.age2$nd, main = "nd", names = c("M1", "M2", "W1", "W2"))
boxplot(CVrate.men.age1$CV.in, CVrate.men.age2$CV.in, 
        CVrate.women.age1$CV.in, CVrate.women.age2$CV.in, main = "CV", names = c("M1", "M2", "W1", "W2"))

par(mfrow = c(1,1))
boxplot(CVrate.men.age1$CV.in, CVrate.men.age1$CV.dir,
        CVrate.men.age2$CV.in, CVrate.men.age2$CV.dir,
        CVrate.women.age1$CV.in, CVrate.women.age1$CV.dir,
        CVrate.women.age2$CV.in, CVrate.women.age2$CV.dir,
        main = "Bootstrap-based CVs", names = c("M1.in", "M1.dir","M2.in","M2.dir", "W1.in", "W1.dir", "W2.in", "W2.dir"))



boxplot(data.frame(CV.design = CV.rate[c(T,F,F,F)], CV.dir = CVrate.men.age1$CV.dir, CV.in = CVrate.men.age1$CV.in), names = c("Desing", "Direct", "Plug-in"))
boxplot(data.frame(CV.design = CV.rate[c(F,T,F,F)], CV.dir = CVrate.men.age2$CV.dir, CV.in = CVrate.men.age2$CV.in), names = c("Desing", "Direct", "Plug-in"))
boxplot(data.frame(CV.design = CV.rate[c(F,F,T,F)], CV.dir = CVrate.women.age1$CV.dir, CV.in = CVrate.women.age1$CV.in), names = c("Desing", "Direct", "Plug-in"))
boxplot(data.frame(CV.design = CV.rate[c(F,F,F,T)], CV.dir = CVrate.women.age2$CV.dir, CV.in = CVrate.women.age2$CV.in), names = c("Desing", "Direct", "Plug-in"))




#### Scatter plots



D <- nrow(data)

df.pred2 <- data %>%
  select(Prov, Age, Sex, nd) %>%      # keep nd if you want point sizes
  slice(rep(1:D, times = 3)) %>%      # repeat each domain 3 times
  mutate(
    Category = factor(rep(c("E","U","I"), each = D), levels = c("E","U","I")),
    Direct   = c(yy[,1],    yy[,2],    yy[,3]),
    Plug_in  = c(mud.in[,1],mud.in[,2],mud.in[,3]),
  )


df.pred2 <- df.pred2 %>%
  mutate(Age = factor(Age), Sex = factor(Sex))



ggplot(df.pred2, aes(x = Plug_in, y = Direct, color = Category)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.8, color = "black") +
  geom_point(aes(size = nd), alpha = 0.55) +
  scale_color_manual(
    name   = "Category",
    values = c(E = "#1b9e77", U = "#d95f02", I = "#7570b3"),
    labels = c(E = "Employed", U = "Unemployed", I = "Inactive")
  ) +
  scale_size_continuous(name = "Sample size", range = c(2.2, 5)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "Plug-in prediction", y = "Direct estimate",
       title = "Labour Categories") +
  theme_classic(base_size = 20) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 20),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 15),
    legend.position = c(0.9, 0.35),
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # bigger points in Category legend
    size  = guide_legend(order = 2)                       # keep size legend separate
  )



df.pred3 <- data %>%
  select(Prov, Age, Sex, nd) %>%      # keep nd if you want point sizes
  mutate(
    Rd.dir = Rd.y,
    Rd.in = Rd.in
  )


df.pred3 <- df.pred3 %>%
  mutate(Age = factor(Age), Sex = factor(Sex))



ggplot(df.pred3, aes(x = Rd.in, y = Rd.dir)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.8, color = "red") +
  geom_point(aes(size = nd), alpha = 0.85, color = "#003366") +
  scale_size_continuous(name = "Sample size", range = c(2.2, 5)) +
  coord_equal(xlim = c(0, 0.5), ylim = c(0, 0.5), expand = FALSE) +
  labs(x = "Plug-in prediction", y = "Direct Estimate",
       title = "Unemployment Rate") +
  theme_classic(base_size = 20) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 20),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 15),
    legend.position = c(0.9, 0.2),
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    size  = guide_legend(order = 2)
  )
