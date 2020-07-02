library(CVXR)
library(MASS)
library(R.matlab)
library(e1071)
library(glmnet)
rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/MSVM_weak_hard.R")
set.seed(1233)
L2normsquare <- function(x){
  sum(x^2)
}

Duchi <- function(w){
  L2normsquare(w)
}

line_sum_norm <- function(beta){
  sum(beta[, 1:2]^2)
}
multinom_ridge_lik <- function(beta,X,Y, lambda){
  m <- nrow(beta)
  sum(apply(((cbind(X,1)%*%t(beta))*Y)%*%matrix(1, nrow = m, ncol = m) - cbind(X,1)%*%t(beta), MARGIN = 1, function(x){log(sum(exp(-x)))})) + lambda/2*sum(beta[,1:2]^2)
}
multinom_lik <- function(beta,X,Y){
  m <- nrow(beta)
  sum(apply(((cbind(X,1)%*%t(beta))*Y)%*%matrix(1, nrow = m, ncol = m) - cbind(X,1)%*%t(beta), MARGIN = 1, function(x){log(sum(exp(-x)))}))
}

minimum_margin <- function(beta, X, Y){
  M <- as.matrix(((cbind(X,1)%*%t(beta))*Y)%*%matrix(1, nrow = m, ncol = m) - cbind(X,1)%*%t(beta))
  diag(M) <- NA
  min(M, na.rm = T)
  # unique(M)
  # M
}

L2square <- function(x){
  sum(x^2)
}

mu1 <- c(-2,1)
mu2 <- c(-2,-1)
mu3 <- c(2,1)
draw_line <- function(beta, lty = 2){
  abline(a = -beta[3]/beta[2], b = -beta[1]/beta[2], lty =lty, col = "red")
}
X1 <- matrix(mvrnorm(20,mu1,diag(0.00001,2)), nrow =20)
X2 <- matrix(mvrnorm(20,mu2,diag(0.00001,2)), nrow =20)
X3 <- matrix(mvrnorm(20,mu3,diag(0.00001,2)), nrow =20)
X <- rbind(X1,X2,X3)

C=1
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
plot(X1, col='red', xlim = c(-3,3), ylim=c(-3,3), xlab = "X1", ylab = "X2", main = "Oracle")
points(X2, col='green')
points(X3, col='blue')

lambda = 1e-4
multi_fit <- glmnet(x=X, y=y, family = "multinomial", alpha = 0, lambda = lambda)
beta <- cbind(t(do.call(cbind,multi_fit$beta)), multi_fit$a0)

class_idx <- sort(unique(y))
m <- length(class_idx)
Y <- sapply(class_idx, function(id){as.numeric(y==id)})

deviance(multi_fit)/2
multinom_lik(beta, X, Y)


X <- rbind(c(-2,1), c(-2,-1), c(2,1))
y = c(1,2,3)
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
n <- nrow(X)
p <- ncol(X)
m <- length(class_idx)

1/minimum_margin(beta/sqrt(sum(beta[,1:2]^2)),X,Y) / sqrt(sum(beta[,1:2]^2))

rescale_par <- Rescale3_beta(X,y,beta,type = "Duchi")
rescale_par$value
rescale_par$gamma
beta_rescale <- rescale_par$gamma*beta



w <- Variable(rows = m, cols = p)
b <- Variable(m)
ml_lik <- sum(do.call(rbind,sapply(1:n, function(i){log_sum_exp(-(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)))[i,])})))
lambda = 1e-5
objective <- Minimize(ml_lik + lambda * sum_squares(w))
constraints <- list( sum_entries(b) == 0,
                     sum_entries(w, axis = 2) ==0)
MN <- Problem(objective, constraints)
CVXR_MN <- solve(MN,  abstol = 1e-16, num_iter = 1e8)
MN_beta <- cbind(CVXR_MN$getValue(w),CVXR_MN$getValue(b))
multinom_ridge_lik(MN_beta, X, Y, lambda)
CVXR_MN$value

rescale_par <- Rescale3_beta(X,y,MN_beta,type = "Duchi")
rescale_par$value
rescale_par$gamma
MN_beta_rescale <- rescale_par$gamma*MN_beta
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
m <- length(class_idx)
  
plot_decision_boundary(X, y, beta_rescale, paste("Multinomial, glmnet:", round(Duchi(beta_rescale[,1:2]),3)))

weak_hard_beta <- MSVM_Weak_Hard_opt(X,y,"Duchi")
weak_hard_beta2 <- MSVM_Weak_Hard_opt(X,y,"Duchi2")
weak_hard_beta3 <- MSVM_Weak_Hard_opt(X,y,"Duchi3")

plot_decision_boundary(X, y, weak_hard_beta, paste("MSVM, Duchi:", round(Duchi(weak_hard_beta[,1:2]),3)))
plot_decision_boundary(X, y, weak_hard_beta2, paste("MSVM, Duchi:", round(Duchi(weak_hard_beta2[,1:2]),3)))
plot_decision_boundary(X, y, weak_hard_beta3, paste("MSVM, Duchi:", round(Duchi(weak_hard_beta3[,1:2]),3)))
draw_line(beta[1,]-beta[2,])
draw_line(beta[1,]-beta[3,])
draw_line(beta[2,]-beta[3,])
    
plot_decision_boundary(X, y, MN_beta_rescale, paste("MN_beta", round(Duchi(MN_beta_rescale[,1:2]),3)))
plot_decision_boundary(X, y, weak_hard_beta, paste("MSVM, Duchi:", round(Duchi(weak_hard_beta[,1:2]),3)))
draw_line(MN_beta_rescale[1,]-MN_beta_rescale[2,])
draw_line(MN_beta_rescale[1,]-MN_beta_rescale[3,])
draw_line(MN_beta_rescale[2,]-MN_beta_rescale[3,])
# 
multinom_lik(25*MN_beta_rescale/sqrt(sum(MN_beta_rescale[,1:2]^2)), X, Y)
# multinom_lik(beta/sqrt(sum(beta[1:2,]^2)), X, Y)
multinom_lik(29.7*beta_rescale/sqrt(sum(beta_rescale[,1:2]^2)), X, Y)
multinom_lik(29.7*weak_hard_beta/sqrt(sum(weak_hard_beta[,1:2]^2)), X, Y)

multinom_ridge_lik(20*beta/sqrt(sum(beta[,1:2]^2)), X, Y, lambda)
multinom_ridge_lik(20*beta_rescale/sqrt(sum(beta_rescale[,1:2]^2)), X, Y, lambda)
multinom_ridge_lik(20*weak_hard_beta/sqrt(sum(weak_hard_beta[,1:2]^2)), X, Y, lambda)
# 
Duchi(weak_hard_beta[,1:2])
Duchi(beta_rescale[,1:2])

line_sum_norm(weak_hard_beta/sqrt(sum(weak_hard_beta[,1:2]^2)))
line_sum_norm(beta/sqrt(sum(beta[,1:2]^2)))

minimum_margin(MN_beta_rescale/sqrt(sum(MN_beta_rescale[,1:2]^2)), X,Y)
minimum_margin(beta_rescale/sqrt(sum(beta_rescale[,1:2]^2)), X,Y)
minimum_margin(weak_hard_beta/sqrt(sum(weak_hard_beta[,1:2]^2)), X,Y)
minimum_margin(weak_hard_beta2/sqrt(sum(weak_hard_beta2[,1:2]^2)), X,Y)
minimum_margin(weak_hard_beta3/sqrt(sum(weak_hard_beta3[,1:2]^2)), X,Y)

minimum_margin(weak_hard_beta, X,Y)
minimum_margin(beta_rescale, X,Y)
sum((weak_hard_beta[,1:2])^2)
sum((beta_rescale[,1:2])^2)

minimum_margin(weak_hard_beta2/sqrt(sum(weak_hard_beta2[1:2,]^2)),X,Y)

1/1.251304^2
#################################################
####################################################
##################################################
# Sigma <- matrix(c(4,2,2,4), nrow = 2)
# eigen(Sigma)$vectors[1,]
# 
# 
# 
# X1 <- mvrnorm(50,c(0,0),Sigma)
# X2 <- mvrnorm(50,5*eigen(Sigma)$vectors[1,], Sigma)
# X <- rbind(X1,X2)
# y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)))
# plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "Oracle")
# points(X2, col='green')

mu1 <- c(-3,5)
mu2 <- c(2,-2)
mu3 <- c(8,4)
mu4 <- c(-4,-1)
v = 0.0001
X1 <- matrix(mvrnorm(20, mu1, diag(v,2)), nrow =20)
X2 <- matrix(mvrnorm(20, mu2, diag(v,2)), nrow =20)
X3 <- matrix(mvrnorm(20, mu3, diag(v,2)), nrow =20)
X4 <- matrix(mvrnorm(20, mu4, diag(v,2)), nrow =20)
X <- rbind(X1,X2,X3,X4)

C=1
y <- c(rep(1,nrow(X1)), rep(2,nrow(X2)), rep(3,nrow(X3)), rep(4,nrow(X4)))

multi_fit <- glmnet(x=X, y=y, family = "multinomial", alpha = 0, lambda = 1e-2)
beta <- cbind(t(do.call(cbind,multi_fit$beta)), multi_fit$a0)

L2square <- function(x){
  sum(x^2)
}
# beta <- beta/(L2square(beta[1,1:2])+L2square(beta[2,1:2])+L2square(beta[3,1:2])+L2square(beta[4,1:2]))

beta - colSums(beta)

rescale_par <- Rescale4_beta(X,y,beta,type = "Duchi")
rescale_par$value
rescale_par$gamma
beta_rescale <- rescale_par$gamma*beta
# 
  
plot_decision_boundary(X, y, beta_rescale, paste("Multinomial, Duchi:", round(Duchi(beta_rescale[,1:2]),3)))

weak_hard_beta <- MSVM_Weak_Hard_opt4(X,y,"Duchi")

plot_decision_boundary(X, y, weak_hard_beta, paste("MSVM, Duchi:", round(Duchi(weak_hard_beta[,1:2]),3)))

# rescale_par <- Rescale4_beta(X,y,weak_hard_beta,type = "Duchi")
# rescale_par$value
# rescale_par$gamma
# weak_hard_beta_rescale <- rescale_par$gamma*weak_hard_beta
# # 
# 
# plot_decision_boundary(X, y, weak_hard_beta_rescale, paste("MSVM, Duchi:", round(Duchi(weak_hard_beta_rescale[,1:2]),3)))

# plot_decision_boundary(X, y, MSVM_Weak_Hard_opt4(X,y,"New1"), title = "MSVM")
draw_line(beta = beta[1,]-beta[2,])
draw_line(beta = beta[1,]-beta[3,])
draw_line(beta = beta[1,]-beta[4,])
draw_line(beta = beta[2,]-beta[3,])
draw_line(beta = beta[2,]-beta[4,])
draw_line(beta = beta[3,]-beta[4,])

X <- rbind(mu1,mu2, mu3,mu4)
v = 0.001
y <- c(1,2,3,4)
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
n <- nrow(X)
p <- ncol(X)
m <- length(class_idx)



w <- Variable(rows = m, cols = p)
b <- Variable(m)
ml_lik <- sum(do.call(rbind,sapply(1:n, function(i){log_sum_exp(-(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)))[i,])})))
lambda = 1e-3
objective <- Minimize(ml_lik + lambda * sum_squares(w))
constraints <- list(sum_entries(b) == 0,
                    sum_entries(w, axis = 2) ==0
                    )
MN <- Problem(objective, constraints)
CVXR_MN <- solve(MN,  abstol = 1e-8)
# ?solve
MN_beta <- cbind(CVXR_MN$getValue(w),CVXR_MN$getValue(b))
# multinom_ridge_lik(MN_beta, X, Y, lambda)
multinom_lik(MN_beta, X, Y)
CVXR_MN$value
line_sum_norm(MN_beta)

multinom_lik(10*MN_beta/sqrt(sum(MN_beta[,1:2]^2)), X, Y)
multinom_lik(10*weak_hard_beta/sqrt(sum(weak_hard_beta[,1:2]^2)), X, Y)
multinom_lik(10*beta/sqrt(sum(beta[,1:2]^2)), X, Y)
# multinom_lik(10*NN_beta/sqrt(sum(NN_beta[1:2,]^2)), X, Y)


line_sum_norm(MN_beta/sqrt(sum(MN_beta[,1:2]^2)))
line_sum_norm(weak_hard_beta/sqrt(sum(weak_hard_beta[,1:2]^2)))
line_sum_norm(beta/sqrt(sum(beta[,1:2]^2)))
# line_sum_norm(NN_beta/sqrt(sum(NN_beta[,1:2]^2)))
minimum_margin(MN_beta/sqrt(sum(MN_beta[,1:2]^2)), X,Y)
minimum_margin(beta_rescale/sqrt(sum(beta_rescale[,1:2]^2)), X,Y)
minimum_margin(weak_hard_beta/sqrt(sum(weak_hard_beta[,1:2]^2)), X,Y)


# 
# rescale_par <- Rescale4_beta(X,y,MN_beta,type = "Duchi")
# rescale_par$value
# rescale_par$gamma
# MN_beta_rescale <- rescale_par$gamma*MN_beta

plot_decision_boundary(X, y, MN_beta, paste("Multinomial:", round(Duchi(MN_beta[,1:2]),3)))

CVXR_MN$value

