library(CVXR)
library(MASS)
library(R.matlab)
library(e1071)
rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/MSVM_weak_hard.R")
set.seed(1233)



########################
line_norm <- function(beta){
  A <- matrix(0, nrow = nrow(beta), ncol=nrow(beta))
  for(i in 1:nrow(beta)){
    for(j in 1:nrow(beta)){
      A[i,j] <- L2norm(beta[i,1:2]-beta[j,1:2])
    }
  }
  return(A)
}

L2norm <- function(x){
  sqrt(sum(x^2))
}
L2normsquare <- function(x){
  sum(x^2)
}
draw_line <- function(beta, lty = 2){
  abline(a = -beta[3]/beta[2], b = -beta[1]/beta[2], lty =lty, col = "red")
}
Bayes_rule_Gaussian <- function(mu1,mu2){
  c((mu1-mu2), 1/2*sum(mu2^2)-1/2*sum(mu1^2))
}

max_of_all <- function(w){
  max(L2normsquare(w[2,] - w[1,]), L2normsquare(w[3,] - w[1,]), L2normsquare(w[4,] - w[1,]),
      L2normsquare(w[3,] - w[2,]), L2normsquare(w[4,] - w[2,]), 
      L2normsquare(w[4,] - w[3,]))
}
sum_of_max <- function(w){
  sum(max(L2normsquare(w[2,] - w[1,]), L2normsquare(w[3,] - w[1,]), L2normsquare(w[4,] - w[1,])), 
      max(L2normsquare(w[1,] - w[2,]),L2normsquare(w[3,] - w[2,]), L2normsquare(w[4,] - w[2,])),
      max(L2normsquare(w[1,] - w[3,]), L2normsquare(w[2,] - w[3,]), L2normsquare(w[4,] - w[3,])), 
      max(L2normsquare(w[4,] - w[1,]), L2normsquare(w[4,] - w[2,]), L2normsquare(w[4,] - w[3,])))
}
max_of_sum <- function(w){
  max(L2normsquare(w[2,] - w[1,]) + L2normsquare(w[3,] - w[1,]) + L2normsquare(w[4,] - w[1,]),
      L2normsquare(w[1,] - w[2,]) + L2normsquare(w[3,] - w[2,]) + L2normsquare(w[4,] - w[2,]),
      L2normsquare(w[1,] - w[3,]) + L2normsquare(w[2,] - w[3,]) + L2normsquare(w[4,] - w[3,]), 
      L2normsquare(w[4,] - w[1,]) + L2normsquare(w[4,] - w[2,]) + L2normsquare(w[4,] - w[3,]))
}
Duchi <- function(w){
  L2normsquare(w)
}
New1 <- function(w){
  sum(L2normsquare(w[1,] - w[4,]), L2normsquare(w[2,] - w[4,]), L2normsquare(w[3,] - w[4,]))
}

X <- rbind(c(-3,5),c(-4,-1),c(2,-2),c(8,4))
y = c(1,2,3,4)
par(mfrow=c(1,2))
########Verify the nearest neighbour line are indeed perpendicular lines using Bayes rule of Bivariate Gaussian###########
plot_nearest_neighbour_decision_boundary(X,y, title = "NN")
draw_line(Bayes_rule_Gaussian(X[1,,drop= F], X[2,,drop=F]))
draw_line(Bayes_rule_Gaussian(X[1,,drop= F], X[3,,drop=F]))
draw_line(Bayes_rule_Gaussian(X[1,,drop= F], X[4,,drop=F]))
draw_line(Bayes_rule_Gaussian(X[2,,drop= F], X[3,,drop=F]))
draw_line(Bayes_rule_Gaussian(X[2,,drop= F], X[4,,drop=F]))
draw_line(Bayes_rule_Gaussian(X[3,,drop= F], X[4,,drop=F]))


A <- rbind(c(1,-1,0,0),
           c(1,0,-1,0),
           c(1,0,0,-1),
           c(0,1,-1,0),
           c(0,1,0,-1),
           c(0,0,1,-1))

dif_beta <- rbind(Bayes_rule_Gaussian(X[1,,drop= F], X[2,,drop=F]),
                  Bayes_rule_Gaussian(X[1,,drop= F], X[3,,drop=F]),
                  Bayes_rule_Gaussian(X[1,,drop= F], X[4,,drop=F]),
                  Bayes_rule_Gaussian(X[2,,drop= F], X[3,,drop=F]),
                  Bayes_rule_Gaussian(X[2,,drop= F], X[4,,drop=F]),
                  Bayes_rule_Gaussian(X[3,,drop= F], X[4,,drop=F]))
  
NN_beta <- ginv(A)%*%dif_beta
# NN_beta
# round(colSums(NN_beta),4)
# 
# plot_decision_boundary(X,y,beta = NN_beta, title = "Bayes rule")
constraints_value <- ((cbind(X,1)%*%t(NN_beta)*diag(1,4))%*%matrix(1, ncol = 4,nrow=4) - cbind(X,1)%*%t(NN_beta))
gamma <- min(constraints_value/line_norm(NN_beta), na.rm = T)

NN_beta <- NN_beta/gamma*min(1/line_norm(NN_beta))

((cbind(X,1)%*%t(NN_beta)*diag(1,4))%*%matrix(1, ncol = 4,nrow=4) - cbind(X,1)%*%t(NN_beta))

