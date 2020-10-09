library(CVXR)
library(MASS)
library(R.matlab)
library(e1071)
rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/MSVM_weak_hard.R")
set.seed(1233)

########################
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

################################################
# n = 400
# mu = list(c(-3,5),c(-4,-1),c(2,-2),c(8,4))
# C=100
# X <- do.call(rbind,lapply(mu, FUN = function(x){ t(apply(mvrnorm(n,c(0,0),diag(1,2)), MARGIN =1, FUN = function(y,z){y/sqrt(sum(y^2)) + z},z=x))}))
# y = c(rep(1,n),rep(2,n),rep(3,n),rep(4,n))
X <- rbind(c(-3,5),c(2,-2),c(8,4),c(-4,-1))
y = c(1,2,3,4)
par(mfrow=c(1,1))
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
NN_beta
round(colSums(NN_beta),4)

max_of_all <- function(w){
  max(L2normsquare(w[2,] - w[1,]), L2normsquare(w[3,] - w[1,]), L2normsquare(w[4,] - w[1,]),
      L2normsquare(w[3,] - w[2,]), L2normsquare(w[4,] - w[2,]), 
      L2normsquare(w[4,] - w[3,]))
}
par(mfrow=c(1,2))
##########max_of_all
rescale_par <- Rescale4_beta(X,y,NN_beta,type = "max_of_all")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta
# 
# 1/L2normsquare(NN_beta_rescale[2,1:2] - NN_beta_rescale[1,1:2])
# 1/L2normsquare(NN_beta_rescale[3,1:2] - NN_beta_rescale[1,1:2])
# 1/L2normsquare(NN_beta_rescale[4,1:2] - NN_beta_rescale[1,1:2])
# 1/L2normsquare(NN_beta_rescale[3,1:2] - NN_beta_rescale[2,1:2])
# 1/L2normsquare(NN_beta_rescale[4,1:2] - NN_beta_rescale[2,1:2]) 
# 1/L2normsquare(NN_beta_rescale[4,1:2] - NN_beta_rescale[3,1:2])


plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, max_of_all:",round(max_of_all(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt4(X,y,type = "max_of_all")
draw_line(NN_beta_rescale[4,] - NN_beta_rescale[2,]+c(0,0,-1))
draw_line(NN_beta_rescale[4,] - NN_beta_rescale[2,], lty = 1)
draw_line(NN_beta_rescale[4,] - NN_beta_rescale[2,]+c(0,0,1))
draw_line(NN_beta_rescale[3,] - NN_beta_rescale[1,]+c(0,0,-1))
draw_line(NN_beta_rescale[3,] - NN_beta_rescale[1,], lty = 1)
draw_line(NN_beta_rescale[3,] - NN_beta_rescale[1,]+c(0,0,1))
draw_line(NN_beta_rescale[4,] - NN_beta_rescale[3,]+c(0,0,1))
draw_line(NN_beta_rescale[4,] - NN_beta_rescale[3,], lty =1)

plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, max_of_all:", round(max_of_all(weak_hard_beta[,1:2]),3)))
max_of_all(weak_hard_beta[,1:2])
# 
# 1/L2normsquare(weak_hard_beta[2,1:2] - weak_hard_beta[1,1:2])
# 1/L2normsquare(weak_hard_beta[3,1:2] - weak_hard_beta[1,1:2])
# 1/L2normsquare(weak_hard_beta[4,1:2] - weak_hard_beta[1,1:2])
# 1/L2normsquare(weak_hard_beta[3,1:2] - weak_hard_beta[2,1:2])
# 1/L2normsquare(weak_hard_beta[4,1:2] - weak_hard_beta[2,1:2]) 
# 1/L2normsquare(weak_hard_beta[4,1:2] - weak_hard_beta[3,1:2])

draw_line(weak_hard_beta[4,] - weak_hard_beta[2,])
draw_line(weak_hard_beta[4,] - weak_hard_beta[2,]+c(0,0,-1))
draw_line(weak_hard_beta[4,] - weak_hard_beta[2,]+c(0,0,1))

draw_line(weak_hard_beta[3,] - weak_hard_beta[1,]+c(0,0,-1))
draw_line(weak_hard_beta[3,] - weak_hard_beta[1,]+c(0,0,1))


##########sum_of_max  
rescale_par <- Rescale4_beta(X,y,NN_beta,type = "sum_of_max")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, sum_of_max:",round(sum_of_max(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt4(X,y,type = "sum_of_max")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, sum_of_max:", round(sum_of_max(weak_hard_beta[,1:2]),3)))



##########max_of_sum
rescale_par <- Rescale4_beta(X,y,NN_beta,type = "max_of_sum")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, max_of_sum:",round(max_of_sum(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt4(X,y,type = "max_of_sum")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, max_of_sum:", round(max_of_sum(weak_hard_beta[,1:2]),3)))

##########Duchi
rescale_par <- Rescale4_beta(X,y,NN_beta,type = "Duchi")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, Duchi:",round(Duchi(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <- MSVM_Weak_Hard_opt4(X,y,type = "Duchi")
# weak_hard_beta/sqrt(sum(weak_hard_beta^2))
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, Duchi:", round(Duchi(weak_hard_beta[,1:2]),3)))

##########New1
par(mfrow=c(1,2))
rescale_par <- Rescale4_beta(X,y,NN_beta,type = "New1")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, New1:",round(New1(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt4(X,y,type = "New1")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, New1,0:", round(New1(weak_hard_beta[,1:2]),3)))

plot_decision_boundary(X,y,beta = rbind(NN_beta_rescale[1,] -NN_beta_rescale[4,],NN_beta_rescale[2,] -NN_beta_rescale[4,],NN_beta_rescale[3,] -NN_beta_rescale[4,]),
                       title = paste("NN, New1,dagger", round(New1(NN_beta_rescale[,1:2]),3)), dagger_rule_w = T)
plot_decision_boundary(X,y,beta = rbind(weak_hard_beta[1,] -weak_hard_beta[4,],weak_hard_beta[2,] -weak_hard_beta[4,],weak_hard_beta[3,] -weak_hard_beta[4,]),
                       title = paste("MSVM, New1,dagger", round(New1(weak_hard_beta[,1:2]),3)), dagger_rule_w = T)








######################################
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
  max(L2normsquare(w[2,] - w[1,]), L2normsquare(w[3,] - w[1,]), 
      L2normsquare(w[3,] - w[2,]))
}
sum_of_max <- function(w){
  sum(max(L2normsquare(w[2,] - w[1,]), L2normsquare(w[3,] - w[1,])), 
      max(L2normsquare(w[1,] - w[2,]),L2normsquare(w[3,] - w[2,])),
      max(L2normsquare(w[1,] - w[3,]), L2normsquare(w[2,] - w[3,])))
}
max_of_sum <- function(w){
  max(L2normsquare(w[2,] - w[1,]) + L2normsquare(w[3,] - w[1,]),
      L2normsquare(w[1,] - w[2,]) + L2normsquare(w[3,] - w[2,]),
      L2normsquare(w[1,] - w[3,]) + L2normsquare(w[2,] - w[3,]))
}
Duchi <- function(w){
  L2normsquare(w)
}
New1 <- function(w){
  sum(L2normsquare(w[1,] - w[3,]), L2normsquare(w[2,] - w[3,]))
}

#############3 class triangle
X = rbind(c(-2,1), c(-2,-1), c(2,1))
y = c(1,2,3)
C = 1
  
par(mfrow=c(1,2))
########Verify the nearest neighbour line are indeed perpendicular lines using Bayes rule of Bivariate Gaussian###########
plot_nearest_neighbour_decision_boundary(X,y, title = "NN")

A <- rbind(c(1,-1,0),
           c(1,0,-1),
           c(0,1,-1))

dif_beta <- rbind(Bayes_rule_Gaussian(X[1,,drop= F], X[2,,drop=F]),
                  Bayes_rule_Gaussian(X[1,,drop= F], X[3,,drop=F]),
                  Bayes_rule_Gaussian(X[2,,drop= F], X[3,,drop=F]))

NN_beta <- ginv(A)%*%dif_beta
NN_beta
round(colSums(NN_beta),4)

plot_decision_boundary(X,y,beta = NN_beta, title = "Bayes rule")

## max_of_all

par(mfrow=c(1,2))
rescale_par <- Rescale3_beta(X,y,NN_beta,type = "max_of_all")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, max_of_all:",round(max_of_all(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt(X,y,type = "max_of_all")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, max_of_all:", round(max_of_all(weak_hard_beta[,1:2]),3)))


## sum_of_max
par(mfrow=c(1,2))
rescale_par <- Rescale3_beta(X,y,NN_beta,type = "sum_of_max")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, sum_of_max:",round(sum_of_max(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt(X,y,type = "sum_of_max")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, sum_of_max:", round(sum_of_max(weak_hard_beta[,1:2]),3)))


## max_of_sum
par(mfrow=c(1,2))
rescale_par <- Rescale3_beta(X,y,NN_beta,type = "max_of_sum")
rescale_par$value
rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, max_of_sum:",round(max_of_sum(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt(X,y,type = "max_of_sum")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, max_of_sum:", round(max_of_sum(weak_hard_beta[,1:2]),3)))

## Duchi

par(mfrow=c(1,2))
rescale_par <- Rescale3_beta(X,y,NN_beta,type = "Duchi")
# rescale_par$value
# rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, Duchi:",round(Duchi(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt(X,y,type = "Duchi")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, Duchi:", round(Duchi(weak_hard_beta[,1:2]),3)))



## New1

par(mfrow=c(1,2))
rescale_par <- Rescale3_beta(X,y,NN_beta,type = "New1")
# rescale_par$value
# rescale_par$gamma
NN_beta_rescale <- rescale_par$gamma*NN_beta

plot_decision_boundary(X,y,beta = NN_beta_rescale, title = paste("NN, New1:",round(New1(NN_beta_rescale[,1:2]),3)))
weak_hard_beta <-  MSVM_Weak_Hard_opt(X,y,type = "New1")
plot_decision_boundary(X,y,beta = weak_hard_beta, title = paste("MSVM, New1,0:", round(New1(weak_hard_beta[,1:2]),3)))

plot_decision_boundary(X,y,beta = rbind(NN_beta_rescale[1,] -NN_beta_rescale[3,],NN_beta_rescale[2,] -NN_beta_rescale[3,]),
                       title = paste("NN, New1,dagger:", round(New1(NN_beta_rescale[,1:2]),3)), dagger_rule_w = T)
plot_decision_boundary(X,y,beta = rbind(weak_hard_beta[1,] -weak_hard_beta[3,],weak_hard_beta[2,] -weak_hard_beta[3,]),
                       title = paste("MSVM, New1,dagger:", round(New1(weak_hard_beta[,1:2]),3)), dagger_rule_w = T)


####################### Example in Japenese guys paper #################### 
X = rbind(c(0,2),c(0,1),c(1,0),c(2,0))
y = c(1,2,2,3)
plot_nearest_neighbour_decision_boundary(X, y, dis = "L2", title = NULL, np_resolution = 500)

plot_decision_boundary(X,y,beta = MSVM_Weak_Hard_opt(X,y,C), title = "hard, max of all")

plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C=100), title = "Duchi")
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1,0")
X = rbind(c(0,2),c(2,0),c(0,1),c(1,0))
y = c(1,2,3,3)
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1,dagger", dagger_rule_w = T)

X = rbind(c(0,1),c(1,0),c(0,2),c(2,0))
y = c(1,1,2,3)
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1,0")
source("~/Desktop/Multiclass Classification/MSVM Code/MSVM_weak_hard.R")
plot_decision_boundary(X,y,beta = MSVM_Weak_Hard_opt(X,y,C), title = "Max sum")
source("~/Desktop/Multiclass Classification/MSVM Code/MSVM_weak_hard.R")

plot_decision_boundary(X,y,beta = MSVM_Weak_Hard_opt(X,y,C), title = "Max all")
plot_decision_boundary(X,y,beta = MSVM_Weak_Hard_opt(X,y,C), title = "Sum max")

#######################################################





