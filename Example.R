library(CVXR)
library(MASS)
set.seed(1223)
rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")
p <- 2
m <- 3
X1 <- mvrnorm(20,c(0,5),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(20,c(5,0),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(20,c(0,0),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
C <- 100


Duchi_pri_opt(X,y,C)


plot_decision_boundary(X,y,beta = WW_pri_opt(X,y,C), title = "WW")
plot_decision_boundary(X,y,beta = CS_pri_opt(X,y,C), title = "CS")
plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C), title = "Duchi")
plot_decision_boundary(X,y,beta = MDuchi_pri_opt(X,y,C), title = "Mduchi")
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1")
plot_decision_boundary(X,y,beta = New3_pri_opt(X,y,C), title = "New3")

plot_decision_boundary(X,y,beta = OVA_pri_opt(X,y,C), title = "OVA")
plot_decision_boundary(X,y,beta = LLW_pri_opt(X,y,C), title = "LLW")
plot_decision_boundary(X,y,beta = MSVM7_pri_opt(X,y,C), title = "MSVM7")
plot_decision_boundary(X,y,beta = MSVM8_pri_opt(X,y,C), title = "MSVM8")


LLW_pri_opt(X,y,C)

###########4 class 
set.seed(1223)
# set.seed(1234)
p <- 2
m <- 4
Sigma <- diag(4,nrow=2)
X1 <- mvrnorm(10,c(0,2),Sigma)
X2 <- mvrnorm(10,c(2,0),Sigma)
X3 <- mvrnorm(10,c(-2,0),Sigma)
X4 <- mvrnorm(10,c(0,-2),Sigma)
C <- 1
X <- rbind(X1,X2,X3,X4)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X3)))



plot_decision_boundary(X,y,beta = WW_pri_opt(X,y,C), title = "WW")
plot_decision_boundary(X,y,beta = CS_pri_opt(X,y,C), title = "CS")
plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C), title = "Duchi")
plot_decision_boundary(X,y,beta = MDuchi_pri_opt(X,y,C), title = "Mduchi")
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1")
plot_decision_boundary(X,y,beta = New3_pri_opt(X,y,C), title = "New3")

plot_decision_boundary(X,y,beta = OVA_pri_opt(X,y,C), title = "OVA")
plot_decision_boundary(X,y,beta = LLW_pri_opt(X,y,C), title = "LLW")
plot_decision_boundary(X,y,beta = MSVM7_pri_opt(X,y,C), title = "MSVM7", xlim = c(-5,5), ylim = c(-8,8))
plot_decision_boundary(X,y,beta = MSVM8_pri_opt(X,y,C), title = "MSVM8")

