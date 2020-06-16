library(CVXR)
library(MASS)
library(R.matlab)
set.seed(1223)
rm(list = ls())

setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")
p <- 2
m <- 3
X1 <- mvrnorm(20,c(0,5),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(20,c(5,0),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(20,c(-5,-5),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
C <- 1
n <- nrow(X)
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})

############################################# 
w1 <- Variable(p)
w2 <- Variable(p)
b1 <- Variable(1)
b2 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 2)
slack2 <- Variable(rows = nrow(X2), cols = 2)
slack3 <- Variable(rows = nrow(X3), cols = 2)
xi1 <- Variable(rows = nrow(X1))
xi2 <- Variable(rows = nrow(X2))


objective <- Minimize(sum_squares(vstack(w1,w2))/2 +  C*(sum(xi1)+sum(xi2)+sum_entries(slack3)))
constraints <- list(
                    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
                    X1 %*% w1 + b1 >= 1 - xi1,
                    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                    X2 %*% w2 + b2 >= 1 - xi2,
                    X3 %*% w1 + b1 <= -1 + slack3[,1],
                    X3 %*% w2 + b2 <= -1 + slack3[,2],
                    xi1 >= max_entries(slack1, axis = 1)/2,
                    xi2 >= max_entries(slack2, axis = 1)/2,
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

New3 <- Problem(objective, constraints)
CVXR_New3 <- solve(New3, solver = "MOSEK")


beta1 <- c(CVXR_New3$getValue(w1),CVXR_New3$getValue(b1))
beta2 <- c(CVXR_New3$getValue(w2),CVXR_New3$getValue(b2))

New3_primary_beta <- rbind(beta1,beta2)

New3_primary_beta

New3_beta1_beta2 <- beta1 - beta2
New3_beta1_beta3 <- beta1 - 0
New3_beta2_beta3 <- beta2 - 0

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "New3")
points(X2, col='green')
points(X3, col='blue')

abline(a = -New3_beta1_beta2[3]/New3_beta1_beta2[2], b = -New3_beta1_beta2[1]/New3_beta1_beta2[2], lty =2, col = "red")
abline(a = -New3_beta1_beta3[3]/New3_beta1_beta3[2], b = -New3_beta1_beta3[1]/New3_beta1_beta3[2], lty =2, col = "red")
abline(a = -New3_beta2_beta3[3]/New3_beta2_beta3[2], b = -New3_beta2_beta3[1]/New3_beta2_beta3[2], lty =2)


############################################# 
set.seed(1234)
p <- 2
m <- 4
X1 <- mvrnorm(20,c(0,4),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(20,c(4,-1),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(20,c(-4,-1),matrix(c(5,0,0,5),nrow=2))
X4 <- mvrnorm(20,c(0,-4),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3,X4)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X4)))
C <- 1
n <- nrow(X)
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
# writeMat(con="./matlab code/New3-4class.mat", X1 = X1, X2 = X2, X3 = X3, X4 = X4, p=p, m=m, C = C)
w1 <- Variable(p)
w2 <- Variable(p)
w3 <- Variable(p)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)

slack1 <- Variable(rows = nrow(X1), cols = 3)
slack2 <- Variable(rows = nrow(X2), cols = 3)
slack3 <- Variable(rows = nrow(X3), cols = 3)
slack4 <- Variable(rows = nrow(X4), cols = 3)
xi1 <- Variable(rows = nrow(X1))
xi2 <- Variable(rows = nrow(X2))
xi3 <- Variable(rows = nrow(X2))


objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum(xi1)+sum(xi2)+sum(xi3)+sum_entries(slack4)))
constraints <- list(
  X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
  X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,3],
  X1 %*% w1 + b1 >= 1 - xi1,
  X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
  X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,3],
  X2 %*% w2 + b2 >= 1 - xi2,
  X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
  X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
  X3 %*% w3 + b3 >= 1 - xi3,
  X4 %*% w1 + b1 <= -1 + slack4[,1],
  X4 %*% w2 + b2 <= -1 + slack4[,2],
  X4 %*% w3 + b3 <= -1 + slack4[,3],
  xi1 >= max_entries(slack1, axis = 1)/2,
  xi2 >= max_entries(slack2, axis = 1)/2,
  xi3 >= max_entries(slack3, axis = 1)/2,
  
  xi1 >= sum_entries(slack1, axis = 1)/3,
  xi2 >= sum_entries(slack2, axis = 1)/3,
  xi3 >= sum_entries(slack3, axis = 1)/3,
  
  slack1 >=0,
  slack2 >=0,
  slack3 >=0,
  slack4 >=0)

New3 <- Problem(objective, constraints)
CVXR_New3 <- solve(New3, solver = "MOSEK")


beta1 <- c(CVXR_New3$getValue(w1),CVXR_New3$getValue(b1))
beta2 <- c(CVXR_New3$getValue(w2),CVXR_New3$getValue(b2))
beta3 <- c(CVXR_New3$getValue(w3),CVXR_New3$getValue(b3))
New3_primary_beta <- rbind(beta1,beta2,beta3)

New3_primary_beta

####################################

class_idx <- sort(unique(y))
m = length(class_idx)
X_p1 <- X[y!=class_idx[m],,drop = F]
X_p2 <- X[y==class_idx[m],,drop = F]
n_p1 <- nrow(X_p1)
n_p2 <- nrow(X_p2)
Y_p1 <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y[1:n_p1]==id)})

w <- Variable(rows = m-1, cols = p)
b <- Variable(m-1)


slack_p1 <- Variable(rows = n_p1, cols = m-1)
slack_p2 <- Variable(rows = n_p2, cols = m-1)
epsilon_p1 <- Variable(n_p1)
t <- Variable(n_p1*(m-2))
u <- Variable(rows = n_p1*(m-2), cols = m-1)


objective <- Minimize(sum_squares(w)/2 +  C*(sum(epsilon_p1)+sum_entries(slack_p2)))
constraints <- list(((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1)%*%matrix(1, nrow = m-1, ncol = m-1) - (X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b)) >= (1-Y_p1)*(1-slack_p1),
                    vec(t(epsilon_p1%*%matrix(1,ncol=m-2))) >=  rep(seq(1,m-2),n_p1)/rep(seq(2,m-1),n_p1)*t + 1/rep(seq(2,m-1),n_p1)*sum_entries(u,axis=1),
                    t%*%matrix(1,nrow=1,ncol=m-1) + u >=  (diag(1,n_p1)%x%matrix(1,nrow=m-2))%*%slack_p1,
                    u>=0,
                    sum_entries((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1, axis = 1) >= 1-epsilon_p1,
                    X_p2%*%t(w) + matrix(1,nrow=n_p2,ncol =1)%*%t(b) <= -1 + slack_p2,
                    slack_p1 >=0,
                    slack_p2 >=0)


New3 <- Problem(objective, constraints)
CVXR_New3 <- solve(New3, solver = "MOSEK")
New3_primary_beta
cbind(CVXR_New3$getValue(w),CVXR_New3$getValue(b))
