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

