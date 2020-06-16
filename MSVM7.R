library(CVXR)
library(MASS)

rm(list = ls())

source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")
# set.seed(153)
d <- 2
m <- 3
X1 <- mvrnorm(10, c(-2,3), matrix(c(6,0,0,6),nrow=2))
X2 <- mvrnorm(10, c(3,-2), matrix(c(6,0,0,6),nrow=2))
X3 <- mvrnorm(10, c(-3,-3), matrix(c(6,0,0,6),nrow=2))
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
X <- rbind(X1,X2,X3)


w1 <- Variable(d)
w2 <- Variable(d)
b1 <- Variable(1)
b2 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 2)
slack2 <- Variable(rows = nrow(X2), cols = 2)
slack3 <- Variable(rows = nrow(X3), cols = 2)


C <- 1
objective <- Minimize(sum_squares(vstack(w1,w2))/2 +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))
constraints <- list(X1 %*% w1  + b1  >= 1 - sum_entries(slack1, axis=1),
                    X1 %*% w2  + b2 <= -1 + slack1[,2],
                    X2 %*% w1  + b1 <= -1 + slack2[,1],
                    X2 %*% w2  + b2  >= 1 - sum_entries(slack2, axis=1),
                    X3 %*% w1  + b1 <= -1 + slack3[,1],
                    X3 %*% w2  + b2 <= -1 +  slack3[,2],
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

MSVM7 <- Problem(objective, constraints)


CVXR_MSVM7 <- solve(MSVM7)

CVXR_MSVM7_w1 <- CVXR_MSVM7$getValue(w1)
CVXR_MSVM7_b1 <- CVXR_MSVM7$getValue(b1)
CVXR_MSVM7_w2 <- CVXR_MSVM7$getValue(w2)
CVXR_MSVM7_b2 <- CVXR_MSVM7$getValue(b2)

MSVM7_beta1 <- c(CVXR_MSVM7_w1,CVXR_MSVM7_b1)
MSVM7_beta2 <- c(CVXR_MSVM7_w2,CVXR_MSVM7_b2)

MSVM7_primary_beta <- rbind(MSVM7_beta1,MSVM7_beta2)
plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10))

points(X2, col='green')

points(X3, col='blue')

abline(a = -MSVM7_beta1[3]/MSVM7_beta1[2], b = -MSVM7_beta1[1]/MSVM7_beta1[2])
abline(a = -MSVM7_beta2[3]/MSVM7_beta2[2], b = -MSVM7_beta2[1]/MSVM7_beta2[2])

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10))

points(X2, col='green')

points(X3, col='blue')
# 
# MSVM7_beta1_beta2 <- MSVM7_beta1 - MSVM7_beta2
# MSVM7_beta1_beta3 <- MSVM7_beta1 - MSVM7_beta3
# MSVM7_beta2_beta3 <- MSVM7_beta2 - MSVM7_beta3
# 
# 
# 
# abline(a = -MSVM7_beta1_beta2[3]/MSVM7_beta1_beta2[2], b = -MSVM7_beta1_beta2[1]/MSVM7_beta1_beta2[2], lty =2, col = "red")
# abline(a = -MSVM7_beta1_beta3[3]/MSVM7_beta1_beta3[2], b = -MSVM7_beta1_beta3[1]/MSVM7_beta1_beta3[2], lty =2, col = "red")
# abline(a = -MSVM7_beta2_beta3[3]/MSVM7_beta2_beta3[2], b = -MSVM7_beta2_beta3[1]/MSVM7_beta2_beta3[2], lty =2)


base_class=2

class_idx <- sort(unique(y))
m <- length(class_idx)
if(is.null(base_class)){
  Y <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y==id)})
}else{
  Y <- sapply(class_idx[class_idx!=base_class], function(id){as.numeric(y==id)})
}


n <- nrow(X)
p <- ncol(X)

w <- Variable(rows = m-1, cols = p)
b <- Variable(rows = m-1)
slack <- Variable(rows = n, cols = m-1)

objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(slack))


constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(1- slack),
                    sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >=(1- sum_entries(slack, axis=1))*sum_entries(Y,axis=1),
                    slack >=0)

MSVM7 <- Problem(objective, constraints)
CVXR_MSVM7 <- solve(MSVM7, solver = "MOSEK")

cbind(CVXR_MSVM7$getValue(w),CVXR_MSVM7$getValue(b))
MSVM7_primary_beta

