library(CVXR)
library(MASS)

rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
set.seed(125)
par(mfrow=c(1,1))
d <- 2
m <- 3
mu1 <- c(-3,3)
mu2 <- c(3,-3)
mu3 <- c(-3,-2)
X1 <- mvrnorm(10,mu1,matrix(c(10,0,0,10),nrow=2))
X2 <- mvrnorm(10,mu2,matrix(c(10,0,0,10),nrow=2))
X3 <- mvrnorm(10,mu3,matrix(c(10,0,0,10),nrow=2))
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
X <- rbind(X1,X2,X3)
# X1 = matrix(c(-2,1), nrow = 1)
# X2 = matrix(c(-2,-1), nrow =1)
# X3 = matrix(c(2,1), nrow = 1)
# y = c(1,2,3)

#
# Oracle_beta1_beta2 <- c(2*(mu1-mu2), sum(mu1^2)-sum(mu2^2))
# Oracle_beta1_beta3 <- c(2*(mu1-mu3), sum(mu1^2)-sum(mu3^2))
# Oracle_beta2_beta3 <- c(2*(mu2-mu3), sum(mu2^2)-sum(mu3^2))
#
# plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "Oracle")
#
# points(X2, col='green')
#
# points(X3, col='blue')
#
# abline(a = -Oracle_beta1_beta2[3]/Oracle_beta1_beta2[2], b = -Oracle_beta1_beta2[1]/Oracle_beta1_beta2[2], lty =2, col = "red")
# abline(a = -Oracle_beta1_beta3[3]/Oracle_beta1_beta3[2], b = -Oracle_beta1_beta3[1]/Oracle_beta1_beta3[2], lty =2, col = "red")
# abline(a = -Oracle_beta2_beta3[3]/Oracle_beta2_beta3[2], b = -Oracle_beta2_beta3[1]/Oracle_beta2_beta3[2], lty =2)
#


w1 <- Variable(d)
w2 <- Variable(d)
w3 <- Variable(d)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1))
slack2 <- Variable(rows = nrow(X2))
slack3 <- Variable(rows = nrow(X3))

# 
# g <- function(alpha){
#   max(max(alpha, 1)-1, sum_largest(alpha, 2)/2 - 1/2 , sum(alpha)/3-1/3)
# }
# 
# gslack1 <- lapply(seq(1,nrow(X1)), FUN = function(i){g(slack1[i,])})
# gslack2 <- lapply(seq(1,nrow(X2)), FUN = function(i){g(slack2[i,])})
# gslack3 <- lapply(seq(1,nrow(X3)), FUN = function(i){g(slack3[i,])})


C <- 1
objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))

constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1,
                    X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1,
                    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2,
                    X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2,
                    X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3,
                    X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3,
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

CS <- Problem(objective, constraints)
CVXR_CS <- solve(CS)

CVXR_CS_w1 <- CVXR_CS$getValue(w1)
CVXR_CS_b1 <- CVXR_CS$getValue(b1)
CVXR_CS_w2 <- CVXR_CS$getValue(w2)
CVXR_CS_b2 <- CVXR_CS$getValue(b2)
CVXR_CS_w3 <- CVXR_CS$getValue(w3)
CVXR_CS_b3 <- CVXR_CS$getValue(b3)

beta1 <- c(CVXR_CS_w1,CVXR_CS_b1)
beta2 <- c(CVXR_CS_w2,CVXR_CS_b2)
beta3 <- c(CVXR_CS_w3,CVXR_CS_b3)

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "CS")

points(X2, col='green')

points(X3, col='blue')

abline(a = -beta1[3]/beta1[2], b = -beta1[1]/beta1[2])
abline(a = -beta2[3]/beta2[2], b = -beta2[1]/beta2[2])
abline(a = -beta3[3]/beta3[2], b = -beta3[1]/beta3[2])
rbind(beta1,beta2,beta3)

# CS_pri_opt(X,y,C)
 

CS_beta1_beta2 <- beta1 - beta2
CS_beta1_beta3 <- beta1 - beta3
CS_beta2_beta3 <- beta2 - beta3


abline(a = -CS_beta1_beta2[3]/CS_beta1_beta2[2], b = -CS_beta1_beta2[1]/CS_beta1_beta2[2], lty =2, col = "red")
abline(a = -CS_beta1_beta3[3]/CS_beta1_beta3[2], b = -CS_beta1_beta3[1]/CS_beta1_beta3[2], lty =2, col = "red")
abline(a = -CS_beta2_beta3[3]/CS_beta2_beta3[2], b = -CS_beta2_beta3[1]/CS_beta2_beta3[2], lty =2)


################################################

class_idx <- sort(unique(y))

Y <- sapply(class_idx, function(id){as.numeric(y==id)})
n <- nrow(X)
p <- ncol(X)
m <- length(class_idx)

K <- kernelMatrix(vanilladot(), X)
suppressWarnings(class(K) <- "matrix")

v <- Variable(n,m)
b <- Variable(m)
slack <- Variable(rows = n, cols = m)

objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2  + C*sum_entries(max_entries((1-Y)*slack, axis = 1)))
constraints <- list( ((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (K%*%v  + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                     sum_entries(b)*matrix(1,nrow=n) + K %*% sum_entries(v, axis = 1) == 0,
                     slack >=0)


CS <- Problem(objective, constraints)
CVXR_CS <- solve(CS, solver = "MOSEK")
CVXR_CS_v <- CVXR_CS$getValue(v)
CVXR_CS_b <- CVXR_CS$getValue(b)

model <- list(v = CVXR_CS_v, b = as.numeric(CVXR_CS_b), kernel = kernel, X = X, y= y)
class(model) <- "msvm_kernel"

CS_fit <- CS_kernel_pri_opt(X,y,C=1, kernel = vanilladot())

cbind(t(CS_fit$v)%*%X,CS_fit$b)


