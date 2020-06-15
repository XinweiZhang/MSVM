rm(list = ls())
# set.seed(1211)


setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
d <- 2
m <- 3
X1 <- mvrnorm(10, c(-2,3), matrix(c(6,0,0,6),nrow=2))
X2 <- mvrnorm(10, c(3,-2), matrix(c(6,0,0,6),nrow=2))
X3 <- mvrnorm(10, c(-3,-3), matrix(c(6,0,0,6),nrow=2))

y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
X <- rbind(X1,X2,X3)

w1 <- Variable(d)
w2 <- Variable(d)
w3 <- Variable(d)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 2)
slack2 <- Variable(rows = nrow(X2), cols = 2)
slack3 <- Variable(rows = nrow(X3), cols = 2)


C <- 1

objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))
constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                    X1 %*% w2  + b2 <= -1 + slack1[,1],
                    X1 %*% w3  + b3 <= -1 + slack1[,2],
                    X2 %*% w1  + b1 <= -1 + slack2[,1],
                    X2 %*% w3  + b3 <= -1 + slack2[,2],
                    X3 %*% w1  + b1 <= -1 + slack3[,1],
                    X3 %*% w2  + b2 <= -1 +  slack3[,2],
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

LLW <- Problem(objective, constraints)


CVXR_LLW <- solve(LLW)

CVXR_LLW_w1 <- CVXR_LLW$getValue(w1)
CVXR_LLW_b1 <- CVXR_LLW$getValue(b1)
CVXR_LLW_w2 <- CVXR_LLW$getValue(w2)
CVXR_LLW_b2 <- CVXR_LLW$getValue(b2)
CVXR_LLW_w3 <- CVXR_LLW$getValue(w3)
CVXR_LLW_b3 <- CVXR_LLW$getValue(b3)

LLW_beta1 <- c(CVXR_LLW_w1,CVXR_LLW_b1)
LLW_beta2 <- c(CVXR_LLW_w2,CVXR_LLW_b2)
LLW_beta3 <- c(CVXR_LLW_w3,CVXR_LLW_b3)

LLW_primary_beta <- rbind(LLW_beta1,LLW_beta2,LLW_beta3)


#########################################Upgraded
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
n <- nrow(X)
p <- ncol(X)
m <- length(class_idx)

w <- Variable(rows = m, cols = p)
b <- Variable(rows = m)
slack <- Variable(rows = n, cols = m)


objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(slack))

constraints <- list(sum_entries(b) == 0,
                    sum_entries(w, axis = 2) ==0,
                    ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1- slack) * (1-Y),
                    slack >=0)

LLW <- Problem(objective, constraints)
CVXR_LLW <- solve(LLW)
LLW_primary_beta
cbind(CVXR_LLW$getValue(w),CVXR_LLW$getValue(b))

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "LLW")

points(X2, col='green')

points(X3, col='blue')

abline(a = -LLW_beta1[3]/LLW_beta1[2], b = -LLW_beta1[1]/LLW_beta1[2])
abline(a = -LLW_beta2[3]/LLW_beta2[2], b = -LLW_beta2[1]/LLW_beta2[2])
abline(a = -LLW_beta3[3]/LLW_beta3[2], b = -LLW_beta3[1]/LLW_beta3[2])

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "LLW")

points(X2, col='green')

points(X3, col='blue')

LLW_beta1_beta2 <- LLW_beta1 - LLW_beta2
LLW_beta1_beta3 <- LLW_beta1 - LLW_beta3
LLW_beta2_beta3 <- LLW_beta2 - LLW_beta3



abline(a = -LLW_beta1_beta2[3]/LLW_beta1_beta2[2], b = -LLW_beta1_beta2[1]/LLW_beta1_beta2[2], lty =2, col = "red")
abline(a = -LLW_beta1_beta3[3]/LLW_beta1_beta3[2], b = -LLW_beta1_beta3[1]/LLW_beta1_beta3[2], lty =2, col = "red")
abline(a = -LLW_beta2_beta3[3]/LLW_beta2_beta3[2], b = -LLW_beta2_beta3[1]/LLW_beta2_beta3[2], lty =2)
