library(CVXR)
library(MASS)
library(Rmosek)


rm(list = ls())
set.seed(153)
p <- 2
m <- 3
v <- 3
X1 <- mvrnorm(10, c(-2,3), diag(v,p))
X2 <- mvrnorm(10, c(3,-2), diag(v,p))
X3 <- mvrnorm(10, c(-3,-3), diag(v,p))
# 
# m <- 3
# X1 <- mvrnorm(10, c(-2,3,3), matrix(c(8,0,0,0,8,0,0,0,8),nrow=3))
# X2 <- mvrnorm(10, c(3,-2,2),  matrix(c(8,0,0,0,8,0,0,0,8),nrow=3))
# X3 <- mvrnorm(10, c(-3,-3,3),  matrix(c(8,0,0,0,8,0,0,0,8),nrow=3))
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
# 

w1 <- Variable(p)
w2 <- Variable(p)
w3 <- Variable(p)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 3)
slack2 <- Variable(rows = nrow(X2), cols = 3)
slack3 <- Variable(rows = nrow(X3), cols = 3)


C <- 1

objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))
constraints <- list(X1 %*% w1  + b1  >= 1 - slack1[,1],
                    X1 %*% w2  + b2 <= -1 + slack1[,2],
                    X1 %*% w3  + b3 <= -1 + slack1[,3],
                    X2 %*% w1  + b1 <= -1 + slack2[,1],
                    X2 %*% w2  + b2  >= 1 - slack2[,2],
                    X2 %*% w3  + b3 <= -1 + slack2[,3],
                    X3 %*% w1  + b1 <= -1 + slack3[,1],
                    X3 %*% w2  + b2 <= -1 +  slack3[,2],
                    X3 %*% w3  + b3  >= 1- slack3[,3], 
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

OVA <- Problem(objective, constraints)


CVXR_OVA <- solve(OVA, solver="MOSEK")

CVXR_OVA_w1 <- CVXR_OVA$getValue(w1)
CVXR_OVA_b1 <- CVXR_OVA$getValue(b1)
CVXR_OVA_w2 <- CVXR_OVA$getValue(w2)
CVXR_OVA_b2 <- CVXR_OVA$getValue(b2)
CVXR_OVA_w3 <- CVXR_OVA$getValue(w3)
CVXR_OVA_b3 <- CVXR_OVA$getValue(b3)

CVXR_OVA_w <- t(cbind(CVXR_OVA_w1,CVXR_OVA_w2,CVXR_OVA_w3))
OVA_beta1 <- c(CVXR_OVA_w1,CVXR_OVA_b1)
OVA_beta2 <- c(CVXR_OVA_w2,CVXR_OVA_b2)
OVA_beta3 <- c(CVXR_OVA_w3,CVXR_OVA_b3)


##############################

X <- rbind(X1,X2,X3)
n <- nrow(X)
m <- length(unique(y))
p <- ncol(X)

Y <- sapply(unique(y), function(id){as.numeric(y==id)})
Delta <- matrix(1, nrow  = n, ncol = m) - 2*Y
# Delta <-  2*Y - matrix(1, nrow  = n, ncol = m) 
alpha <- Variable(rows = n, cols = m)

objective <- Maximize((-sum_squares(t(alpha*Delta)%*%X)/2 +  sum_entries(alpha)))

C = 1
constraints <- list(sum_entries(alpha*Delta, axis = 2) ==0,
                    alpha<=C,
                    # alpha*(1-Y) + sum_entries(alpha*Y,axis=1)%*%matrix(1,ncol=m,nrow=1)<=C,
                    alpha>=0) 


OVA_dual <- Problem(objective, constraints)
CVXR_OVA_dual <- solve(OVA_dual)
alpha_sol <- CVXR_OVA_dual$getValue(alpha)


cbind(CVXR_OVA$getDualValue(constraints(OVA)[[1]]),CVXR_OVA$getDualValue(constraints(OVA)[[2]]),CVXR_OVA$getDualValue(constraints(OVA)[[3]]))

round(alpha_sol[1:10,] ,2)

w <- -t(alpha_sol*Delta)%*%X
CVXR_OVA_w
w

CVXR_OVA_w[1,2]/CVXR_OVA_w[1,1]

round(Y*alpha_sol,1)
round((1-Y)*alpha_sol,1)

  lpred <- X%*%t(w)

b <- sapply(unique(y), function(id){-(max(lpred[((C - Y*alpha_sol) > .01 & Y*alpha_sol > 0.01)[,id] == T,id]) + min(lpred[((C - (1-Y)*alpha_sol) > .01 & (1-Y)*alpha_sol > 0.01)[,id] == T,id])) / 2})

rbind(OVA_beta1,OVA_beta2,OVA_beta3)

beta <- cbind(w,b)
beta

par(mfrow=c(2,2))
plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "OVA")
points(X2, col='green')
points(X3, col='blue')
abline(a = -OVA_beta1[3]/OVA_beta1[2], b = -OVA_beta1[1]/OVA_beta1[2])
abline(a = -OVA_beta2[3]/OVA_beta2[2], b = -OVA_beta2[1]/OVA_beta2[2])
abline(a = -OVA_beta3[3]/OVA_beta3[2], b = -OVA_beta3[1]/OVA_beta3[2])
plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "OVA")
points(X2, col='green')
points(X3, col='blue')
OVA_beta1_beta2 <- OVA_beta1 - OVA_beta2
OVA_beta1_beta3 <- OVA_beta1 - OVA_beta3
OVA_beta2_beta3 <- OVA_beta2 - OVA_beta3
abline(a = -OVA_beta1_beta2[3]/OVA_beta1_beta2[2], b = -OVA_beta1_beta2[1]/OVA_beta1_beta2[2], lty =2, col = "red")
abline(a = -OVA_beta1_beta3[3]/OVA_beta1_beta3[2], b = -OVA_beta1_beta3[1]/OVA_beta1_beta3[2], lty =2, col = "red")
abline(a = -OVA_beta2_beta3[3]/OVA_beta2_beta3[2], b = -OVA_beta2_beta3[1]/OVA_beta2_beta3[2], lty =2)

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "OVA")
points(X2, col='green')
points(X3, col='blue')
abline(a = -beta[1,3]/beta[1,2], b = -beta[1,1]/beta[1,2])
abline(a = -beta[2,3]/beta[2,2], b = -beta[2,1]/beta[2,2])
abline(a = -beta[3,3]/beta[3,2], b = -beta[3,1]/beta[3,2])
plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "OVA")
points(X2, col='green')
points(X3, col='blue')
beta1_beta2 <- beta[1,] - beta[2,]
beta1_beta3 <- beta[1,] - beta[3,]
beta2_beta3 <- beta[2,] - beta[3,]
abline(a = -beta1_beta2[3]/beta1_beta2[2], b = -beta1_beta2[1]/beta1_beta2[2], lty =2, col = "red")
abline(a = -beta1_beta3[3]/beta1_beta3[2], b = -beta1_beta3[1]/beta1_beta3[2], lty =2, col = "red")
abline(a = -beta2_beta3[3]/beta2_beta3[2], b = -beta2_beta3[1]/beta2_beta3[2], lty =2)




