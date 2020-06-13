library(nnet)


setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")

# set.seed(122321)
# set.seed(1234)
p <- 2
m <- 4
Sigma <- diag(1,nrow=2)
  n <- 10
X1 <- mvrnorm(n,c(0,2),Sigma)
X2 <- mvrnorm(n,c(2,0),Sigma)
X3 <- mvrnorm(n,c(-2,0),Sigma)
X4 <- mvrnorm(n,c(0,-2),Sigma)
C <- 1
X <- rbind(X1,X2,X3,X4)
y <- c(rep(1,nrow(X1)), rep(2,nrow(X2)), rep(3,nrow(X3)), rep(4,nrow(X4)))
oracle <- multinom(y ~ X, trace = F)
mean(predict(oracle, X) == y)

plot(X1,  col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "Data")
points(X2, col='green')
points(X3, col='blue')
points(X4, col='yellow')

writeMat(con="Duchi-4class.mat", X1 = X1, X2 = X2, X3 = X3, X4 = X4, p=p, m=m, C = C)
writeMat(con="MDuchi-4class.mat", X1 = X1, X2 = X2, X3 = X3, X4 = X4, p=p, m=m, C = C)

Duchi_dual_opt(X,y,C)
MDuchi_dual_opt(X,y,C)


WW_dual_opt(X,y,C)

