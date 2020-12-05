library(e1071)
library(CVXR)
library(kernlab)
library(MASS)
rm(list=ls())
rm(list = ls())
par(mfrow=c(1,1))
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
# source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/SDCA functions.R")


PL_WW_Loss <- function(W, K, Y, lambda){
  s <- W%*%K
  PL_WW_V <- s - t(replicate(m, s[Y == 1])) + (1-Y)
  PL_WW_v <- mean(apply(PL_WW_V, MARGIN = 2, FUN = sum))
  PL_WW_v <- PL_WW_v + lambda/2*sum(diag(W%*%K%*%t(W)))
  return(PL_WW_v)
}


DL_WW_Loss <- function(A, K, Y, lambda){
  DL_WW_v <- -sum(A[Y==1])/ncol(K) - 1/(2*lambda*ncol(X)^2)*sum(diag(A%*%K%*%t(A)))
  return(DL_WW_v)
}


data_generate <- function(n,mu_list, sep =  1,  v = 1.5^2, m){
  y <- sort(sample(seq(1,m), size = n-m, replace = T, prob =  1/rep(m,m)))
  n_list <- sapply(seq(1:m), function(t){sum(y==t)+1})
  X <- apply(cbind(mu_list, n_list), MARGIN = 1, FUN = function(mu){matrix(mvrnorm(mu[3], sep*mu[1:2], diag(v,nrow=2)), nrow = mu[3])})
  X <- do.call(rbind,X)
  y <- unlist(sapply(seq(1:m), function(t){rep(t,n_list[t])}))
  return(list(X = X,y = y))
}

set.seed(153)


data_generate <- function(n,mu_list, sep =  1,  v = 1.5^2, m){
  y <- sort(sample(seq(1,m), size = n-m, replace = T, prob =  1/rep(m,m)))
  n_list <- sapply(seq(1:m), function(t){sum(y==t)+1})
  X <- apply(cbind(mu_list, n_list), MARGIN = 1, FUN = function(mu){matrix(mvrnorm(mu[3], sep*mu[1:2], diag(v,nrow=2)), nrow = mu[3])})
  X <- do.call(rbind,X)
  y <- unlist(sapply(seq(1:m), function(t){rep(t,n_list[t])}))
  return(list(X = X,y = y))
}


p <- 2
m <- 4
mu_list <- matrix(rnorm(2*m, mean =0, sd = 10), m)
data <- data_generate(200, mu_list, sep = 3, v = 4.5^2, m)
X <- t(data$X)
y <- data$y

m <- length(unique(y))

n <- ncol(X)
Y <- t(sapply(sort(unique(y)), function(id){as.numeric(y==id)}))
start_time = Sys.time()
W <- matrix(0, m, n)
A <- matrix(0, m, n)
# A <- matrix(rnorm(m*n), m, n)
# K <- t(X)%*%X
# s <- t(W)%*%X
kernel <- rbfdot(sigma = 2^-4)
K <- as.matrix(kernelMatrix(kernel, t(X))) 

lambda <- 1e-2


bt <- max_iter_num <- 10000
P_obj <- rep(0, max_iter_num)
D_obj <- rep(0, length(P_obj))

P_obj[1] <- (m-1)/m
D_obj[1] <- 0


############## Cyclic Update###################

for(t in 2:max_iter_num){
  perm_idx <- sample(n, replace = F)
  j=1
  for(j in perm_idx){
    
    ##Option 6
    q <- c(A%*%K[,j] - K[j,j]*A[,j])
    # q2 <- Y[,j] - (lambda*n*Y[,j] + q)/K[j,j]
    q2 <- -(lambda*n*Y[,j] + q)/K[j,j]
    
    values.cand <- c(Weights%*%sort(b, decreasing = T) - 1/seq(1,m))
    x_res6 <- b - max(values.cand)
    x_res6[x_res6 <= 0] <- 0
    
    A[y[j],j] <- x_res6[y[j]] - 1
    A[-y[j],j] <- x_res6[-y[j]]
  }
  
  W <- -A/lambda/n
  
  P_obj[t] <- PL_WW_Loss(W,K,Y,lambda)
  D_obj[t] <- DL_WW_Loss(A,K,Y,lambda)
  cat("P:", P_obj[t],"D:",  D_obj[t],"Gap:",  P_obj[t] - D_obj[t],"\n")
  
  
  if(P_obj[t]-D_obj[t] < 1e-2){
    bt <- t
    break
  }
}

