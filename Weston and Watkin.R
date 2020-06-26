library(CVXR)
library(MASS)
library(R.matlab)
rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
set.seed(1322)
par(mfrow=c(3,2))
p <- 2
m <- 3
mu1 <- c(-5,5)
mu2 <- c(5,-5)
mu3 <- c(-5,-5)
X1 <- matrix(mvrnorm(1,mu1,matrix(c(5,0,0,5),nrow=2)), nrow =1)
X2 <- matrix(mvrnorm(1,mu2,matrix(c(5,0,0,5),nrow=2)), nrow =1)
X3 <- matrix(mvrnorm(1,mu3,matrix(c(5,0,0,5),nrow=2)), nrow =1)
X <- rbind(X1,X2,X3)

C=1
n <- nrow(X)


y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
Y <- sapply(unique(y), function(id){as.numeric(y==id)})

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "Oracle")

points(X2, col='green')

points(X3, col='blue')

writeMat(con="./matlab code/WW-3class.mat", X =X, Y = y, Y_mat = Y, p=p, m=m, C = C)



w1 <- Variable(p)
w2 <- Variable(p)
w3 <- Variable(p)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 2)
slack2 <- Variable(rows = nrow(X2), cols = 2)
slack3 <- Variable(rows = nrow(X3), cols = 2)

# 
# g <- function(alpha){
#   max(max(alpha, 1)-1, sum_largest(alpha, 2)/2 - 1/2 , sum(alpha)/3-1/3)
# }
# 
# gslack1 <- lapply(seq(1,nrow(X1)), FUN = function(i){g(slack1[i,])})
# gslack2 <- lapply(seq(1,nrow(X2)), FUN = function(i){g(slack2[i,])})
# gslack3 <- lapply(seq(1,nrow(X3)), FUN = function(i){g(slack3[i,])})


C <- 1
objective <- Minimize(sum_squares(vstack(w1,w2,w3)) +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))

constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,1],
                    X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,2],
                    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                    X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,2],
                    X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
                    X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

WW <- Problem(objective, constraints)
WW_primary_time <- system.time(CVXR_WW <- solve(WW, solver = "MOSEK"))


CVXR_WW_w1 <- CVXR_WW$getValue(w1)
CVXR_WW_b1 <- CVXR_WW$getValue(b1)
CVXR_WW_w2 <- CVXR_WW$getValue(w2)
CVXR_WW_b2 <- CVXR_WW$getValue(b2)
CVXR_WW_w3 <- CVXR_WW$getValue(w3)
CVXR_WW_b3 <- CVXR_WW$getValue(b3)

beta1 <- c(CVXR_WW_w1,CVXR_WW_b1)
beta2 <- c(CVXR_WW_w2,CVXR_WW_b2)
beta3 <- c(CVXR_WW_w3,CVXR_WW_b3)
CVXR_WW_primary_beta <- rbind(beta1,beta2,beta3)
CVXR_WW_primary_beta
WW_pri_opt(X,y,C)

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "WW")

points(X2, col='green')

points(X3, col='blue')


# abline(a = -beta1[3]/beta1[2], b = -beta1[1]/beta1[2])
# abline(a = -beta2[3]/beta2[2], b = -beta2[1]/beta2[2])
# abline(a = -beta3[3]/beta3[2], b = -beta3[1]/beta3[2])


WW_beta1_beta2 <- beta1 - beta2
WW_beta1_beta3 <- beta1 - beta3
WW_beta2_beta3 <- beta2 - beta3


abline(a = -WW_beta1_beta2[3]/WW_beta1_beta2[2], b = -WW_beta1_beta2[1]/WW_beta1_beta2[2], lty =2, col = "red")
abline(a = -WW_beta1_beta3[3]/WW_beta1_beta3[2], b = -WW_beta1_beta3[1]/WW_beta1_beta3[2], lty =2, col = "red")
abline(a = -WW_beta2_beta3[3]/WW_beta2_beta3[2], b = -WW_beta2_beta3[1]/WW_beta2_beta3[2], lty =2)
# 
# ###############Cross comparison############################
# source("CS.R")
# 
# source("Duchi.R")
# 
# source("Modified Duchi.R")
# 
# source("MSVM8.R")
###########################################

#######################4 class######
set.seed(14)
X1 <- mvrnorm(20,c(0,2),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(20,c(2,0),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(20,c(-2,0),matrix(c(5,0,0,5),nrow=2))
X4 <- mvrnorm(20,c(0,-2),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3,X4)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X4)))
Y <- sapply(unique(y), function(id){as.numeric(y==id)})


w1 <- Variable(p)
w2 <- Variable(p)
w3 <- Variable(p)
w4 <- Variable(p)

b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)
b4 <- Variable(1)

slack1 <- Variable(rows = nrow(X1), cols = 3)
slack2 <- Variable(rows = nrow(X2), cols = 3)
slack3 <- Variable(rows = nrow(X3), cols = 3)
slack4 <- Variable(rows = nrow(X4), cols = 3)

objective <- Minimize(sum_squares(vstack(w1,w2,w3,w4))/2 +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3) + sum_entries(slack4)))

constraints <- list(w1+w2+w3+w4 == 0, b1+b2+b3+b4 ==0,
                    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,1],
                    X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,2],
                    X1 %*% (w1 - w4) + b1 - b4 >= 1-slack1[,3],
                    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                    X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,2],
                    X2 %*% (w2 - w4) + b2 - b4 >= 1-slack2[,3],
                    X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
                    X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
                    X3 %*% (w3 - w4) + b3 - b4 >= 1-slack3[,3],
                    X4 %*% (w4 - w1) + b4 - b1 >= 1-slack4[,1],
                    X4 %*% (w4 - w2) + b4 - b2 >= 1-slack4[,2],
                    X4 %*% (w4 - w3) + b4 - b3 >= 1-slack4[,3],
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0,
                    slack4 >=0)

WW <- Problem(objective, constraints)
WW_primary_time <- system.time(CVXR_WW <- solve(WW, solver = "MOSEK"))

CVXR_WW_w1 <- CVXR_WW$getValue(w1)
CVXR_WW_b1 <- CVXR_WW$getValue(b1)
CVXR_WW_w2 <- CVXR_WW$getValue(w2)
CVXR_WW_b2 <- CVXR_WW$getValue(b2)
CVXR_WW_w3 <- CVXR_WW$getValue(w3)
CVXR_WW_b3 <- CVXR_WW$getValue(b3)
CVXR_WW_w4 <- CVXR_WW$getValue(w4)
CVXR_WW_b4 <- CVXR_WW$getValue(b4)

beta1 <- c(CVXR_WW_w1,CVXR_WW_b1)
beta2 <- c(CVXR_WW_w2,CVXR_WW_b2)
beta3 <- c(CVXR_WW_w3,CVXR_WW_b3)
beta4 <- c(CVXR_WW_w4,CVXR_WW_b4)
CVXR_WW_primary_beta <- rbind(beta1,beta2,beta3,beta4)
WW_pri_opt(X,y,C)
CVXR_WW_primary_beta

writeMat(con="WW-4class.mat", X =X, Y = y, Y_mat = Y, p=p, m=m, C = C)


LLW_pri_opt(X,y,C)
###########Verifying the objectives ##################################
# alpha_sol
# 
# 
# res1 <- 0
# for(k in 1:m){
#   for( i in 1: n){
#     for( j in 1:n)
#     {
#       res1 <- res1 + Delta[i,k]*Delta[j,k]*alpha_sol[i,k]*alpha_sol[j,k]*t(X[i,])%*%X[j,]
#     }
#   }
# }
# res1
# sum((t(alpha_sol*Delta)%*%X)^2)
# 
# 
# res2 <- 0
# for(k in 1:m){
#   for( i in 1: n){
#     for( j in 1:n){
#       for(l in 1:m){
#         for(s in 1:m)  {
#           res2 <- res2 + Delta[i,l]*Delta[j,s]*Y[i,k]*Y[j,k]*alpha_sol[i,l]*alpha_sol[j,s]*t(X[i,])%*%X[j,]
#           
#         }
#       }
#     }
#   }
# }
# res2
# sum((t(rowSums(alpha_sol*Delta)%*%matrix(1,nrow=1,ncol = m)*Y)%*%X)^2) 
# 
# rowSums(alpha_sol*Delta)%*%matrix(1,nrow=1,ncol = m)*Y
# 
# 
# res3 <- 0
# for(k in 1:m){
#   for( i in 1: n){
#     for( j in 1:n){
#       for(l in 1:m){
#           res3 <- res3 + Delta[i,k]*Delta[j,l]*Y[j,k]*alpha_sol[i,k]*alpha_sol[j,l]*t(X[i,])%*%X[j,]
#        }
#     }
#   }
# }
# res3
# 
# -1/2*res1 + res3 - 1/2*res2
# 
# -sum((t(rowSums(alpha_sol)%*%matrix(1,nrow=1,ncol = m)*Y - alpha_sol)%*%X)^2)/2

##############################################################

n <- nrow(X)
m <- length(unique(y))
p <- ncol(X)

Y <- sapply(unique(y), function(id){as.numeric(y==id)})

Delta <- matrix(1, nrow  = n, ncol = m) - Y
alpha <- Variable(rows = n, cols = m)
beta <- Variable(rows = n, cols = m)

# objective <- Maximize(-1/2*sum_entries((t(alpha*Delta)%*%X)^2) - 1/2*sum_squares(t(sum_entries(alpha*Delta, axis = 1)%*%matrix(1,nrow=1,ncol = m)*Y)%*%X) 
#                       +  sum_entries(alpha))

objective <- Maximize(-sum_squares(t(sum_entries(alpha, axis=1)%*%matrix(1,nrow=1,ncol = m)*Y - alpha)%*%X)/2 +  sum_entries(alpha))


# C = 1
constraints <- list(sum_entries(Y*(sum_entries(Delta*alpha, axis = 1)%*%matrix(1,nrow=1,ncol = m)), axis = 2) == sum_entries(Delta*alpha, axis = 2),
                    sum_entries(alpha*Y, axis = 1) ==0,
                    alpha>=0, 
                    alpha<=C) 


WW_dual <- Problem(objective, constraints)


WW_dual_time <- system.time(WW_dual <- solve(WW_dual, solver = "MOSEK"))

alpha_sol <- WW_dual$getValue(alpha)


solve_w_and_b <- function(X, Y, alpha_sol, C){
  # 
  # sum((round(alpha_sol,6)>0 & round(tau_sum - alpha_sol,6) >0))
  zero_tol <- 4
  w <- t(rowSums(alpha_sol)%*%matrix(1,nrow=1,ncol = m)*Y - alpha_sol)%*%X
  lpred_dif <- X%*%t(w) - (Y*X%*%t(w))%*%matrix(1,nrow = m, ncol=m)

  b_system_mat <- rep(1,m)
  b_system_v <- 0
  m = ncol(Y)
  for(i in 1:m){
    for( j in (i+1):m)
    {
      if(i==m){
        break
      }
      lpred_dif_b <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(alpha_sol,zero_tol) < C) * (y==i)%*%matrix(1,nrow=1, ncol = m)
      lpred_dif_b[lpred_dif_b == 0] <- NA
      b_tmp <- colMeans(lpred_dif_b, na.rm = T)[j]
      
      lpred_dif_b_inv <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(alpha_sol,zero_tol) < C) * (y==j)%*%matrix(1,nrow=1, ncol = m)
      lpred_dif_b_inv[lpred_dif_b_inv == 0] <- NA
      b_tmp <- c(b_tmp ,-colMeans(lpred_dif_b_inv, na.rm = T)[i])
      
      # print(i)
      # print(j)
      # print(mean(b_tmp,na.rm = T))
      if(is.na(mean(b_tmp,na.rm = T))!=T){
        tmp <- rep(0,m)
        tmp[i] <- 1
        tmp[j] <- -1
        b_system_mat <- rbind(b_system_mat, tmp)
        b_system_v <- c(b_system_v, mean(b_tmp,na.rm = T))
      }
    }
  }
  
  b <- ginv(b_system_mat)%*%b_system_v
  return(cbind(w,b))
}
# 
solve_w_and_b(X,Y,alpha_sol, C)
CVXR_WW_primary_beta
