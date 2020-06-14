library(CVXR)
library(MASS)
library(R.matlab)
set.seed(12223)
rm(list = ls())

setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")
p <- 2
m <- 3
X1 <- mvrnorm(4,c(0,5),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(4,c(5,0),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(4,c(0,0),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
C <- 1
n <- nrow(X)
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})

# 
# p <- 900
# m <- 3
# v <- 200
# X1 <- mvrnorm(300, rep(1,p), diag(v,p))
# X2 <- mvrnorm(300, rep(-1,p), diag(v,p))
# X3 <- mvrnorm(300, c(rep(1,p/2),rep(-1,p/2)), diag(v,p))
# y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
# X4 <- mvrnorm(300, c(rep(-1,p/2),rep(1,p/2)), diag(v,p))
# y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X3)))
C <- 1
writeMat(con="Duchi-3class.mat", X1 = X1, X2 = X2, X3 = X3, X=X, Y_mat = Y, p=p, m=m, C = C)
  
w1 <- Variable(p)
w2 <- Variable(p)
w3 <- Variable(p)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 3)
slack2 <- Variable(rows = nrow(X2), cols = 3)
slack3 <- Variable(rows = nrow(X3), cols = 3)
xi1 <- Variable(rows = nrow(X1))
xi2 <- Variable(rows = nrow(X2))
xi3 <- Variable(rows = nrow(X3))


  # g <- function(alpha){
  #   max(max(alpha)-1, sum_largest(alpha, 2)/2 - 1/2 , sum(alpha)/3-1/3)
  # }
  #
  # gslack1 <- lapply(seq(1,nrow(X1)), FUN = function(i){g(slack1[i,])})
  # gslack2 <- lapply(seq(1,nrow(X2)), FUN = function(i){g(slack2[i,])})
  # gslack3 <- lapply(seq(1,nrow(X3)), FUN = function(i){g(slack3[i,])})


objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum(xi1)+sum(xi2)+sum(xi3)))
constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                      slack1[,1] == 1,
                      X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
                      X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,3],
                      slack2[,2] == 1,
                      X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                      X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,3],
                      slack3[,3] == 1,
                      X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
                      X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
                      xi1 >= max_entries(slack1, axis = 1)-1,
                      xi2 >= max_entries(slack2, axis = 1)-1,
                      xi3 >= max_entries(slack3, axis = 1)-1,
                      xi1 >= (sum_entries(slack1, axis = 1) - min_entries(slack1, axis = 1))/2-1/2,
                      xi2 >= (sum_entries(slack2, axis = 1) - min_entries(slack2, axis = 1))/2-1/2,
                      xi3 >= (sum_entries(slack3, axis = 1) - min_entries(slack3, axis = 1))/2-1/2,
                      xi1 >= sum_entries(slack1, axis = 1)/3 - 1/3,
                      xi2 >= sum_entries(slack2, axis = 1)/3 - 1/3,
                      xi3 >= sum_entries(slack3, axis = 1)/3 - 1/3,
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)

Duchi <- Problem(objective, constraints)
CVXR_Duchi <- solve(Duchi, solver = "MOSEK")


######################################################2


CVXR_Duchi_w1 <- CVXR_Duchi$getValue(w1)
CVXR_Duchi_b1 <- CVXR_Duchi$getValue(b1)
CVXR_Duchi_w2 <- CVXR_Duchi$getValue(w2)
CVXR_Duchi_b2 <- CVXR_Duchi$getValue(b2)
CVXR_Duchi_w3 <- CVXR_Duchi$getValue(w3)
CVXR_Duchi_b3 <- CVXR_Duchi$getValue(b3)

beta1 <- c(CVXR_Duchi_w1,CVXR_Duchi_b1)
beta2 <- c(CVXR_Duchi_w2,CVXR_Duchi_b2)
beta3 <- c(CVXR_Duchi_w3,CVXR_Duchi_b3)

Duchi_primary_beta <- rbind(beta1,beta2,beta3)

slack <- Variable(rows = n, cols = 3)
xi <- Variable(rows = n, cols = 1)
cn1 <-nrow(X1)
cn2 <-nrow(X1) + nrow(X2)
cn3 <-nrow(X1) + nrow(X2) + nrow(X3)

objective2 <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum(xi)))
constraints2 <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                    slack[1:cn1,1] == 1,
                    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack[1:cn1,2],
                    X1 %*% (w1 - w3) + b1 - b3 >= 1-slack[1:cn1,3],
                    slack[(cn1+1):cn2,2] == 1,
                    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack[(cn1+1):cn2,1],
                    X2 %*% (w2 - w3) + b2 - b3 >= 1-slack[(cn1+1):cn2,3],
                    slack[(cn2+1):cn3,3] == 1,
                    X3 %*% (w3 - w1) + b3 - b1 >= 1-slack[(cn2+1):cn3,1],
                    X3 %*% (w3 - w2) + b3 - b2 >= 1-slack[(cn2+1):cn3,2],
                    xi >= max_entries(slack, axis = 1)-1,
                    xi >= (sum_entries(slack, axis = 1) - min_entries(slack, axis = 1))/2-1/2,
                    xi >= sum_entries(slack, axis = 1)/3 - 1/3,
                    slack >=0)



Duchi2 <- Problem(objective2, constraints2)
CVXR_Duchi2 <- solve(Duchi2, solver = "MOSEK")



CVXR_Duchi2_w1 <- CVXR_Duchi2$getValue(w1)
CVXR_Duchi2_b1 <- CVXR_Duchi2$getValue(b1)
CVXR_Duchi2_w2 <- CVXR_Duchi2$getValue(w2)
CVXR_Duchi2_b2 <- CVXR_Duchi2$getValue(b2)
CVXR_Duchi2_w3 <- CVXR_Duchi2$getValue(w3)
CVXR_Duchi2_b3 <- CVXR_Duchi2$getValue(b3)

beta1 <- c(CVXR_Duchi2_w1,CVXR_Duchi2_b1)
beta2 <- c(CVXR_Duchi2_w2,CVXR_Duchi2_b2)
beta3 <- c(CVXR_Duchi2_w3,CVXR_Duchi2_b3)

Duchi_primary_beta2 <- rbind(beta1,beta2,beta3)

Duchi_primary_beta
Duchi_primary_beta2

w <- t(cbind(CVXR_Duchi2_w1, CVXR_Duchi2_w2, CVXR_Duchi2_w3))
b <- c(CVXR_Duchi2_b1, CVXR_Duchi2_b2, CVXR_Duchi2_b3)



######################################################3
# ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))
# 
# X1 %*% (CVXR_Duchi_w1 - CVXR_Duchi_w2) + CVXR_Duchi_b1 - CVXR_Duchi_b2


slack <- Variable(rows = n, cols = 3)
xi <- Variable(rows = n, cols = 1)
t <- Variable(n*m)
u <- Variable(rows = n*m, cols = m)

w <- Variable(rows = m, cols = p)
b <- Variable(rows = m, cols = 1)



objective3 <- Minimize(sum_squares(w)/2 +  C*sum(xi))
constraints3 <- list( sum_entries(b) == 0,
                     sum_entries(w, axis = 2) ==0,
                     ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                    xi >= max_entries(slack, axis = 1)-1,
                    xi >= (sum_entries(slack, axis = 1) - min_entries(slack, axis = 1))/2-1/2,
                    xi >= sum_entries(slack, axis = 1)/3 - 1/3,
                    # vec(t(xi%*%matrix(1,ncol=m))) >= t - sum_entries(rep(1/seq(1,m),n)%*%matrix(1,ncol = m)*u,axis=1) - rep(1/seq(1,m),n),
                    # t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m))%*%slack,
                    # u>=0,
                    slack >=0)



Duchi3 <- Problem(objective3, constraints3)
CVXR_Duchi3 <- solve(Duchi3, solver = "MOSEK")
Duchi_primary_beta3 <- cbind(CVXR_Duchi3$getValue(w),CVXR_Duchi3$getValue(b))
  
  
  Duchi_primary_beta
  Duchi_primary_beta2
  Duchi_primary_beta3

######################################################4
# ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))
# 
# X1 %*% (CVXR_Duchi_w1 - CVXR_Duchi_w2) + CVXR_Duchi_b1 - CVXR_Duchi_b2


slack <- Variable(rows = n, cols = 3)
xi <- Variable(rows = n, cols = 1)
t <- Variable(n*m)
u <- Variable(rows = n*m, cols = m)

w <- Variable(rows = m, cols = p)
b <- Variable(rows = m, cols = 1)



objective4 <- Minimize(sum_squares(w)/2 +  C*sum(xi))
constraints4 <- list( sum_entries(b) == 0,
                      sum_entries(w, axis = 2) ==0,
                      ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                      # xi >= max_entries(slack, axis = 1)-1,
                      # xi >= (sum_entries(slack, axis = 1) - min_entries(slack, axis = 1))/2-1/2,
                      # xi >= sum_entries(slack, axis = 1)/3 - 1/3,
                      vec(t(xi%*%matrix(1,ncol=m))) >= t - sum_entries(rep(1/seq(1,m),n)%*%matrix(1,ncol = m)*u,axis=1) - rep(1/seq(1,m),n),
                      t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m))%*%slack,
                      u>=0,
                      slack >=0)



Duchi4 <- Problem(objective4, constraints4)
CVXR_Duchi4 <- solve(Duchi4, solver = "MOSEK")
Duchi_primary_beta4 <- cbind(CVXR_Duchi4$getValue(w),CVXR_Duchi4$getValue(b))


Duchi_primary_beta
Duchi_primary_beta2
Duchi_primary_beta3
Duchi_primary_beta4


  
#######################

w <- t(cbind(CVXR_Duchi_w1,CVXR_Duchi_w2,CVXR_Duchi_w3))
b <- c(CVXR_Duchi_b1,CVXR_Duchi_b2,CVXR_Duchi_b3)

class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
n <- nrow(X)
p <- ncol(X)
m <- length(class_idx)
((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))


w <- Variable(rows = m, cols = p)
b <- Variable(m)
slack <- Variable(rows = n, cols = m)
epsilon <- Variable(n)
t <- Variable(n*m)
u <- Variable(rows = n*m, cols = m)



objective <- Minimize(sum_squares(w)/2 +  C*sum(epsilon))
constraints <- list( sum_entries(b) == 0,
                     sum_entries(w, axis = 2) ==0,
                     ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                     vec(t(epsilon%*%matrix(1,ncol=m))) >= t - sum_entries(rep(1/seq(1,m),n)%*%matrix(1,ncol = m)*u,axis=1) - rep(1/seq(1,m),n),
                     # (diag(1,n)%x%matrix(1,nrow=m))%*%epsilon >=  rep(1/seq(1,m),n)*(rep(seq(1,m),n)*t - sum_entries(u,axis=1)) - rep(1/seq(1,m),n),
                     t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m))%*%slack,
                     slack*Y == 1,
                     u>=0,
                     # epsilon >= max_entries((1-Y)*slack, axis = 1),
                     slack >=0)

Duchi2 <- Problem(objective, constraints)
CVXR_Duchi2 <- solve(Duchi2, solver = "MOSEK")
# CS_pri_opt(X,y,C)
CVXR_Duchi2$getValue(w)
CVXR_Duchi2$getValue(b)

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "Duchi")

points(X2, col='green')

points(X3, col='blue')


# abline(a = -beta1[3]/beta1[2], b = -beta1[1]/beta1[2])
# abline(a = -beta2[3]/beta2[2], b = -beta2[1]/beta2[2])
# abline(a = -beta3[3]/beta3[2], b = -beta3[1]/beta3[2])


Duchi_beta1_beta2 <- beta1 - beta2
Duchi_beta1_beta3 <- beta1 - beta3
Duchi_beta2_beta3 <- beta2 - beta3


abline(a = -Duchi_beta1_beta2[3]/Duchi_beta1_beta2[2], b = -Duchi_beta1_beta2[1]/Duchi_beta1_beta2[2], lty =2, col = "red")
abline(a = -Duchi_beta1_beta3[3]/Duchi_beta1_beta3[2], b = -Duchi_beta1_beta3[1]/Duchi_beta1_beta3[2], lty =2, col = "red")
abline(a = -Duchi_beta2_beta3[3]/Duchi_beta2_beta3[2], b = -Duchi_beta2_beta3[1]/Duchi_beta2_beta3[2], lty =2)

##########################################
set.seed(1223)
# set.seed(1234)
p <- 2
m <- 4
Sigma <- diag(10,nrow=2)
X1 <- mvrnorm(10,c(0,2),Sigma)
X2 <- mvrnorm(10,c(2,0),Sigma)
X3 <- mvrnorm(10,c(-2,0),Sigma)
X4 <- mvrnorm(10,c(-2,-2),Sigma)
C <- 1
X <- rbind(X1,X2,X3,X4)

y <- c(rep(1,nrow(X1)), rep(2,nrow(X2)), rep(3,nrow(X3)), rep(4,nrow(X4)))
writeMat(con="Duchi-4class.mat", X1 = X1, X2 = X2, X3 = X3, X4 = X4, p=p, m=m, C = C)


################################################################
n <- nrow(X)
m <- length(unique(y))
p <- ncol(X)
Y <- sapply(unique(y), function(id){as.numeric(y==id)})

alpha <- Variable(rows = n, cols = m)
beta <- Variable(rows = n, cols = m)
tau <- Variable(rows = n*m, cols = m)

objective <- Maximize(-sum_squares(t(sum_entries(alpha, axis=1)%*%matrix(1,nrow=1,ncol = m)*Y - alpha)%*%X)/2 + 
                        sum_entries(alpha) - sum_entries(beta*sapply(1:m, function(j){rep(1/j,n)})))

constraints <- list(sum_entries(Y*(sum_entries(alpha, axis = 1)%*%matrix(1,nrow=1,ncol = m)), axis = 2) == sum_entries(alpha, axis = 2),
                    # sum_entries(alpha*Y, axis = 1) ==0,
                    C - sum_entries(beta, axis =1) ==0, 
                    vec(beta*sapply(1:m, function(j){rep(1/j,n)}))%*%matrix(1,nrow=1,ncol=m) - tau >=0, 
                    beta - reshape_expr(sum_entries(tau, axis =1), c(n,m))==0,
                    # sum_entries(reshape_expr(tau, c(n, m*m))[,1:3], axis=1) - alpha[,1] >=0,
                    # sum_entries(reshape_expr(tau, c(n, m*m))[,4:6], axis=1) - alpha[,2] >=0,
                    # sum_entries(reshape_expr(tau, c(n, m*m))[,7:9], axis=1) - alpha[,3] >=0,
                    alpha>=0,
                    beta>=0,
                    tau>=0) 

constraints <- c(constraints,lapply(0:(m-1), function(id){
  sum_entries(reshape_expr(tau, c(n, m*m))[,(id*m+1):((id+1)*m)], axis=1) - alpha[,id+1] >=0}))


Duchi_dual <- Problem(objective, constraints)


Duchi_dual_time <- system.time(Duchi_dual <- solve(Duchi_dual, solver = "MOSEK"))

alpha_sol <- Duchi_dual$getValue(alpha)
tau_sol <- Duchi_dual$getValue(tau)
# i=1;j=4
solve_w_and_b <- function(X, Y, alpha_sol, tau_sol){
  min_rm_zero <- function(x){
    ifelse(sum(x!=0)!=0, min(x[x!=0]),NA)
  }
  max_rm_zero <- function(x){
    ifelse(sum(x!=0)!=0, max(x[x!=0]),NA)
  }
  zero_tol <- 4
  m = ncol(Y)
  w <- t(rowSums(alpha_sol)%*%matrix(1,nrow=1,ncol = m)*Y - alpha_sol)%*%X
  lpred_dif <- X%*%t(w) - (Y*X%*%t(w))%*%matrix(1,nrow = m, ncol=m)
  tau_sum <- do.call(cbind,lapply(0:(m-1), function(id,tau){
    rowSums(tau[,(id*m+1):((id+1)*m)])}, tau = matrix(c(tau_sol), n, m*m, byrow = F))) 
  b_system_mat <- rep(1,m)
  b_system_v <- 0
  b_system_mat2 <- NULL
  b_system_v2 <- NULL
  for(i in 1:m){
    for( j in (i+1):m)
    {
        if(i==m){
          break
        }
        lpred_dif_b <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(tau_sum - alpha_sol,zero_tol) >0) * (y==i)%*%matrix(1,nrow=1, ncol = m)
        lpred_dif_b[lpred_dif_b == 0] <- NA
        b_tmp <- colMeans(lpred_dif_b, na.rm = T)[j]
        
        lpred_dif_b_inv <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(tau_sum - alpha_sol,zero_tol) >0) * (y==j)%*%matrix(1,nrow=1, ncol = m)
        lpred_dif_b_inv[lpred_dif_b_inv == 0] <- NA
        b_tmp <- c(b_tmp ,-colMeans(lpred_dif_b_inv, na.rm = T)[i])

        if(is.na(mean(b_tmp,na.rm = T))==T){
          b_tmp_upper <- min(apply((1+lpred_dif)*((round(alpha_sol,zero_tol) == round(tau_sum,zero_tol)) & round(tau_sum, zero_tol) >0) * (y==i)%*%matrix(1,nrow=1, ncol = m), MARGIN = 2, FUN = min_rm_zero)[j],
          -apply((1+lpred_dif)*(round(alpha_sol,zero_tol) ==0 & round(tau_sum, zero_tol) >0) * (y==j)%*%matrix(1,nrow=1, ncol = m), MARGIN = 2, FUN = max_rm_zero)[i], na.rm = T)
          
          b_tmp_lower <- max(-apply((1+lpred_dif)*((round(alpha_sol,zero_tol) == round(tau_sum,zero_tol)) & round(tau_sum, zero_tol) >0) * (y==j)%*%matrix(1,nrow=1, ncol = m), MARGIN = 2, FUN = min_rm_zero)[i],
                    apply((1+lpred_dif)*(round(alpha_sol,zero_tol) ==0 & round(tau_sum, zero_tol) >0) * (y==i)%*%matrix(1,nrow=1, ncol = m), MARGIN = 2, FUN = max_rm_zero)[j],na.rm = T)
        }
        
        
        if(is.na(mean(b_tmp,na.rm = T))!=T){
          tmp <- rep(0,m)
          tmp[i] <- 1
          tmp[j] <- -1
          b_system_mat <- rbind(b_system_mat, tmp)
          b_system_v <- c(b_system_v, mean(b_tmp,na.rm = T))
        }else if(is.na(mean(c(b_tmp_upper,b_tmp_lower),na.rm = T))!=T){
            tmp <- rep(0,m)
            tmp[i] <- 1
            tmp[j] <- -1
            b_system_mat2 <- rbind(b_system_mat2, tmp)
            b_system_v2 <- c(b_system_v, mean(c(b_tmp_upper,b_tmp_lower),na.rm = T))
        }
    }
  }
  if(sum(svd(b_system_mat)$d>0)==m){
    b <- ginv(b_system_mat)%*%b_system_v
  }else{
    b <- ginv(rbind(b_system_mat,b_system_mat2))%*%c(b_system_v,b_system_v2)
  }
  
  return(cbind(w,b))
}
# 
solve_w_and_b(X,Y,alpha_sol, tau_sol)
# b_system_mat
# b_system_v
# Duchi_pri_opt(X,y,C)

Duchi_dual_opt(X,y,C)



apply((1+lpred_dif)*((round(alpha_sol,zero_tol) == round(tau_sum,zero_tol)) & round(tau_sum, zero_tol) >0) * (y==1)%*%matrix(1,nrow=1, ncol = m), MARGIN = 2, FUN = min_rm_zero)

