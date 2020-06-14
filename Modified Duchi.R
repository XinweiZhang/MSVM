library(CVXR)
library(MASS)
set.seed(136)
rm(list = ls())

p <- 2
m <- 3
X1 <- mvrnorm(10,c(0,5),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(10,c(5,0),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(10,c(0,0),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
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
w1 <- Variable(p)
w2 <- Variable(p)
w3 <- Variable(p)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 2)
slack2 <- Variable(rows = nrow(X2), cols = 2)
slack3 <- Variable(rows = nrow(X3), cols = 2)

xi1 <- Variable(rows = nrow(X1))
xi2 <- Variable(rows = nrow(X2))
xi3 <- Variable(rows = nrow(X3))

C <- 1
objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum(xi1)+sum(xi2)+sum(xi3)))
constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,1],
                    X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,2],
                    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                    X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,2],
                    X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
                    X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
                    xi1 >= max_entries(slack1, axis = 1)/2,
                    xi2 >= max_entries(slack2, axis = 1)/2,
                    xi3 >= max_entries(slack3, axis = 1)/2,
                    xi1 >= sum_entries(slack1, axis = 1)/3,
                    xi2 >= sum_entries(slack2, axis = 1)/3,
                    xi3 >= sum_entries(slack3, axis = 1)/3,
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

M_Duchi <- Problem(objective, constraints)
M_Duchi_primary_time <- system.time(CVXR_M_Duchi <- solve(M_Duchi, solver = "MOSEK"))

CVXR_M_Duchi_w1 <- CVXR_M_Duchi$getValue(w1)
CVXR_M_Duchi_b1 <- CVXR_M_Duchi$getValue(b1)
CVXR_M_Duchi_w2 <- CVXR_M_Duchi$getValue(w2)
CVXR_M_Duchi_b2 <- CVXR_M_Duchi$getValue(b2)
CVXR_M_Duchi_w3 <- CVXR_M_Duchi$getValue(w3)
CVXR_M_Duchi_b3 <- CVXR_M_Duchi$getValue(b3)

beta1 <- c(CVXR_M_Duchi_w1,CVXR_M_Duchi_b1)
beta2 <- c(CVXR_M_Duchi_w2,CVXR_M_Duchi_b2)
beta3 <- c(CVXR_M_Duchi_w3,CVXR_M_Duchi_b3)

CVXR_M_Duchi_primary_beta <- rbind(beta1,beta2,beta3)


class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
n <- nrow(X)
p <- ncol(X)
m <- length(class_idx)

w <- Variable(rows = m, cols = p)
b <- Variable(m)
slack <- Variable(rows = n, cols = m)
epsilon <- Variable(n)
t <- Variable(n*(m-1))
u <- Variable(rows = n*(m-1), cols = m)

dim(t%*%matrix(1,nrow=1,ncol=m-1))

objective <- Minimize(sum_squares(w)/2 +  C*sum(epsilon))
constraints <- list( sum_entries(b) == 0,
                     sum_entries(w, axis = 2) ==0,
                     ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-Y)*(1-slack),
                     # vec(t(epsilon%*%matrix(1,ncol=m))) >= t - rep(1/seq(1,m),n)*sum_entries(u,axis=1) - rep(1/seq(1,m),n),
                     (diag(1,n)%x%matrix(1,nrow=m-1))%*%epsilon >=  rep(seq(1,m-1),n)/rep(seq(2,m),n)*t - 1/rep(seq(2,m),n)*sum_entries(u,axis=1),
                     t%*%matrix(1,nrow=1,ncol=m-1) + u >=  (diag(1,n)%x%matrix(1,nrow=m-1))%*%slack,
                     u>=0,
                     slack >=0)


########plot##############################
plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "Modified Duchi")

points(X2, col='green')

points(X3, col='blue')
# abline(a = -beta1[3]/beta1[2], b = -beta1[1]/beta1[2])
# abline(a = -beta2[3]/beta2[2], b = -beta2[1]/beta2[2])
# abline(a = -beta3[3]/beta3[2], b = -beta3[1]/beta3[2])



M_Duchi_beta1_beta2 <- beta1 - beta2
M_Duchi_beta1_beta3 <- beta1 - beta3
M_Duchi_beta2_beta3 <- beta2 - beta3


abline(a = -M_Duchi_beta1_beta2[3]/M_Duchi_beta1_beta2[2], b = -M_Duchi_beta1_beta2[1]/M_Duchi_beta1_beta2[2], lty =2, col = "red")
abline(a = -M_Duchi_beta1_beta3[3]/M_Duchi_beta1_beta3[2], b = -M_Duchi_beta1_beta3[1]/M_Duchi_beta1_beta3[2], lty =2, col = "red")
abline(a = -M_Duchi_beta2_beta3[3]/M_Duchi_beta2_beta3[2], b = -M_Duchi_beta2_beta3[1]/M_Duchi_beta2_beta3[2], lty =2)


##################### 4 class########################
set.seed(1323)
X1 <- mvrnorm(20,c(0,2),matrix(c(8,0,0,8),nrow=2))
X2 <- mvrnorm(20,c(2,0),matrix(c(8,0,0,8),nrow=2))
X3 <- mvrnorm(20,c(-2,-2),matrix(c(8,0,0,8),nrow=2))
X4 <- mvrnorm(20,c(0,-2),matrix(c(8,0,0,8),nrow=2))
X <- rbind(X1,X2,X3,X4)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X4)))
m <- 4
p <- 2
C <- 1
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})

writeMat(con="MDuchi-4class.mat", X1 = X1, X2 = X2, X3 = X3, X4 = X4, X=X, Y_mat =Y, p=p, m=m, C = C)

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

xi1 <- Variable(rows = nrow(X1))
xi2 <- Variable(rows = nrow(X2))
xi3 <- Variable(rows = nrow(X3))
xi4 <- Variable(rows = nrow(X4))

C <- 1
objective <- Minimize(sum_squares(vstack(w1,w2,w3,w4))/2 +  C*(sum(xi1)+sum(xi2)+sum(xi3)+ sum(xi4)) )
constraints <- list(w1+w2+w3 + w4 == 0, b1+b2+b3 +b4==0,
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
                    xi1 >= max_entries(slack1, axis = 1)/2,
                    xi2 >= max_entries(slack2, axis = 1)/2,
                    xi3 >= max_entries(slack3, axis = 1)/2,
                    xi4 >= max_entries(slack4, axis = 1)/2,
                    xi1 >= (sum_entries(slack1, axis = 1)-min_entries(slack1, axis = 1))/3,
                    xi2 >= (sum_entries(slack2, axis = 1)-min_entries(slack2, axis = 1))/3,
                    xi3 >= (sum_entries(slack3, axis = 1)-min_entries(slack3, axis = 1))/3,
                    xi4 >= (sum_entries(slack4, axis = 1)-min_entries(slack4, axis = 1))/3,
                    xi1 >= sum_entries(slack1, axis = 1)/4,
                    xi2 >= sum_entries(slack2, axis = 1)/4,
                    xi3 >= sum_entries(slack3, axis = 1)/4,
                    xi4 >= sum_entries(slack4, axis = 1)/4,
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0,
                    slack4 >=0)

M_Duchi <- Problem(objective, constraints)
M_Duchi_primary_time <- system.time(CVXR_M_Duchi <- solve(M_Duchi, solver = "MOSEK"))


CVXR_M_Duchi_w1 <- CVXR_M_Duchi$getValue(w1)
CVXR_M_Duchi_b1 <- CVXR_M_Duchi$getValue(b1)
CVXR_M_Duchi_w2 <- CVXR_M_Duchi$getValue(w2)
CVXR_M_Duchi_b2 <- CVXR_M_Duchi$getValue(b2)
CVXR_M_Duchi_w3 <- CVXR_M_Duchi$getValue(w3)
CVXR_M_Duchi_b3 <- CVXR_M_Duchi$getValue(b3)
CVXR_M_Duchi_w4 <- CVXR_M_Duchi$getValue(w4)
CVXR_M_Duchi_b4 <- CVXR_M_Duchi$getValue(b4)

beta1 <- c(CVXR_M_Duchi_w1,CVXR_M_Duchi_b1)
beta2 <- c(CVXR_M_Duchi_w2,CVXR_M_Duchi_b2)
beta3 <- c(CVXR_M_Duchi_w3,CVXR_M_Duchi_b3)
beta4 <- c(CVXR_M_Duchi_w4,CVXR_M_Duchi_b4)

CVXR_M_Duchi_primary_beta <- rbind(beta1,beta2,beta3,beta4)


##################
# ##################### 5 class########################
set.seed(32)
X1 <- mvrnorm(20,c(0,2),matrix(c(8,0,0,8),nrow=2))
X2 <- mvrnorm(20,c(2,0),matrix(c(8,0,0,8),nrow=2))
X3 <- mvrnorm(20,c(-2,-2),matrix(c(8,0,0,8),nrow=2))
X4 <- mvrnorm(20,c(0,-2),matrix(c(8,0,0,8),nrow=2))
X5 <- mvrnorm(20,c(-2,0),matrix(c(8,0,0,8),nrow=2))
X <- rbind(X1,X2,X3,X4,X5)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X4)),rep(5,nrow(X5)))
m <- 5
p <- 2
C <- 1
writeMat(con="MDuchi-5class.mat", X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, p=p, m=m, C = C)


#########################################################################

# X <- rbind(X1,X2,X3,X4)
n <- nrow(X)
m <- length(unique(y))
p <- ncol(X)

Y <- sapply(unique(y), function(id){as.numeric(y==id)})

Delta <- matrix(1, nrow  = n, ncol = m) - Y
alpha <- Variable(rows = n, cols = m)
beta <- Variable(rows = n, cols = m-1)
tau <- Variable(rows = n*(m-1), cols = m)

objective <- Maximize(-sum_squares(t(sum_entries(alpha, axis=1)%*%matrix(1,nrow=1,ncol = m)*Y - alpha)%*%X)/2 +  sum_entries(alpha))
C=1
constraints <- list(sum_entries(Y*(sum_entries(alpha, axis = 1)%*%matrix(1,nrow=1,ncol = m)), axis = 2) == sum_entries(alpha, axis = 2),
                    # sum_entries(Y*(sum_entries(Delta*a  lpha, axis = 1)%*%matrix(1,nrow=1,ncol = m)), axis = 2) == sum_entries(Delta*alpha, axis = 2),
                    sum_entries(alpha*Y, axis = 1) ==0,
                    C - sum_entries(beta, axis =1) ==0, 
                    vec(beta*sapply(1:(m-1), function(j){rep(1/(j+1),n)}))%*%matrix(1,nrow=1,ncol=m) - tau >=0, 
                    beta*sapply(1:(m-1), function(j){rep(j/(j+1),n)}) - reshape_expr(sum_entries(do.call(rbind,lapply(1:(m-1),
                                                 function(x){1-Y}))*tau, axis =1), c(n,m-1))==0,
                    # sum_entries(reshape_expr(tau, c(n, (m-1)*m))[,1:3], axis=1) - alpha[,1] >=0,
                    # sum_entries(reshape_expr(tau, c(n, (m-1)*m))[,4:6], axis=1) - alpha[,2] >=0,
                    # sum_entries(reshape_expr(tau, c(n, (m-1)*m))[,7:9], axis=1) - alpha[,3] >=0,
                    # sum_entries(reshape_expr(tau, c(n, (m-1)*m))[,10:12], axis=1) - alpha[,4] >=0,
                    sum_entries(do.call(rbind,lapply(1:(m-1), function(x){Y}))*tau, axis =1) == 0,
                    alpha>=0,
                    # alpha*Y==0,
                    beta>=0,
                    tau>=0) 
constraints <- c(constraints,lapply(0:(m-1), function(id){sum_entries(reshape_expr(tau, c(n, (m-1)*m))[,(id*(m-1)+1):((id+1)*(m-1))], axis=1) - alpha[,id+1] >=0}))


M_Duchi_dual <- Problem(objective, constraints)


M_Duchi_dual_time <- system.time(M_Duchi_dual <- solve(M_Duchi_dual, solver = "MOSEK"))

alpha_sol <- M_Duchi_dual$getValue(alpha)
tau_sol <- M_Duchi_dual$getValue(tau)

solve_w_and_b <- function(X, Y, alpha_sol, tau_sol){
  zero_tol <- 4
  w <- t(rowSums(alpha_sol)%*%matrix(1,nrow=1,ncol = m)*Y - alpha_sol)%*%X
  lpred_dif <- X%*%t(w) - (Y*X%*%t(w))%*%matrix(1,nrow = m, ncol=m)
  m = ncol(Y)
  tau_sum <- do.call(cbind,lapply(0:(m-1), function(id,tau){
    rowSums(tau[,(id*(m-1)+1):((id+1)*(m-1))])}, tau = matrix(c(tau_sol), n, (m-1)*m, byrow = F))) 
  b_system_mat <- rep(1,m)
  b_system_v <- 0
  
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

solve_w_and_b(X,Y,alpha_sol, tau_sol)
# CVXR_M_Duchi_primary_beta
