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
X1 <- mvrnorm(20,c(0,4),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(20,c(4,-1),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(20,c(-4,-1),matrix(c(5,0,0,5),nrow=2))
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
constraints <- list(X1 %*% w1 + b1 >= 1-xi1,
                    slack1[,1] == 1,
                    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
                    X2 %*% w2 + b2 >= 1-xi2,
                    slack2[,2] == 1,
                    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                    X3 %*% w1 + b1 <= -1 + slack3[,1],
                    X3 %*% w2 + b2 <= -1 + slack3[,2],
                    xi1 >= max_entries(slack1, axis = 1)-1,
                    xi2 >= max_entries(slack2, axis = 1)-1,
                    xi1 >= sum_entries(slack1, axis = 1)/2 - 1/2,
                    xi2 >= sum_entries(slack2, axis = 1)/2 - 1/2,
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

New1 <- Problem(objective, constraints)
CVXR_New1 <- solve(New1, solver = "MOSEK")
# CVXR_New1$getValue(slack3)

beta1 <- c(CVXR_New1$getValue(w1),CVXR_New1$getValue(b1))
beta2 <- c(CVXR_New1$getValue(w2),CVXR_New1$getValue(b2))

New1_primary_beta <- rbind(beta1,beta2)

New1_primary_beta

New1_beta1_beta2 <- beta1 - beta2
New1_beta1_beta3 <- beta1 - 0
New1_beta2_beta3 <- beta2 - 0

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "New1")
points(X2, col='green')
points(X3, col='blue')

abline(a = -New1_beta1_beta2[3]/New1_beta1_beta2[2], b = -New1_beta1_beta2[1]/New1_beta1_beta2[2], lty =2, col = "red")
abline(a = -New1_beta1_beta3[3]/New1_beta1_beta3[2], b = -New1_beta1_beta3[1]/New1_beta1_beta3[2], lty =2, col = "red")
abline(a = -New1_beta2_beta3[3]/New1_beta2_beta3[2], b = -New1_beta2_beta3[1]/New1_beta2_beta3[2], lty =2)


############################################# 
set.seed(1234)
p <- 2
m <- 4
X1 <- mvrnorm(20,c(0,4),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(20,c(4,-1),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(20,c(-4,-1),matrix(c(5,0,0,5),nrow=2))
X4 <- mvrnorm(20,c(0,-4),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3,X4)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X4)))
C <- 1
n <- nrow(X)
class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
writeMat(con="./matlab code/New1-4class.mat", X1 = X1, X2 = X2, X3 = X3, X4 = X4, p=p, m=m, C = C)

############################################# 


class_idx <- sort(unique(y))
m = length(class_idx)
X_p1 <- X[y!=class_idx[m],,drop = F]
X_p2 <- X[y==class_idx[m],,drop = F]
n_p1 <- nrow(X_p1)
n_p2 <- nrow(X_p2)
Y_p1 <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y[1:n_p1]==id)})

w <- Variable(rows = m-1, cols = p)
b <- Variable(m-1)

slack_p1 <- Variable(rows = n_p1, cols = m-1)
slack_p2 <- Variable(rows = n_p2, cols = m-1)
epsilon_p1 <- Variable(n_p1)
t <- Variable(n_p1*(m-1))
u <- Variable(rows = n_p1*(m-1), cols = m-1)

objective <- Minimize(sum_squares(w)/2 +  C*(sum(epsilon_p1)+sum_entries(slack_p2)))
constraints <- list(((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1)%*%matrix(1, nrow = m-1, ncol = m-1) - (X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b)) >= 1-slack_p1,
                    vec(t(epsilon_p1%*%matrix(1,ncol=m-1))) >= t + sum_entries(rep(1/seq(1,m-1),n_p1)%*%matrix(1,ncol = m-1)*u,axis=1) - rep(1/seq(1,m-1),n_p1),
                    t%*%matrix(1,nrow=1,ncol=m-1) + u >=  (diag(1,n_p1)%x%matrix(1,nrow=m-1))%*%slack_p1,
                    u>=0,
                    sum_entries((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1, axis = 1) >= 1-epsilon_p1,
                    X_p2%*%t(w) + matrix(1,nrow=n_p2,ncol =1)%*%t(b) <= -1 + slack_p2,
                    slack_p1 >=0,
                    slack_p2 >=0)

New1 <- Problem(objective, constraints)
CVXR_New1 <- solve(New1, solver = "MOSEK")

cbind(CVXR_New1$getValue(w),CVXR_New1$getValue(b))
New1_primary_beta




###########################
par(mfrow=c(1,1))
X1 <- matrix(c(2,1),ncol=2)
X2 <- matrix(c(-2,-1),ncol=2)
X3 <- matrix(c(-2,1),ncol=2)
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
C <- 100


New1_primary_beta <- New1_pri_opt(X,y,C)
beta1 <- round(New1_primary_beta[1,],2)
beta2 <- round(New1_primary_beta[2,],2)
New1_beta1_beta2 <- beta1 - beta2
New1_beta1_beta3 <- beta1 - 0
New1_beta2_beta3 <- beta2 - 0

plot(X1, col='red', xlim = c(-4,4), ylim=c(-4,4), xlab = "X1", ylab = "X2", main = "New1")
points(X2, col='green')
points(X3, col='blue')

abline(a = -New1_beta1_beta2[3]/New1_beta1_beta2[2], b = -New1_beta1_beta2[1]/New1_beta1_beta2[2], lty =2, col = "red")
abline(a = -New1_beta1_beta3[3]/New1_beta1_beta3[2], b = -New1_beta1_beta3[1]/New1_beta1_beta3[2], lty =2, col = "red")
abline(a = -New1_beta2_beta3[3]/New1_beta2_beta3[2], b = -New1_beta2_beta3[1]/New1_beta2_beta3[2], lty =2)



###########################
par(mfrow=c(1,1))
X1 <- matrix(c(-2,1),ncol=2)
X2 <- matrix(c(2,1),ncol=2)
X3 <- matrix(c(-2,-1),ncol=2)
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
C <- 100000

Duchi_primary_beta <- Duchi_pri_opt(X,y,C)
beta = Duchi_primary_beta

plot_decision_boundary <- function(X, y, beta, title = NULL, np_reslution = 500){
  
  if(nrow(beta) == length(unique(y))){
    X_dat <- as.data.frame(X)
    nd.x = seq(from = floor(min(X_dat$V1))-1, to =  ceiling(max(X_dat$V1))+1, length.out = np_reslution)
    nd.y = seq(from = floor(min(X_dat$V1))-1, to =  ceiling(max(X_dat$V2))+1, length.out = np_reslution)
    nd = expand.grid(Var1 = nd.x, Var2 = nd.y)
    prd = apply(as.matrix(cbind(nd,1))%*%t(beta), MARGIN = 1, FUN = which.max)
  }else{
    X_dat <- as.data.frame(X)
    nd.x = seq(from = floor(min(X_dat$V1))-1, to =  ceiling(max(X_dat$V1))+1, length.out = np_reslution)
    nd.y = seq(from = floor(min(X_dat$V1))-1, to =  ceiling(max(X_dat$V2))+1, length.out = np_reslution)
    nd = expand.grid(Var1 = nd.x, Var2 = nd.y)
    prd_p = as.matrix(cbind(nd,1))%*%t(beta)
    prd_p =  cbind(prd_p, 1 - apply(as.matrix(cbind(prd_p[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(prd_p[,2]+1,0)), MARGIN = 1, FUN = max))
    prd = apply(prd_p, MARGIN = 1, FUN = which.max)
  }
  par(mar=c(5.1, 4.1, 4.1, 7), xpd=TRUE)
  plot(X_dat$V1, X_dat$V2, col = as.factor(y), 
       ylim=c( floor(min(X_dat$V2))-1, ceiling(max(X_dat$V2))+1), xlim=c( floor(min(X_dat$V1))-1, ceiling(max(X_dat$V1))+1), xlab ="X", ylab = "Y", main = title)
  
  contour(x = nd.x, y = nd.y, z = matrix(prd, nrow = np_reslution, ncol = np_reslution), 
          levels = unique(y), add = TRUE, drawlabels = FALSE)
  
  legend("topright", inset=c(-0.3,0),legend = sapply(unique(y), function(x){paste("Class ",x)}), col= unique(y), pch = 1)
}


plot_decision_boundary(X,y,Duchi_primary_beta, title = "Duchi")
plot_decision_boundary(X,y,New1_primary_beta, title = "New1")
