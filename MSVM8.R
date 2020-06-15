library(CVXR)
library(MASS)
# set.seed(12453)
par(mfrow=c(2,2))
rm(list = ls())

source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")
p <- 2
m <- 3
v <- 1
X1 <- mvrnorm(10, c(-2,3), diag(v,p))
X2 <- mvrnorm(10, c(3,-2), diag(v,p))
X3 <- mvrnorm(10, c(-3,-3), diag(v,p))
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
X <- rbind(X1,X2,X3)

w1 <- Variable(p)
w2 <- Variable(p)
w3 <- Variable(p)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 3)
slack2 <- Variable(rows = nrow(X2), cols = 3)
slack3 <- Variable(rows = nrow(X3), cols = 3) 

xi1 <- Variable(nrow(X1))
xi2 <- Variable(nrow(X2))
xi3 <- Variable(nrow(X3)) 


C <- 1

objective <- Minimize(1/2*sum_squares(vstack(w1,w2,w3)) +  C*(sum_entries(xi1) + sum_entries(xi2) + sum_entries(xi3)))
constraints <- list(X1 %*% w1  + b1  >= 1 - slack1[,1],
                    X1 %*% w2  + b2 <= -1 + slack1[,2],
                    X1 %*% w3  + b3 <= -1 + slack1[,3],
                    X2 %*% w1  + b1 <= -1 + slack2[,1],
                    X2 %*% w2  + b2  >= 1 - slack2[,2],
                    X2 %*% w3  + b3 <= -1 + slack2[,3],
                    X3 %*% w1  + b1 <= -1 + slack3[,1],
                    X3 %*% w2  + b2 <= -1 + slack3[,2],
                    X3 %*% w3  + b3  >= 1 - slack3[,3] , 
                    xi1 >= slack1[,1],
                    xi1 >= sum_entries(slack1[,2:3], axis=1),
                    xi2 >= slack2[,2],
                    xi2 >= sum_entries(slack2[,c(1,3)], axis=1),
                    xi3 >= slack3[,3],
                    xi3 >= sum_entries(slack3[,1:2], axis=1),
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

system.time(MSVM8 <- Problem(objective, constraints))


CVXR_MSVM8 <- solve(MSVM8, solver = "MOSEK")

CVXR_MSVM8_w1 <- CVXR_MSVM8$getValue(w1)
CVXR_MSVM8_b1 <- CVXR_MSVM8$getValue(b1)
CVXR_MSVM8_w2 <- CVXR_MSVM8$getValue(w2)
CVXR_MSVM8_b2 <- CVXR_MSVM8$getValue(b2)
CVXR_MSVM8_w3 <- CVXR_MSVM8$getValue(w3)
CVXR_MSVM8_b3 <- CVXR_MSVM8$getValue(b3)

MSVM8_primary_beta <- rbind(c(CVXR_MSVM8_w1,CVXR_MSVM8_b1),
      c(CVXR_MSVM8_w2,CVXR_MSVM8_b2),
      c(CVXR_MSVM8_w3,CVXR_MSVM8_b3))


##############################3

class_idx <- sort(unique(y))
Y <- sapply(class_idx, function(id){as.numeric(y==id)})
n <- nrow(X)
p <- ncol(X)
m <- length(class_idx)

w <- Variable(rows = m, cols = p)
b <- Variable(rows = m)
slack <- Variable(rows = n, cols = m)

objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(slack))

constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(1- slack),
                    sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >=1- sum_entries(slack, axis=1),
                    slack >=0)

MSVM8 <- Problem(objective, constraints)
CVXR_MSVM8 <- solve(MSVM8, solver = "MOSEK")

MSVM8_primary_beta
cbind(CVXR_MSVM8$getValue(w),CVXR_MSVM8$getValue(b))
MSVM8_dual_opt(X,y,C)
MSVM8_pri_opt(X,y,C)




plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "MSVM8")

points(X2, col='green')

points(X3, col='blue')

abline(a = -MSVM8_primary_beta[1,3]/MSVM8_primary_beta[1,2], b = -MSVM8_primary_beta[1,1]/MSVM8_primary_beta[1,2])
abline(a = -MSVM8_primary_beta[2,3]/MSVM8_primary_beta[2,2], b = -MSVM8_primary_beta[2,1]/MSVM8_primary_beta[2,2])
abline(a = -MSVM8_primary_beta[3,3]/MSVM8_primary_beta[3,2], b = -MSVM8_primary_beta[3,1]/MSVM8_primary_beta[3,2])


plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "MSVM8")

points(X2, col='green')

points(X3, col='blue')

MSVM8_beta1_beta2 <- MSVM8_primary_beta[1,] - MSVM8_primary_beta[2,]
MSVM8_beta1_beta3 <- MSVM8_primary_beta[1,] - MSVM8_primary_beta[3,]
MSVM8_beta2_beta3 <- MSVM8_primary_beta[2,] - MSVM8_primary_beta[3,]



abline(a = -MSVM8_beta1_beta2[3]/MSVM8_beta1_beta2[2], b = -MSVM8_beta1_beta2[1]/MSVM8_beta1_beta2[2], lty =2, col = "red")
abline(a = -MSVM8_beta1_beta3[3]/MSVM8_beta1_beta3[2], b = -MSVM8_beta1_beta3[1]/MSVM8_beta1_beta3[2], lty =2, col = "red")
abline(a = -MSVM8_beta2_beta3[3]/MSVM8_beta2_beta3[2], b = -MSVM8_beta2_beta3[1]/MSVM8_beta2_beta3[2], lty =2)

# source("LLW.R")
# source("OVA.R")

# rm(list=ls())
#########################################################################

X <- rbind(X1,X2,X3)
n <- nrow(X)
m <- length(unique(y))
p <- ncol(X)

Y <- sapply(unique(y), function(id){as.numeric(y==id)})

Delta <- matrix(1, nrow  = n, ncol = m) - 2*Y
alpha <- Variable(rows = n, cols = m)

objective <- Maximize(-1/2*sum_entries((t(alpha*Delta)%*%X)^2) +  sum_entries(alpha))

# C = 1
constraints <- list(sum_entries(alpha*Delta, axis = 2) ==0,
                    alpha>=0, 
                    (alpha*(1-Y) + (alpha*Y)%*%matrix(1,nrow=3,ncol=3))<=C) 


MSVM8_dual <- Problem(objective, constraints)


system.time(CVXR_MSVM8_dual <- solve(MSVM8_dual, solver = "MOSEK"))

alpha_sol <- CVXR_MSVM8_dual$getValue(alpha)
w <- -t(alpha_sol*Delta)%*%X

lpred <- X%*%t(w)
alpha_sol <- round(alpha_sol,3)
balpha <- (alpha_sol*(1-Y) + (alpha_sol*Y)%*%matrix(1,nrow=3,ncol=3))
KKT1 <- (0 <= alpha_sol*(1-Y) & 0 <  balpha*(1-Y) &   balpha*(1-Y) < C)
KKT2 <- (0 <= alpha_sol*Y & alpha_sol*Y <= C & matrix(rep(apply(balpha*(1-Y), MARGIN = 1, FUN = function(x){max(x) <= C}),m), nrow=n,ncol=m))

b <- sapply(unique(y), function(id){ 
  bb <- NULL
  if(sum(KKT1[,id]) > 0){
    bb <- -1 - max(lpred[KKT1[,id],id])
  }
  if(is.null(bb) & sum(KKT2[,id]) > 0){
    bb <- 1 - min(lpred[KKT2[,id],id])
  }
  return(bb)
  }
)
b

MSVM8_dual_beta <- cbind(w,b)

MSVM8_primary_beta
MSVM8_dual_beta




  plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "MSVM8")

points(X2, col='green')

points(X3, col='blue')

abline(a = -MSVM8_dual_beta[1,3]/MSVM8_dual_beta[1,2], b = -MSVM8_dual_beta[1,1]/MSVM8_dual_beta[1,2])
abline(a = -MSVM8_dual_beta[2,3]/MSVM8_dual_beta[2,2], b = -MSVM8_dual_beta[2,1]/MSVM8_dual_beta[2,2])
abline(a = -MSVM8_dual_beta[3,3]/MSVM8_dual_beta[3,2], b = -MSVM8_dual_beta[3,1]/MSVM8_dual_beta[3,2])


plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "MSVM8")

points(X2, col='green')

points(X3, col='blue')

MSVM8_beta1_beta2 <- MSVM8_dual_beta[1,] - MSVM8_dual_beta[2,]
MSVM8_beta1_beta3 <- MSVM8_dual_beta[1,] - MSVM8_dual_beta[3,]
MSVM8_beta2_beta3 <- MSVM8_dual_beta[2,] - MSVM8_dual_beta[3,]



abline(a = -MSVM8_beta1_beta2[3]/MSVM8_beta1_beta2[2], b = -MSVM8_beta1_beta2[1]/MSVM8_beta1_beta2[2], lty =2, col = "red")
abline(a = -MSVM8_beta1_beta3[3]/MSVM8_beta1_beta3[2], b = -MSVM8_beta1_beta3[1]/MSVM8_beta1_beta3[2], lty =2, col = "red")
abline(a = -MSVM8_beta2_beta3[3]/MSVM8_beta2_beta3[2], b = -MSVM8_beta2_beta3[1]/MSVM8_beta2_beta3[2], lty =2)






# 
# res <- 0
# for(k in 1:m){
#   for( i in 1: n){
#     for( j in 1:n)
#     {
#       res <- res + Delta[i,k]*Delta[j,k]*alpha_sol[i,k]*alpha_sol[j,k]*t(X[i,])%*%X[j,]
#     }
#   }
# }
# res
# sum((t(alpha_sol*Delta)%*%X)^2)