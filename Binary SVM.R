rm(list = ls())
set.seed(123)
p <- 2
m <- 3
v <- 3
X1 <- mvrnorm(10, c(-2,3), diag(v,p))
X2 <- mvrnorm(10, c(3,-2), diag(v,p))

y <- c(rep(1,nrow(X1)),rep(-1,nrow(X2)))

w <- Variable(p)
b <- Variable(1)

slack <- Variable(rows=nrow(rbind(X1,X2)), cols=1)
                   
slack1 <- Variable(nrow(X1))
slack2 <- Variable(nrow(X2))

C <- 1
X <-rbind(X1,X2)
objective <- Minimize(1/2*sum_squares(w) +  C*sum_entries(slack))
constraints <- list( 1- y*(X %*% w  + b)  <=  slack,
                    slack >=0)

SVM <- Problem(objective, constraints)
SVM_sol <- solve(SVM)

# SVM_sol$getValue(w)
# SVM_sol$getValue(b)
# SVM_sol$getValue(slack)
SVM_sol$getDualValue(constraints(SVM)[[1]])

beta <- c(SVM_sol$getValue(w),SVM_sol$getValue(b))

par(mfrow=c(1,1))
plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "OVA")
points(X2, col='green')
abline(a = -beta[3]/beta[2], b = -beta[1]/beta[2])

# t(SVM_sol$getDualValue(constraints(SVM)[[1]]))
# round(t(alpha_sol[,1]),2)
# 



X <- rbind(X1,X2)
n <- nrow(X)
m <- length(unique(y))
p <- ncol(X)

Y <- sapply(unique(y), function(id){as.numeric(y==id)})
Delta <- matrix(1, nrow  = n, ncol = m) - 2*Y
# Delta <-  2*Y - matrix(1, nrow  = n, ncol = m)
alpha <- Variable(rows = n, cols = m)
objective <- Maximize((-sum_entries((t(alpha*Delta)%*%X)^2)/2 +  sum_entries(alpha)))

C = 1
constraints <- list( 
                    # sum(alpha[1:10,1])-sum(alpha[11:20,1]) ==0,
                    #  sum(alpha[1:10,2])-sum(alpha[11:20,2]) ==0,
                      t(alpha*Delta)%*%matrix(1, nrow=n, ncol = 1) ==0,
                      alpha<=C,
                      # alpha*(1-Y) + sum_entries(alpha*Y,axis=1)%*%matrix(1,ncol=m,nrow=1)<=C,
                      alpha>=0) 


OVA_dual <- Problem(objective, constraints)
CVXR_OVA_dual <- solve(OVA_dual)
alpha_sol <- CVXR_OVA_dual$getValue(alpha)

# CVXR_OVA_dual$getDualValue(constraints)

# round(alpha_sol ,1)
w2 <- colSums(alpha_sol*cbind(y,y)*X)
lpred <- X%*%w2

w2
beta

round(alpha_sol[,1],2)

b2 <- -max(lpred[((C - Y*alpha_sol) > .01 & (Y*alpha_sol > 0.01) == T)[,1]])
b2 <- b2 - min(lpred[((C - (1-Y)*alpha_sol) > .01 & ((1-Y)*alpha_sol > 0.01) == T)[,1]])

b2 <- b2/2
beta[3]
