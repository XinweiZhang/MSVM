library(CVXR)
library(MASS)
library(glmnet)
set.seed(1223)
rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("primary form functions.R")
source("Dual form functions.R")



#########Toy Example
X = rbind(c(-3,5),c(-4,-1),c(2,-2),c(5,4))
y = c(1,2,3,4)
C = 20
plot_nearest_neighbour_decision_boundary(X, y, dis = "L2", title = "Nearest Neighbour in L2")
plot_nearest_neighbour_decision_boundary(X, y, dis = "L1", title = "Nearest Neighbour in L1")
plot_decision_boundary(X,y,beta = WW_pri_opt(X,y,C), title = "WW")
plot_decision_boundary(X,y,beta = CS_pri_opt(X,y,C), title = "CS")
plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C), title = "Duchi")
plot_decision_boundary(X,y,beta = MDuchi_pri_opt(X,y,C), title = "Mduchi")
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1, 0", dagger_rule = F)
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1, dagger", dagger_rule = T)
plot_decision_boundary(X,y,beta = New3_pri_opt(X,y,C), title = "New3, 0", dagger_rule = F)
plot_decision_boundary(X,y,beta = OVA_pri_opt(X,y,C), title = "OVA")
# debugonce(plot_decision_boundary)
plot_decision_boundary(X,y,beta = LLW_pri_opt(X,y,C), title = "LLW", xlim=c(-30,30),  ylim=c(-30,30))


plot_decision_boundary(X,y,beta = MSVM7_pri_opt(X,y,C), title = "Modifed LLW")
plot_decision_boundary(X,y,beta = MSVM8_pri_opt(X,y,C), title = "MSVM8")




#########Toy Example

X = rbind(c(-2,1),c(-2,-1),c(2,1))
y = c(1,2,3)
C = 1
# Duchi_pri_opt(X,y,C)
w =  t(rbind(c(-0.249850,     0.250150,    -0.000300), 
c(0.499886,     0.000227,    -0.500113)))
b <- c(-0.332920, 0.166740, 0.166180)
w
b
beta <- cbind(w,b)
beta <- WW_pri_opt(X,y,C)
w <- beta[,1:2]
b <- beta[,3]

# fit = glmnet(X, y, family = "multinomial", alpha = 1, lambda = .000001)
# call(rbind,coef(fit))
# coef(fit)

# beta <- cbind(w,b)
plot_decision_boundary(X,y,beta = WW_pri_opt(X,y,C), title = "WW")
# dev.off()
plot_decision_boundary(X,y,beta = beta, title = "WW")


beta = WW_pri_opt(X,y,C)

(c(-2,1)%*%(w[1,1:2]-w[2,1:2]) + b[1] - b[2])/sqrt(sum((w[1,1:2]-w[2,1:2])^2))
(c(-2,1)%*%(w[1,1:2]-w[3,1:2]) + b[1] - b[3])/sqrt(sum((w[1,1:2]-w[2,1:2])^2))
(c(2,1)%*%(w[3,1:2]-w[1,1:2]) + b[3] - b[1])/sqrt(sum((w[3,1:2]-w[1,1:2])^2))
(c(2,1)%*%(w[3,1:2]-w[2,1:2]) + b[3] - b[2])/sqrt(sum((w[3,1:2]-w[2,1:2])^2))

(c(-2,-1)%*%(w[2,1:2]-w[1,1:2]) + b[2] - b[1])/sqrt(sum((w[2,1:2]-w[1,1:2])^2))
(c(-2,-1)%*%(w[2,1:2]-w[3,1:2]) + b[2] - b[3])/sqrt(sum((w[2,1:2]-w[3,1:2])^2))



plot_decision_boundary(X,y,beta = Duchi_dual_opt(X,y,C), title = "WW")

plot_decision_boundary(X,y,beta = CS_pri_opt(X,y,C), title = "CS")
plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C), title = "Duchi")
plot_decision_boundary(X,y,beta = MDuchi_pri_opt(X,y,C), title = "Mduchi")
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1", dagger_rule = F)


X = rbind(c(-2,1),c(-2,-1),c(2,1),c(2,-2))
y = c(1,2,3,4)
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1", dagger_rule = T)
X = rbind(c(-2,1),c(-2,-1),c(2,-2),c(2,1))
y = c(1,2,3,4)
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1", dagger_rule = T)
plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C), title = "Duchi")

plot_decision_boundary(X,y,beta = New3_pri_opt(X,y,C), title = "New3")

plot_decision_boundary(X,y,beta = OVA_pri_opt(X,y,C), title = "OVA")
plot_decision_boundary(X,y,beta = LLW_pri_opt(X,y,C), title = "LLW")
plot_decision_boundary(X,y,beta = MSVM7_pri_opt(X,y,C), title = "Modifed LLW")
plot_decision_boundary(X,y,beta = MSVM8_pri_opt(X,y,C=100), title = "MSVM8")

###########3 class Gaussian
set.seed(123)
p <- 2
m <- 3
X1 <- mvrnorm(5,c(0,5),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(5,c(5,0),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(5,c(0,0),matrix(c(5,0,0,5),nrow=2))
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
C <- 1
plot_decision_boundary(X,y,beta = WW_pri_opt(X,y,C), title = "WW")
plot_decision_boundary(X,y,beta = CS_pri_opt(X,y,C), title = "CS")
plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C), title = "Duchi")
plot_decision_boundary(X,y,beta = MDuchi_pri_opt(X,y,C), title = "Mduchi")
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1")
plot_decision_boundary(X,y,beta = New3_pri_opt(X,y,C), title = "New3")

plot_decision_boundary(X,y,beta = OVA_pri_opt(X,y,C), title = "OVA")
plot_decision_boundary(X,y,beta = LLW_pri_opt(X,y,C), title = "LLW")
plot_decision_boundary(X,y,beta = MSVM7_pri_opt(X,y,C), title = "Modified LLW")
plot_decision_boundary(X,y,beta = MSVM8_pri_opt(X,y,C), title = "MSVM8")


###########4 class 
set.seed(1223)
# set.seed(1234)
p <- 2
m <- 4
Sigma <- diag(4,nrow=2)
X1 <- mvrnorm(10,c(0,2),Sigma)
X2 <- mvrnorm(10,c(2,0),Sigma)
X3 <- mvrnorm(10,c(-2,0),Sigma)
X4 <- mvrnorm(10,c(0,-2),Sigma)
C <- 1
X <- rbind(X1,X2,X3,X4)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X3)))

plot_decision_boundary(X,y,beta = WW_pri_opt(X,y,C), title = "WW")
plot_decision_boundary(X,y,beta = CS_pri_opt(X,y,C), title = "CS")
plot_decision_boundary(X,y,beta = Duchi_pri_opt(X,y,C), title = "Duchi")
plot_decision_boundary(X,y,beta = MDuchi_pri_opt(X,y,C), title = "Mduchi")
plot_decision_boundary(X,y,beta = New1_pri_opt(X,y,C), title = "New1")
plot_decision_boundary(X,y,beta = New3_pri_opt(X,y,C), title = "New3")

plot_decision_boundary(X,y,beta = OVA_pri_opt(X,y,C), title = "OVA")
plot_decision_boundary(X,y,beta = LLW_pri_opt(X,y,C), title = "LLW")
plot_decision_boundary(X,y,beta = MSVM7_pri_opt(X,y,C), title = "MSVM7", xlim = c(-5,5), ylim = c(-8,8))
plot_decision_boundary(X,y,beta = MSVM8_pri_opt(X,y,C), title = "MSVM8")


###########3 class 
set.seed(1223)
# set.seed(1234)
p <- 2
Sigma <- diag(10,nrow=2)
X1 <- mvrnorm(1000,c(-2,1),Sigma)
X2 <- mvrnorm(1000,c(-2,-1),Sigma)
X3 <- mvrnorm(1000,c(2,1),Sigma)
C <- 1
X <- rbind(X1,X2,X3)
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))

beta = WW_pri_opt(X,y,C)
beta = beta/sum(apply(beta, MARGIN =1, function(x){sqrt(sum(x[1:2]^2))}))
plot_decision_boundary(X,y,beta = Duchi_dual_opt(X,y,C), title = "WW")

w <- beta[,1:2]
b <- beta[,3]
# (c(-2,1)%*%(w[1,1:2]-w[2,1:2]) + b[1] - b[2])/sqrt(sum((w[1,1:2]-w[2,1:2])^2))
# (c(-2,1)%*%(w[1,1:2]-w[3,1:2]) + b[1] - b[3])/sqrt(sum((w[1,1:2]-w[2,1:2])^2))
# (c(2,1)%*%(w[3,1:2]-w[1,1:2]) + b[3] - b[1])/sqrt(sum((w[3,1:2]-w[1,1:2])^2))
# (c(2,1)%*%(w[3,1:2]-w[2,1:2]) + b[3] - b[2])/sqrt(sum((w[3,1:2]-w[2,1:2])^2))
# 
# (c(-2,-1)%*%(w[2,1:2]-w[1,1:2]) + b[2] - b[1])/sqrt(sum((w[2,1:2]-w[1,1:2])^2))
# (c(-2,-1)%*%(w[2,1:2]-w[3,1:2]) + b[2] - b[3])/sqrt(sum((w[2,1:2]-w[3,1:2])^2))

(c(-2,1)%*%(w[1,1:2]-w[2,1:2]) + b[1] - b[2])
(c(-2,1)%*%(w[1,1:2]-w[3,1:2]) + b[1] - b[3])
(c(2,1)%*%(w[3,1:2]-w[1,1:2]) + b[3] - b[1])
(c(2,1)%*%(w[3,1:2]-w[2,1:2]) + b[3] - b[2])

(c(-2,-1)%*%(w[2,1:2]-w[1,1:2]) + b[2] - b[1])
(c(-2,-1)%*%(w[2,1:2]-w[3,1:2]) + b[2] - b[3])


a <- c(0,1,0)
(c(-2,1,1)%*%a)/sqrt(sum((a[1:2])^2))
a <- -c(0,1,0)
(c(-2,-1,1)%*%a)/sqrt(sum((a[1:2])^2))
a <- c(2,0,0)
(c(2,1,1)%*%a)/sqrt(sum((a[1:2])^2))
a <- -c(2,0,0)
(c(-2,1,1)%*%a)/sqrt(sum((a[1:2])^2))
a <- c(2,1,0)
(c(2,1,1)%*%a)/sqrt(sum((a[1:2])^2))
a <- -c(2,1,0)
(c(-2,-1,1)%*%a)/sqrt(sum((a[1:2])^2))


mat <- solve(rbind(c(1,-1,0), c(-1,0,1), c(1,1,1)))%*%rbind(c(0,1,0),c(2,0,0),c(0,0,0))
# sum(apply(mat, MARGIN =1, function(x){sqrt(sum(x[1:2]^2))}))
# sum(apply(beta, MARGIN =1, function(x){sqrt(sum(x[1:2]^2))}))
mat <- mat/sum(apply(mat, MARGIN =1, function(x){sqrt(sum(x[1:2]^2))}))
plot_decision_boundary(X,y,beta = mat, title = "perpendicular")



w <- mat[,1:2]
b <- mat[,3]

(c(-2,1)%*%(w[1,1:2]-w[2,1:2]) + b[1] - b[2])
(c(-2,1)%*%(w[1,1:2]-w[3,1:2]) + b[1] - b[3])
(c(2,1)%*%(w[3,1:2]-w[1,1:2]) + b[3] - b[1])
(c(2,1)%*%(w[3,1:2]-w[2,1:2]) + b[3] - b[2])

(c(-2,-1)%*%(w[2,1:2]-w[1,1:2]) + b[2] - b[1])
(c(-2,-1)%*%(w[2,1:2]-w[3,1:2]) + b[2] - b[3])

