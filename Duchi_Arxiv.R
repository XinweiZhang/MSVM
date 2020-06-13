library(CVXR)
library(MASS)
set.seed(123)
rm(list = ls())
d <- 2
m <- 3
X1 <- mvrnorm(2,c(0,5),matrix(c(5,0,0,5),nrow=2))
X2 <- mvrnorm(2,c(5,0),matrix(c(5,0,0,5),nrow=2))
X3 <- mvrnorm(2,c(0,0),matrix(c(5,0,0,5),nrow=2))
y <- c(rep(1,5),rep(2,5),rep(3,5))

profvis({
w1 <- Variable(d)
w2 <- Variable(d)
w3 <- Variable(d)
b1 <- Variable(1)
b2 <- Variable(1)
b3 <- Variable(1)


slack1 <- Variable(rows = nrow(X1), cols = 3)
slack2 <- Variable(rows = nrow(X2), cols = 3)
slack3 <- Variable(rows = nrow(X3), cols = 3)


g <- function(alpha){
  max(max(alpha)-1, sum_largest(alpha, 2)/2 - 1/2 , sum(alpha)/3-1/3)
}

gslack1 <- lapply(seq(1,nrow(X1)), FUN = function(i){g(slack1[i,])})
gslack2 <- lapply(seq(1,nrow(X2)), FUN = function(i){g(slack2[i,])})
gslack3 <- lapply(seq(1,nrow(X3)), FUN = function(i){g(slack3[i,])})


C <- 1
objective <- Minimize(sum_squares(vstack(w1,w2,w3)) +  C*(Reduce(f = sum, x = gslack1) + Reduce(f = sum, x = gslack2) + Reduce(f = sum, x = gslack3)))
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
                    slack1 >=0,
                    slack2 >=0,
                    slack3 >=0)

Duchi <- Problem(objective, constraints)
CVXR_Duchi <- solve(Duchi)
})


CVXR_Duchi_w1 <- CVXR_Duchi$getValue(w1)
CVXR_Duchi_b1 <- CVXR_Duchi$getValue(b1)
CVXR_Duchi_w2 <- CVXR_Duchi$getValue(w2)
CVXR_Duchi_b2 <- CVXR_Duchi$getValue(b2)
CVXR_Duchi_w3 <- CVXR_Duchi$getValue(w3)
CVXR_Duchi_b3 <- CVXR_Duchi$getValue(b3)

beta1 <- c(CVXR_Duchi_w1,CVXR_Duchi_b1)
beta2 <- c(CVXR_Duchi_w2,CVXR_Duchi_b2)
beta3 <- c(CVXR_Duchi_w3,CVXR_Duchi_b3)

plot(X1, col='red', xlim = c(-10,10), ylim=c(-10,10))

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


