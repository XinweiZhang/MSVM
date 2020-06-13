library(CVXfromR)
set.seed(136)
library(MASS)
rm(list = ls())
p <- 2
m <- 4
X1 <- mvrnorm(10,c(0,5),matrix(c(10,0,0,10),nrow=2))
X2 <- mvrnorm(10,c(5,0),matrix(c(10,0,0,10),nrow=2))
X3 <- mvrnorm(10,c(-5,0),matrix(c(10,0,0,10),nrow=2))
X4 <- mvrnorm(10,c(0,-4),matrix(c(10,0,0,10),nrow=2))
setup.dir <- "/usr/local/MATLAB/R2020a/toolbox/cvx"
C=1
n1 = nrow(X1);
n2 = nrow(X2);
n3 = nrow(X3);
n4 = nrow(X4);

cvxcode <- paste("variables w1(p) w2(p) w3(p) w4(p) b1 b2 b3 b4",
      "variables slack1(n1,m) slack2(n2,m) slack3(n3,m) slack4(n4,m)",
      "variables xi1(n1) xi2(n2) xi3(n3) xi4(n4);minimize( sum_square([w1;w2;w3;w4])/2 + C* sum([xi1;xi2;xi3;xi4]))",
      "subject to","w1 + w2 + w3 + w4 == 0","b1 + b2 + b3 + b4 == 0",
        "slack1(:,1) == 1",
        "X1*(w1 - w2) + b1 - b2 >= 1-slack1(:,2)",
      "X1*(w1 - w3) + b1 - b3 >= 1-slack1(:,3)",
        "X1*(w1 - w4) + b1 - b4 >= 1-slack1(:,4)",
      "slack2(:,2) == 1",
       " X2*(w2 - w1) + b2 - b1 >= 1-slack2(:,1)",
      "X2*(w2 - w3) + b2 - b3 >= 1-slack2(:,3)",
        "X2*(w2 - w4) + b2 - b4 >= 1-slack2(:,4)",
      "slack3(:,3) == 1",
        "X3*(w3 - w1) + b3 - b1 >= 1-slack3(:,1)",
      "X3*(w3 - w2) + b3 - b2 >= 1-slack3(:,2)",
        "X3*(w3 - w4) + b3 - b4 >= 1-slack3(:,4)",
      "slack4(:,4) == 1",
        "X4*(w4 - w1) + b4 - b1 >= 1-slack4(:,1)",
      "X4*(w4 - w2) + b4 - b2 >= 1-slack4(:,2)",
        "X4*(w4 - w3) + b4 - b3 >= 1-slack4(:,3)",
      "slack1 >= 0",
        "slack2 >= 0",
      "slack3 >= 0",
        "slack4 >= 0",
      "xi1 >= norms_largest(slack1,1,2)-1", 
       " xi2 >= norms_largest(slack2,1,2)-1",
      "xi3 >= norms_largest(slack3,1,2)-1",
       " xi4 >= norms_largest(slack4,1,2)-1",
      "xi1 >= norms_largest(slack1,2,2)/2-1/2",
       " xi2 >= norms_largest(slack2,2,2)/2-1/2",
      "xi3 >= norms_largest(slack3,2,2)/2-1/2",
       " xi4 >= norms_largest(slack4,2,2)/2-1/2",
      "xi1 >= norms_largest(slack1,3,2)/3-1/3",
       " xi2 >= norms_largest(slack2,3,2)/3-1/3",
      "xi3 >= norms_largest(slack3,3,2)/3-1/3",
       " xi4 >= norms_largest(slack4,3,2)/3-1/3",
      "xi1 >= norms_largest(slack1,4,2)/4-1/4",
       " xi2 >= norms_largest(slack2,4,2)/4-1/4",
      "xi3 >= norms_largest(slack3,4,2)/4-1/4",
       " xi4 >= norms_largest(slack4,4,2)/4-1/4", sep = ";")




# lasso <- CallCVX(cvxcode, const.vars=list(X1=X1, X2=X2, X3=X3, X4=X4, p=p, m=m, C = C, n1=n1,n2=n2,n3=n3,n4=n4),
#                  opt.var.names=c("w1","w2","w3","w4","b1","b2","b3","b4","xi1","xi2","xi3","xi4","slack1","slack2","slack3","slack4"), setup.dir=setup.dir)
# 
lasso <- CallCVX(cvxcode, const.vars=list(X1=X1, X2=X2, X3=X3, X4=X4, p=p, m=m, C = C, n1=n1,n2=n2,n3=n3,n4=n4),
                 opt.var.names=c("w1","w2","w3","w4","b1","b2","b3","b4"), setup.dir=setup.dir, cvx.modifiers="quiet")

