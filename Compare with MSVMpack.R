###### Test knab package ################3
library(CVXR)
library(MASS)
library(e1071)
library(kernlab)
library(R.matlab)
library(parallel)

rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/kernel primary form functions.R")

source("SDCA functions.R")

set.seed(123)
data_generate <- function(n, sep =  1,  v = 1.5^2, m){
  y <- sort(sample(seq(1,m), size = n-m, replace = T, prob =  1/rep(m,m)))
  n_list <- sapply(seq(1:m), function(t){sum(y==t)+1})
  X <- apply(cbind(mu_list, n_list), MARGIN = 1, FUN = function(mu){matrix(mvrnorm(mu[3], sep*mu[1:2], diag(v,nrow=2)), nrow = mu[3])})
  X <- do.call(rbind,X)
  y <- unlist(sapply(seq(1:m), function(t){rep(t,n_list[t])}))
  return(list(X = X,y = y))
}
m <- 3
mu_list <- matrix(rnorm(2*m, mean =0, sd = 10), m)
data <- data_generate(100, sep = 3, v = 10.5^2, m)
data_test <- data_generate(4000, sep = 3, v = 10.5^2, m)
X_train <- data$X
y_train <- data$y
X_test <- data_test$X
y_test <- data_test$y

X_train <- scale(X_train, center = T, scale = T)
X_test <- scale(X_test, center = attr(X_train,"scaled:center"), scale = attr(X_train,"scaled:scale"))
train_data <- as.data.frame(cbind(y = y_train, X = X_train))
test_data <- as.data.frame(cbind(y = y_test, X = X_test))





#####
ww_fit <- ksvm(factor(y) ~., data=train_data, kernel ="vanilladot", type = "kbb-svc", C = 1/400, scaled = F)
t(ww_fit@alpha)%*%X_train
# sum(c(0, ww_fit@b) - sum(ww_fit@b)/6)
# c(0, ww_fit@b) - sum(ww_fit@b)/6
ww_fit@b
WW_CVXR_fit <- WW_dp_pri_opt(X = X_train, y = y_train, C = 1, intercept = T)
WW_CVXR_fit$beta
sum(WW_CVXR_fit$beta[,3])

WW_CVXR_fit2 <- WW_pri_opt(X = X_train, y = y_train, C = 1, intercept = T)
WW_CVXR_fit2$beta

ww_fit <- ksvm(y ~., data=train_data, kernel ="vanilladot", type = "kbb-svc", C = 1/nrow(train_data), scaled = F)
ww_fit@b

cs_fit <- ksvm(y~., data=train_data,  type = "spoc-svc", C = nrow(train_data), scaled = F)


# for unscaled data


CS_CVXR_fit <- CS_pri_opt(X = X_train, y = y_train, C = 1, intercept = F)
CS_CVXR_fit$beta

mean(predict(CS_CVXR_fit, X_test) == y_test)

CS_SDCA_fit <- SDCA_kernel_CS(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(CS_SDCA_fit, X_test) == y_test)

CS_ksvm_fit <- ksvm(x = X_train, y = as.factor(y_train),kernel ="vanilladot",  type = "spoc-svc", C = 1)
mean(predict(CS_ksvm_fit, X_test) == y_test)

MDuchi_CVXR_fit <- MDuchi_pri_opt(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(MDuchi_CVXR_fit, X_test) == y_test)

MDuchi_SDCA_fit <- SDCA_kernel_MDuchi(X = X_train, y = y_train, C = 1, intercept = F, gap_cvg_tol = 1e-5)
mean(predict(MDuchi_SDCA_fit, X_test) == y_test)

Duchi_CVXR_fit <- Duchi_pri_opt(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(Duchi_CVXR_fit, X_test) == y_test)

Duchi_SDCA_fit <- SDCA_kernel_Duchi(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(Duchi_SDCA_fit, X_test) == y_test)

Y <- sapply(sort(unique(y_train)), function(id){as.numeric(y_train==id)})


CS_SDCA_fit <- SDCA_kernel_CS(X = X_train, y = y_train, C = 1, kernel = rbfdot(.5), intercept = F)
writeMat(con="./matlab code/data.mat", X = X_train, y = y_train, Y = Y, p = 2, m=m, C = 1, K = kernelMatrix(CS_SDCA_fit$kernel, X_train))


head(CS_SDCA_fit$v)

CS_CVXR_fit <- CS_kernel_pri_opt(X = X_train, y= y_train, C=1, intercept = F, kernel = rbfdot(.5))
t(rbind(t(CS_SDCA_fit$v)%*%kernelMatrix(CS_SDCA_fit$kernel, X_train)
,t(CS_CVXR_fit$v)%*%kernelMatrix(CS_SDCA_fit$kernel, X_train)))


t(CS_SDCA_fit$v)%*%X_train
