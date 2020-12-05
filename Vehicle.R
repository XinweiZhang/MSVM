library(e1071)
library(CVXR)
library(MASS)
library(kernlab)
rm(list=ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("primary form functions.R")
source("SDCA functions.R")
source("kernel primary form functions.R")
set.seed(123)
Vehicle <- read.table("~/Desktop/Multiclass Classification/MSVM Code/Vehicle Data/tmp.dat", quote="\"", comment.char="")
X_dat <- as.matrix(Vehicle[,1:18])
y_dat <- as.numeric(Vehicle[,19])
s_idx <- sample(1:nrow(X_dat), round(nrow(X_dat)*0.7))

X <- X_dat[s_idx,]
y <- y_dat[s_idx]

X_test <- X_dat[-s_idx,]
y_test <- y_dat[-s_idx]


MSVM8_kernel_pri_opt(X, y, lambda = 1, intercept = T, kernel = rbfdot(sigma = .5))
