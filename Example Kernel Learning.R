library(CVXR)
library(MASS)
library(ggplot2)
library(kernlab)

rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code/")
source("~/Desktop/Multiclass Classification/MSVM Code/kernel primary form functions.R")
# source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
set.seed(123)
C=1

kernel = vanilladot()

# kernel = rbfdot(sigma = 1)
# kernel = laplacedot(sigma = 1)
# kernel = polydot(degree = 2, scale = 1, offset = 1)

# dat_train <- three_class_data_generate(50, sep = 1)
# dat_test <- three_class_data_generate(1000, sep = 1)

# dat_train <- four_class_data_generate(100, sep = 1.5)
# dat_test <- four_class_data_generate(1000, sep = 1.5)

dat_train <- five_class_data_generate(50, sep = 1.5)
dat_test <- five_class_data_generate(1000, sep = 1.5)

# debugonce(WW_kernel_pri_opt)
system.time(WW_fit <- WW_kernel_pri_opt(dat_train$X, dat_train$y, C = C, kernel))
system.time(CS_fit <- CS_kernel_pri_opt(dat_train$X, dat_train$y, C = C, kernel))
system.time(Duchi_fit <- Duchi_kernel_pri_opt(dat_train$X, dat_train$y, C = C,  kernel))
system.time(MDuchi_fit <- MDuchi_kernel_pri_opt(dat_train$X, dat_train$y, C = C,  kernel))
system.time(New1_fit <- New1_kernel_pri_opt(dat_train$X, dat_train$y, C = C,  kernel))
system.time(New3_fit <- New3_kernel_pri_opt(dat_train$X, dat_train$y, C = C,  kernel))
system.time(MSVM8_fit <- MSVM8_kernel_pri_opt(dat_train$X, dat_train$y, C = C, kernel))
system.time(OVA_fit <- OVA_kernel_pri_opt(dat_train$X, dat_train$y, C = C, kernel))
system.time(MSVM7_fit <- MSVM7_kernel_pri_opt(dat_train$X, dat_train$y, C = C,  kernel))
system.time(LLW_fit <- LLW_kernel_pri_opt(dat_train$X, dat_train$y, C = C,  kernel))
#########

mean(predict(WW_fit, dat_test$X, rule = "simple_max") == dat_test$y)
mean(predict(CS_fit, dat_test$X, rule = "simple_max") == dat_test$y)
mean(predict(Duchi_fit, dat_test$X,  rule = "simple_max") == dat_test$y)
mean(predict(MDuchi_fit, dat_test$X,  rule = "simple_max") == dat_test$y)


mean(predict(New1_fit, dat_test$X,  rule = "dagger_New1") == dat_test$y)
mean(predict(New1_fit, dat_test$X,  rule = "0_max") == dat_test$y)
mean(predict(New3_fit, dat_test$X,  rule = "dagger_New1") == dat_test$y)
mean(predict(New3_fit, dat_test$X,  rule = "0_max") == dat_test$y)
mean(predict(MSVM8_fit, dat_test$X, rule = "simple_max") == dat_test$y)
mean(predict(OVA_fit, dat_test$X, rule = "simple_max") == dat_test$y)
mean(predict(MSVM7_fit, dat_test$X,  rule = "dagger_MSVM7") == dat_test$y)
mean(predict(MSVM7_fit, dat_test$X,  rule = "0_max") == dat_test$y)
mean(predict(LLW_fit, dat_test$X,  rule = "simple_max") == dat_test$y)
# 
plot_kernel_decision_boundary(WW_fit, title = "WW")
plot_kernel_decision_boundary(CS_fit, title = "CS" )
plot_kernel_decision_boundary(Duchi_fit, title = "Duchi" )
plot_kernel_decision_boundary(MDuchi_fit, title = "MDuchi" )
plot_kernel_decision_boundary(New1_fit, title = "New1", rule = "dagger_New1")
plot_kernel_decision_boundary(New1_fit, title = "New1", rule = "dagger_New1")
plot_kernel_decision_boundary(New3_fit, title = "New3", rule = "0_max")
plot_kernel_decision_boundary(New3_fit, title = "New3", rule = "0_max")
plot_kernel_decision_boundary(MSVM7_fit, title = "MSVM7", rule = "dagger_MSVM7")
plot_kernel_decision_boundary(MSVM8_fit, title = "MSVM8" )
plot_kernel_decision_boundary(OVA_fit, title = "OVA" )
plot_kernel_decision_boundary(LLW_fit, title = "LLW" )
    
