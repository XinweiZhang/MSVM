library(CVXR)
library(MASS)
library(R.matlab)
library(e1071)
library(kernlab)
library(ggpubr)
library(glmnet)
library(truncnorm)
library(ggplot2)
library(parallel)
rm(list=ls())

setwd("~/Desktop/Multiclass Classification/MSVM Code")
# source("primary form functions constraint.R")
source("primary form functions.R")


generate_data <- function(n, M, beta){
  
  t <- c(rep(1, round(n*(1 - beta - beta - nu))), rep(2, round(n*beta)), rep(3, round(n*beta)), rep(4, round(n*nu)))
  
  dat <- t(sapply(t, function(z){
    if(z==1){
      xz <- c(-1,0)
    }else if(z==2){
      xz <- c(-1+M,0)
    }else if(z==3){
      xz <- c(-1+2*M,0)
    }else if(z==4){
      xz <- c(1,0)
    }
    
    if(z==1){
      y = 1
    }else if(z==2){
      y = 2
    }else if(z==3){
      y = 3
    }else if(z==4){
      y = 1
    }
    return(c(xz,y))
  }))
  
  X <- as.matrix(dat[,c(1,2)])
  y <- dat[,3]
  
  list(X=X,y=y)
}

nu = 0.15
M_list <- seq(1e-3,1,length.out = 500)
beta_list <- seq(nu,(1-nu)/3, length.out = 100)

error_mat <- matrix(nrow = 10, ncol =500)

cl <- makeCluster(15)

clusterExport(cl=cl, varlist=c("WW_pri_opt", "Duchi_pri_opt","MDuchi_pri_opt", "CS_pri_opt","nu",
                               "LLW_pri_opt", "OVA_pri_opt", "MSVM8_pri_opt","MSVM7_pri_opt","New1_pri_opt","New3_pri_opt","ginv",
                               "solve", "Minimize", "sum_squares", "sum_entries", "vstack", "reshape_expr",
                               "Variable","Problem","max_entries","min_entries","predict.msvm","generate_data","mvrnorm","Maximize","vec"), envir=environment())

for(i in 99:length(M_list)){
  print(i)
  err <- parSapply(cl, beta_list, function(beta, M){
    
    tryCatch({
      
      n_train <- 1000
      n_test <- 10000
      train_data <- generate_data(n_train, M, beta)
      test_data <- generate_data(n_test, M, beta)
      intercept = T
      C = 1
      lambda = 0
      WW_fit <- WW_pri_opt(train_data$X,train_data$y, C=C, lambda=lambda,  intercept = intercept)
      predicted_label <- predict(WW_fit, test_data$X)
      err1 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      CS_fit <- CS_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda, intercept = intercept)
      predicted_label <- predict(CS_fit, test_data$X)
      err2 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      
      Duchi_fit <- Duchi_pri_opt(train_data$X,train_data$y,  C= C, lambda=lambda, intercept = intercept)
      predicted_label <- predict(Duchi_fit, test_data$X)
      err3 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      
      MDuchi_fit <- MDuchi_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda,  intercept = intercept)
      predicted_label <- predict(MDuchi_fit, test_data$X)
      err4 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      New1_fit <- New1_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda,  intercept = intercept, base_class = 1)
      predicted_label <- predict(New1_fit, test_data$X, "dagger")
      err_New1_base1 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      New1_fit <- New1_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda,  intercept = intercept, base_class = 2)
      predicted_label <- predict(New1_fit, test_data$X, "dagger")
      err_New1_base2 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      New1_fit <- New1_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda,  intercept = intercept, base_class = 3)
      predicted_label <- predict(New1_fit, test_data$X, "dagger")
      err_New1_base3 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      
      New3_fit <- New3_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda,  intercept = intercept, base_class = 1)
      predicted_label <- predict(New3_fit, test_data$X, "dagger")
      err_New3_base1 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      New3_fit <- New3_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda,  intercept = intercept, base_class = 2)
      predicted_label <- predict(New3_fit, test_data$X, "dagger")
      err_New3_base2 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      New3_fit <- New3_pri_opt(train_data$X,train_data$y,  C=C, lambda=lambda,  intercept = intercept, base_class = 3)
      predicted_label <- predict(New3_fit, test_data$X, "dagger")
      err_New3_base3 <- round((1-mean(predicted_label == test_data$y))*100,3)
      
      
      return(c(err1,err2,err3,err4, err_New1_base1, err_New1_base2, err_New1_base3, err_New3_base1, err_New3_base2, err_New3_base3))
    }, error = function(e){
      return(rep(NA,10))
    })
    
  }, M = M_list[i])
  
  error_mat[,i] <- apply(err, MARGIN = 1, max, na.rm = TRUE)
}

stopCluster(cl)

save(error_mat, file = "Error_rate_trend (nu=0.15).Rdata")
source("./bound draft/bound trend plot.R")

error_mat[4,]

plot(x=4/M_list, y= error_mat[1,], type = "l", xlab = "B", ylab = "error rate", ylim = c(0,100), xlim = c(4,200), cex=3)
lines(x=4/M_list, y= error_mat[2,], col = "red")
lines(x=4/M_list, y= error_mat[3,], col = "blue")
lines(x=4/M_list, y= error_mat[4,], col = "green")
legend(100, 40, legend=c("WW", "CS", "DUCHI", "MDUCHI"),
       col=c("black", "red", "blue","green"), lty=1, cex=1)

plot(x=M_list, y= error_mat[1,], type = "l", xlab = "M", ylab = "error rate", ylim = c(0,100), xlim = c(0,1), cex=3)
lines(x=M_list, y= error_mat[2,], col = "red")
lines(x=M_list, y= error_mat[3,], col = "blue")
lines(x=M_list, y= error_mat[4,], col = "green")
legend(100, 40, legend=c("WW", "CS", "DUCHI", "MDUCHI"),
       col=c("black", "red", "blue","green"), lty=1, cex=1)


plot(x=4/M_list, y= error_mat[4,], type = "l", xlab = "B", ylab = "error rate", ylim = c(0,100))
lines(x=4/M_list, y= error_mat[5,], col = "red")
lines(x=4/M_list, y= error_mat[6,], col = "blue")
lines(x=4/M_list, y= error_mat[7,], col = "green")
legend(100, 40, legend=c("MDuchi", "New1_base1", "New1_base2", "New1_base3"),
       col=c("black", "red", "blue","green"), lty=1, cex=1)



plot(x=4/M_list, y= error_mat[3,], type = "l", xlab = "B", ylab = "error rate", ylim = c(0,100))
lines(x=4/M_list, y= error_mat[8,], col = "red")
lines(x=4/M_list, y= error_mat[9,], col = "blue")
lines(x=4/M_list, y= error_mat[10,], col = "green")
legend(100, 40, legend=c("MDuchi", "New3_base1", "New3_base2", "New3_base3"),
       col=c("black", "red", "blue","green"), lty=1, cex=1)

