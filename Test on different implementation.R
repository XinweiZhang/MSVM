###### Test knab package ################3
library(CVXR)
library(MASS)
library(e1071)
library(kernlab)
library(parallel)

rm(list = ls())
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("SDCA functions.R")

# set.seed(123)
data_generate <- function(n, sep =  1,  v = 1.5^2, m, mu_list){
  y <- sort(sample(seq(1,m), size = n-m, replace = T, prob =  1/rep(m,m)))
  n_list <- sapply(seq(1:m), function(t){sum(y==t)+1})
  X <- apply(cbind(mu_list, n_list), MARGIN = 1, FUN = function(mu){matrix(mvrnorm(mu[3], sep*mu[1:2], diag(v,nrow=2)), nrow = mu[3])})
  X <- do.call(rbind,X)
  y <- unlist(sapply(seq(1:m), function(t){rep(t,n_list[t])}))
  return(list(X = X,y = y))
}
m <- 4
mu_list <- matrix(rnorm(2*m, mean =0, sd = 10), m)
data <- data_generate(200, sep = 3, v = 10.5^2, m, mu_list)
data_test <- data_generate(2000, sep = 3, v = 10.5^2, m, mu_list)
X_train <- data$X
y_train <- data$y
X_test <- data_test$X
y_test <- data_test$y

X_train <- scale(X_train, center = T, scale = T)
X_test <- scale(X_test, center = attr(X_train,"scaled:center"), scale = attr(X_train,"scaled:scale"))
train_data <- as.data.frame(cbind(y = y_train, X = X_train))
test_data <- as.data.frame(cbind(y = y_test, X = X_test))



#####
# ww_fit <- ksvm(y ~., data=train_data, kernel ="vanilladot", type = "kbb-svc", C = 1/400, scaled = F)
# ww_fit@b
# 
# WW_CVXR_fit <- WW_pri_opt(X = X_train, y = y_train, C = 1, intercept = T)
# WW_CVXR_fit$beta
# 
# ww_fit <- ksvm(y ~., data=train_data, kernel ="vanilladot", type = "kbb-svc", C = 1/nrow(train_data), scaled = F)
# ww_fit@b
# 
# cs_fit <- ksvm(y~., data=train_data,  type = "spoc-svc", C = nrow(train_data), scaled = F)
# 

# for unscaled data

WW_CVXR_fit <- WW_pri_opt(X = X_train, y = y_train, C = 1, intercept = T)
mean(predict(WW_CVXR_fit, X_test) == y_test)

mean(MSVMPack_interface("WW", X_train, y_train, X_test, y_test, C = 1, kernel = 1, sigma = 0) == y_test)

WW_ksvm_fit <- ksvm(x = X_train, y = as.factor(y_train), kernel = vanilladot(),  type = "kbb-svc", C = 1)
mean(predict(WW_ksvm_fit, X_test) == y_test)


CS_CVXR_fit <- CS_pri_opt(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(CS_CVXR_fit, X_test) == y_test)

mean(MSVMPack_interface("CS", X_train, y_train, X_test, y_test, C = 1, kernel = 1, sigma = 0) == y_test)

CS_SDCA_fit <- SDCA_kernel_CS(X = X_train, y = y_train, C = 1, intercept = T)
mean(predict(CS_SDCA_fit, X_test) == y_test)

CS_ksvm_fit <- ksvm(x = X_train, y = as.factor(y_train), kernel = vanilladot(),  type = "spoc-svc", C = 1)
mean(predict(CS_ksvm_fit, X_test) == y_test)

MDuchi_CVXR_fit <- MDuchi_pri_opt(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(MDuchi_CVXR_fit, X_test) == y_test)

MDuchi_SDCA_fit <- SDCA_kernel_MDuchi(X = X_train, y = y_train, C = 1, intercept = F, gap_cvg_tol = 1e-5)
mean(predict(MDuchi_SDCA_fit, X_test) == y_test)

Duchi_CVXR_fit <- Duchi_pri_opt(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(Duchi_CVXR_fit, X_test) == y_test)

Duchi_SDCA_fit <- SDCA_kernel_Duchi(X = X_train, y = y_train, C = 1, intercept = F)
mean(predict(Duchi_SDCA_fit, X_test) == y_test)



CS_CVXR_CVXR_fit <- WW_pri_opt(X = X_train, y = y_train, C = 1, intercept = T)


#################### Linear ################################

set.seed(12123)
CV <- function(flds, X_train, y_train, cost, name, intercept){
  flds.acc <- sapply(1:length(flds), function(i){
    X_train_cv <- X_train[flds[[i]], , drop = F]
    y_train_cv <- y_train[flds[[i]]]
    X_test_cv <- X_train[-flds[[i]], , drop = F]
    y_test_cv <- y_train[-flds[[i]]]
    X_train_cv <- scale(X_train_cv, center = T, scale = T)
    X_test_cv <-  scale(X_test_cv, center = attr(X_train_cv,"scaled:center"), scale = attr(X_train_cv,"scaled:scale"))
        
    acc <- NA
    
    try({
      if(name == "WW_CVXR"){
        fit <- WW_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept =  intercept)
      }else if(name == "WW_ksvm"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), kernel = vanilladot(), type = "kbb-svc", C = cost, scaled = F)
      }else if(name == "WW_Pack"){
        return(mean(MSVMPack_interface("WW", X_train_cv, y_train_cv, X_test_cv, y_test_cv, C = cost, kernel = 1, sigma = 0) == y_test_cv))
      }
      else if(name == "CS_CVXR"){
        fit <- CS_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "CS_SDCA"){
        fit <- SDCA_kernel_CS(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "CS_ksvm"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), kernel = vanilladot(), type = "spoc-svc", C = cost, scaled = F)
      }else if(name == "CS_Pack"){
        return(mean(MSVMPack_interface("CS", X_train_cv, y_train_cv, X_test_cv, y_test_cv, C = cost, kernel = 1, sigma = 0) == y_test_cv))
      }else if(name == "Duchi_CVXR"){
        fit <- Duchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      } else if(name == "Duchi_SDCA"){
        fit <- SDCA_kernel_Duchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "MDuchi_CVXR"){
        fit <- MDuchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "MDuchi_SDCA"){
        fit <- SDCA_kernel_MDuchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "LLW_Pack"){
        return(mean(MSVMPack_interface("LLW", X_train_cv, y_train_cv, X_test_cv, y_test_cv, C = cost, kernel = 1, sigma = 0) == y_test_cv))
      }else if(name == "LLW_CVXR"){
        fit <- MDuchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }
      
      pred <- predict(fit,  X_test_cv)
      acc <- mean(pred == y_test_cv)
    })
    return(acc)
  })
  return(mean(flds.acc, na.rm = T))
}

C_list <- 2^seq(3, -4)
set.seed(123)
# C_list <- 2^seq(12, -4)
m <- 6
mu_list <- matrix(rnorm(2*m, mean =0, sd = 10), m)
data <- data_generate(200, sep = 2, v = 10.5^2, m, mu_list)
data_test <- data_generate(1000, sep = 2, v = 10.5^2, m, mu_list)
X_train = data$X 
y_train = data$y
X_test = data_test$X
y_test = data_test$y

# Vehicle <- read.table("~/Desktop/Multiclass Classification/MSVM Code/Vehicle Data/tmp.dat", quote="\"", comment.char="")
# X_dat <- as.matrix(Vehicle[,1:18])
# y_dat <- as.numeric(Vehicle[,19])
rep_n <- 100

predict.ksvm <- getMethod("predict", "ksvm")
cl <- makeCluster(10)

clusterExport(cl=cl, varlist=c("WW_pri_opt", "Duchi_pri_opt","MDuchi_pri_opt", "CS_pri_opt","svm",
                               "LLW_pri_opt", "OVA_pri_opt", "MSVM8_pri_opt","MSVM7_pri_opt","New1_pri_opt","New3_pri_opt","ginv",
                               "solve", "Minimize", "sum_squares", "sum_entries", "vstack", "reshape_expr","CV","predict.msvm",
                               "Variable","Problem","max_entries","min_entries","mvrnorm","Maximize","vec","svm", "CV_kernel", "predict.msvm_kernel", "ksvm", "L_MDuchi", "PL_kernel_MDuchi_Loss", "predict.ksvm",
                               "DL_kernel_MDuchi_Loss", "SDCA_kernel_MDuchi", "PL_kernel_CS_Loss", "DL_kernel_CS_Loss","SDCA_kernel_CS",
                               "PL_kernel_Duchi_Loss", "DL_kernel_Duchi_Loss", "SDCA_kernel_Duchi", "libsvm_kernel_OVA", "rbfdot",
                               "kernelMatrix", "MSVMPack_interface"), envir=environment())

res <- parSapply(cl, 1:rep_n, function(s, C_list, X_train, y_train, X_test, y_test){
  # s_idx <- sample(1:nrow(X_dat), round(nrow(X_dat)*0.7))
  # X_train <- X_dat
  # y_train <- y_dat
  # X_test <- X_dat[-s_idx,]
  # y_test <- y_dat[-s_idx]
  X_train <- scale(X_train, center = T, scale = T)
  X_test <-  scale(X_test, center = attr(X_train,"scaled:center"), scale = attr(X_train,"scaled:scale"))
  
  kfold <- 5
  
  flds <- lapply(1:kfold, function(x){which(x!=((seq(1,length(y_train))%%5)+1))})
  
  
  WW_CVXR_acc <- tryCatch({
    WW_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "WW_CVXR", intercept = T)
    }))
    WW_CVXR_opt_C <- C_list[WW_opt_idx]
    WW_fit <- WW_pri_opt(X = X_train, y = y_train, C = WW_CVXR_opt_C, intercept = T)
    
    mean(predict(WW_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  WW_ksvm_acc <- tryCatch({
    WW_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "WW_ksvm", intercept = T)
    }))
    WW_ksvm_opt_C <- C_list[WW_opt_idx]
    
    WW_fit <-  ksvm(x = X_train, y = as.factor(y_train), kernel =vanilladot(), type = "kbb-svc", C = WW_ksvm_opt_C, scaled = F)
    mean(predict(WW_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  # debugonce(MSVMPack_interface)
  
  WW_Pack_acc <- tryCatch({
    # for(l in 1:length(C_list)){
    #   WW_Pack_opt_idx <-  CV(flds, X_train, y_train, cost = C_list[l], name = "WW_Pack", intercept = T)
    # }
   
    WW_Pack_opt_idx <- which.max( sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "WW_Pack", intercept = T)
    }))
    WW_Pack_opt_C <- C_list[WW_Pack_opt_idx]
    
    mean(MSVMPack_interface("WW", X_train, y_train, X_test, y_test, C = WW_Pack_opt_C, kernel = 1, sigma = 0) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  CS_CVXR_acc <- tryCatch({
    CS_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "CS_CVXR", intercept = F)
    }))
    CS_CVXR_opt_C <- C_list[CS_opt_idx]
    CS_fit <- CS_pri_opt(X = X_train, y = y_train, C = CS_CVXR_opt_C, intercept = F)
    mean(predict(CS_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
 
  CS_ksvm_acc <- tryCatch({
    CS_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "CS_ksvm", intercept = F)
    }))
    CS_ksvm_opt_C <- C_list[CS_opt_idx]
    CS_fit <-  ksvm(x = X_train, y = as.factor(y_train), kernel ="vanilladot", type = "spoc-svc", C = CS_ksvm_opt_C, scaled = F)
    mean(predict(CS_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  CS_Pack_acc <- tryCatch({
    CS_Pack_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "WW_Pack", intercept = T)
    }))
    CS_Pack_opt_C <- C_list[CS_Pack_opt_idx]
    
    mean(MSVMPack_interface("CS", X_train, y_train, X_test, y_test, C = CS_Pack_opt_C, kernel = 1, sigma = 0) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  CS_SDCA_acc <- tryCatch({
    CS_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "CS_SDCA", intercept = F)
    }))
    CS_SDCA_opt_C <- C_list[CS_opt_idx]
    CS_fit <- SDCA_kernel_CS(X = X_train, y = y_train, C = CS_SDCA_opt_C, intercept = F)
    mean(predict(CS_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  Duchi_CVXR_acc <- tryCatch({
    Duchi_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "Duchi_CVXR", intercept = F)
    }))
    Duchi_CVXR_opt_C <- C_list[Duchi_opt_idx]
    Duchi_fit <- Duchi_pri_opt(X = X_train, y = y_train, C = Duchi_CVXR_opt_C, intercept = F)
    mean(predict(Duchi_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  Duchi_SDCA_acc <- tryCatch({
    Duchi_SDCA_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "Duchi_SDCA", intercept = F)
    }))
    Duchi_SDCA_opt_C <- C_list[Duchi_SDCA_opt_idx]
    Duchi_fit <- SDCA_kernel_Duchi(X = X_train, y = y_train, C = Duchi_SDCA_opt_C, intercept = F)
    mean(predict(Duchi_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  MDuchi_CVXR_acc <- tryCatch({
    MDuchi_CVXR_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "MDuchi_CVXR", intercept = F)
    }))
    MDuchi_CVXR_opt_C <- C_list[MDuchi_CVXR_opt_idx]
    MDuchi_fit <- MDuchi_pri_opt(X = X_train, y = y_train, C = MDuchi_CVXR_opt_C, intercept = F)
    mean(predict(MDuchi_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  MDuchi_SDCA_acc <- tryCatch({
    MDuchi_SDCA_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "MDuchi_SDCA", intercept = F)
    }))
    MDuchi_SDCA_opt_C <- C_list[MDuchi_SDCA_opt_idx]
    MDuchi_SDCA_fit <- SDCA_kernel_MDuchi(X = X_train, y = y_train, C = MDuchi_SDCA_opt_C, intercept = F)
    
    mean(predict(MDuchi_SDCA_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
 
  LLW_CVXR_acc <- tryCatch({
    LLW_CVXR_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "LLW_CVXR", intercept = F)
    }))
    LLW_CVXR_opt_C <- C_list[LLW_CVXR_opt_idx]
    LLW_fit <- LLW_pri_opt(X = X_train, y = y_train, C = LLW_CVXR_opt_C, intercept = F)
    mean(predict(LLW_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  
  LLW_Pack_acc <- tryCatch({
    LLW_Pack_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "LLW_Pack", intercept = T)
    }))
    LLW_Pack_opt_C <- C_list[LLW_Pack_opt_idx]
    
    mean(MSVMPack_interface("CS", X_train, y_train, X_test, y_test, C = LLW_Pack_opt_C, kernel = 1, sigma = 0) == y_test)
  }, error = function(e) {
    return(NA)
  })
  # res <- list(WW_acc = WW_acc, Duchi_acc = Duchi_acc, MDuchi_acc = MDuchi_acc, CS_acc = CS_acc, OVA_acc = OVA_acc, LLW_acc = LLW_acc, MSVM8_acc = MSVM8_acc, MSVM7_acc = MSVM7_acc)
  acc_res <-  c(WW_CVXR = WW_CVXR_acc, WW_ksvm = WW_ksvm_acc, CS_CVXR = CS_CVXR_acc, CS_ksvm = CS_ksvm_acc, CS_SDCA = CS_SDCA_acc,
                Duchi_CVXR = Duchi_CVXR_acc, Duchi_SDCA = Duchi_SDCA_acc, MDuchi_CVXR = MDuchi_CVXR_acc, MDuchi_SDCA = MDuchi_SDCA_acc)
  # C_res <- c(OVO_opt_C = OVO_opt_C, WW_opt_C = WW_opt_C, CS_opt_C = CS_opt_C, Duchi_opt_C = Duchi_opt_C, MDuchi_opt_C = MDuchi_opt_C,  OVA_opt_C = OVA_opt_C, LLW_opt_C = LLW_opt_C, MSVM8_opt_C = MSVM8_opt_C)
  
  return(c(acc_res))
}, C_list = C_list, X_train = data$X, y_train = data$y, X_test = data_test$X, y_test = data_test$y)



stopCluster(cl)
rowMeans(res)


#################### Linear II ################################

set.seed(12123)
CV <- function(flds, X_train, y_train, cost, name, intercept){
  flds.acc <- sapply(1:length(flds), function(i){
    X_train_cv <- X_train[flds[[i]], , drop = F]
    y_train_cv <- y_train[flds[[i]]]
    X_test_cv <- X_train[-flds[[i]], , drop = F]
    y_test_cv <- y_train[-flds[[i]]]
    X_train_cv <- scale(X_train_cv, center = T, scale = T)
    X_test_cv <-  scale(X_test_cv, center = attr(X_train_cv,"scaled:center"), scale = attr(X_train_cv,"scaled:scale"))
    
    acc <- NA
    
    try({
      if(name == "WW_CVXR"){
        fit <- WW_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept =  intercept)
      }else if(name == "WW_ksvm"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), kernel = vanilladot(), type = "kbb-svc", C = cost, scaled = F)
      }else if(name == "WW_Pack"){
        return(mean(MSVMPack_interface("WW", X_train_cv, y_train_cv, X_test_cv, y_test_cv, C = cost, kernel = 1, sigma = 0) == y_test_cv))
      }
      else if(name == "CS_CVXR"){
        fit <- CS_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "CS_SDCA"){
        fit <- SDCA_kernel_CS(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "CS_ksvm"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), kernel = vanilladot(), type = "spoc-svc", C = cost, scaled = F)
      }else if(name == "CS_Pack"){
        return(mean(MSVMPack_interface("CS", X_train_cv, y_train_cv, X_test_cv, y_test_cv, C = cost, kernel = 1, sigma = 0) == y_test_cv))
      }else if(name == "Duchi_CVXR"){
        fit <- Duchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      } else if(name == "Duchi_SDCA"){
        fit <- SDCA_kernel_Duchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "MDuchi_CVXR"){
        fit <- MDuchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "MDuchi_SDCA"){
        fit <- SDCA_kernel_MDuchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "LLW_Pack"){
        return(mean(MSVMPack_interface("LLW", X_train_cv, y_train_cv, X_test_cv, y_test_cv, C = cost, kernel = 1, sigma = 0) == y_test_cv))
      }else if(name == "LLW_CVXR"){
        fit <- MDuchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }
      
      pred <- predict(fit,  X_test_cv)
      acc <- mean(pred == y_test_cv)
    })
    return(acc)
  })
  return(mean(flds.acc, na.rm = T))
}

C_list <- 2^seq(3, -4)
set.seed(123)
# C_list <- 2^seq(12, -4)
m <- 6
mu_list <- matrix(rnorm(2*m, mean =0, sd = 10), m)
data <- data_generate(200, sep = 2, v = 10.5^2, m, mu_list)
data_test <- data_generate(1000, sep = 2, v = 10.5^2, m, mu_list)
X_train = data$X 
y_train = data$y
X_test = data_test$X
y_test = data_test$y

# Vehicle <- read.table("~/Desktop/Multiclass Classification/MSVM Code/Vehicle Data/tmp.dat", quote="\"", comment.char="")
# X_dat <- as.matrix(Vehicle[,1:18])
# y_dat <- as.numeric(Vehicle[,19])
rep_n <- 100

predict.ksvm <- getMethod("predict", "ksvm")
cl <- makeCluster(10)

clusterExport(cl=cl, varlist=c("WW_pri_opt", "Duchi_pri_opt","MDuchi_pri_opt", "CS_pri_opt","svm",
                               "LLW_pri_opt", "OVA_pri_opt", "MSVM8_pri_opt","MSVM7_pri_opt","New1_pri_opt","New3_pri_opt","ginv",
                               "solve", "Minimize", "sum_squares", "sum_entries", "vstack", "reshape_expr","CV","predict.msvm",
                               "Variable","Problem","max_entries","min_entries","mvrnorm","Maximize","vec","svm", "CV_kernel", "predict.msvm_kernel", "ksvm", "L_MDuchi", "PL_kernel_MDuchi_Loss", "predict.ksvm",
                               "DL_kernel_MDuchi_Loss", "SDCA_kernel_MDuchi", "PL_kernel_CS_Loss", "DL_kernel_CS_Loss","SDCA_kernel_CS",
                               "PL_kernel_Duchi_Loss", "DL_kernel_Duchi_Loss", "SDCA_kernel_Duchi", "vanilladot", "libsvm_kernel_OVA", "rbfdot",
                               "kernelMatrix", "MSVMPack_interface"), envir=environment())
C=1
res <- parSapply(cl, 1:rep_n, function(s, C, X_train, y_train, X_test, y_test){
  # s_idx <- sample(1:nrow(X_dat), round(nrow(X_dat)*0.7))
  # X_train <- X_dat
  # y_train <- y_dat
  # X_test <- X_dat[-s_idx,]
  # y_test <- y_dat[-s_idx]
  X_train <- scale(X_train, center = T, scale = T)
  X_test <-  scale(X_test, center = attr(X_train,"scaled:center"), scale = attr(X_train,"scaled:scale"))
  
  kfold <- 5
  
  flds <- lapply(1:kfold, function(x){which(x!=((seq(1,length(y_train))%%5)+1))})
  
  
 
  WW_fit <- WW_pri_opt(X = X_train, y = y_train, C = C, intercept = T)
  WW_CVXR_acc <- mean(predict(WW_fit, X_test) == y_test)

  WW_fit <-  ksvm(x = X_train, y = as.factor(y_train), kernel =vanilladot(), type = "kbb-svc", C = C, scaled = F)
  WW_ksvm_acc <- mean(predict(WW_fit, X_test) == y_test)



  CS_fit <- CS_pri_opt(X = X_train, y = y_train, C = C, intercept = F)
  CS_CVXR_acc  <- mean(predict(CS_fit, X_test) == y_test)

  
  CS_fit <-  ksvm(x = X_train, y = as.factor(y_train), kernel = vanilladot(), type = "spoc-svc", C = C, scaled = F)
  CS_ksvm_acc <- mean(predict(CS_fit, X_test) == y_test)

  CS_Pack_acc <-  mean(MSVMPack_interface("CS", X_train, y_train, X_test, y_test, C = C, kernel = 1, sigma = 0) == y_test)

  CS_fit <- SDCA_kernel_CS(X = X_train, y = y_train, C = C, intercept = F)
  CS_SDCA_acc <-   mean(predict(CS_fit, X_test) == y_test)

  Duchi_fit <- Duchi_pri_opt(X = X_train, y = y_train, C = C, intercept = F)
  Duchi_CVXR_acc  <- mean(predict(Duchi_fit, X_test) == y_test)

  Duchi_fit <- SDCA_kernel_Duchi(X = X_train, y = y_train, C = C, intercept = F)
  Duchi_SDCA_acc <- mean(predict(Duchi_fit, X_test) == y_test)

  MDuchi_fit <- MDuchi_pri_opt(X = X_train, y = y_train, C = C, intercept = F)
  MDuchi_CVXR_acc <-  mean(predict(MDuchi_fit, X_test) == y_test)

  MDuchi_SDCA_fit <- SDCA_kernel_MDuchi(X = X_train, y = y_train, C = C, intercept = F)
  MDuchi_SDCA_acc <- mean(predict(MDuchi_SDCA_fit, X_test) == y_test)

  LLW_fit <- LLW_pri_opt(X = X_train, y = y_train, C = C, intercept = F)
  LLW_CVXR_acc <- mean(predict(LLW_fit, X_test) == y_test)

  LLW_Pack_acc <-  mean(MSVMPack_interface("CS", X_train, y_train, X_test, y_test, C = C, kernel = 1, sigma = 0) == y_test)
 
  # res <- list(WW_acc = WW_acc, Duchi_acc = Duchi_acc, MDuchi_acc = MDuchi_acc, CS_acc = CS_acc, OVA_acc = OVA_acc, LLW_acc = LLW_acc, MSVM8_acc = MSVM8_acc, MSVM7_acc = MSVM7_acc)
  acc_res <-  c(WW_CVXR = WW_CVXR_acc, WW_ksvm = WW_ksvm_acc, WW_Pack = WW_Pack_acc,
                CS_CVXR = CS_CVXR_acc, CS_ksvm = CS_ksvm_acc, CS_SDCA = CS_SDCA_acc, CS_Pack =  CS_Pack_acc,
                Duchi_CVXR = Duchi_CVXR_acc, Duchi_SDCA = Duchi_SDCA_acc, 
                MDuchi_CVXR = MDuchi_CVXR_acc, MDuchi_SDCA = MDuchi_SDCA_acc,
                LLW_CVXR = LLW_CVXR_acc, LLW_Pack =  LLW_Pack_acc)
  
  return(c(acc_res))
}, C = 10, X_train = data$X, y_train = data$y, X_test = data_test$X, y_test = data_test$y)


stopCluster(cl)
rowMeans(res)



#################### Kernel ################################
source("~/Desktop/Multiclass Classification/MSVM Code/kernel primary form functions.R")
set.seed(12123)
C_list <- 2^seq(12, -4)


CV <- function(flds, X_train, y_train, cost, name, intercept){
  flds.acc <- sapply(1:length(flds), function(i){
    X_train_cv <- X_train[flds[[i]], , drop = F]
    y_train_cv <- y_train[flds[[i]]]
    X_test_cv <- X_train[-flds[[i]], , drop = F]
    y_test_cv <- y_train[-flds[[i]]]
    X_train_cv <- scale(X_train_cv, center = T, scale = T)
    X_test_cv <-  scale(X_test_cv, center = attr(X_train_cv,"scaled:center"), scale = attr(X_train_cv,"scaled:scale"))
    
    acc <- NA
    
    try({
      if(name == "WW_CVXR"){
        fit <- WW_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept =  intercept)
      }else if(name == "WW_ksvm"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), type = "kbb-svc", C = cost, scaled = F)
      }
      else if(name == "CS_CVXR"){
        fit <- CS_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "CS_SDCA"){
        fit <- SDCA_kernel_CS(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "CS_ksvm"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), type = "spoc-svc", C = cost, scaled = F)
      }
      else if(name == "Duchi_CVXR"){
        fit <- Duchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      } else if(name == "Duchi_SDCA"){
        fit <- SDCA_kernel_Duchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "MDuchi_CVXR"){
        fit <- MDuchi_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }else if(name == "MDuchi_SDCA"){
        fit <- SDCA_kernel_MDuchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept)
      }
      
      pred <- predict(fit,  X_test_cv)
      acc <- mean(pred == y_test_cv)
    })
    return(acc)
  })
  return(mean(flds.acc, na.rm = T))
}

# C_list <- 2^seq(4, -4)




################


data <- data_generate(400, sep = 3, v = 10.5^2, m)
X_dat <- data$X
y_dat <- data$y

rep_n <- 20
predict.ksvm <- getMethod("predict", "ksvm")
cl <- makeCluster(10)
clusterExport(cl=cl, varlist=c("WW_pri_opt", "Duchi_pri_opt","MDuchi_pri_opt", "CS_pri_opt","svm",
                               "LLW_pri_opt", "OVA_pri_opt", "MSVM8_pri_opt","MSVM7_pri_opt","New1_pri_opt","New3_pri_opt","ginv",
                               "solve", "Minimize", "sum_squares", "sum_entries", "vstack", "reshape_expr","CV","predict.msvm",
                               "Variable","Problem","max_entries","min_entries","mvrnorm","Maximize","vec","svm", "CV_kernel", "predict.msvm_kernel", "ksvm", "L_MDuchi", "PL_kernel_MDuchi_Loss", "predict.ksvm",
                               "DL_kernel_MDuchi_Loss", "SDCA_kernel_MDuchi", "PL_kernel_CS_Loss", "DL_kernel_CS_Loss","SDCA_kernel_CS",
                               "PL_kernel_Duchi_Loss", "DL_kernel_Duchi_Loss", "SDCA_kernel_Duchi", "libsvm_kernel_OVA", "rbfdot",
                               "kernelMatrix"), envir=environment())

res <- parSapply(cl, 1:rep_n, function(s, C_list, X_dat, y_dat){
  s_idx <- sample(1:nrow(X_dat), round(nrow(X_dat)*0.7))
  X_train <- X_dat[s_idx,]
  y_train <- y_dat[s_idx]
  X_test <- X_dat[-s_idx,]
  y_test <- y_dat[-s_idx]
  X_train <- scale(X_train, center = T, scale = T)
  X_test <-  scale(X_test, center = attr(X_train,"scaled:center"), scale = attr(X_train,"scaled:scale"))
  
  kfold <- 5
  
  flds <- lapply(1:kfold, function(x){which(x!=((seq(1,length(y_train))%%5)+1))})
  
  
  WW_CVXR_acc <- tryCatch({
    WW_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "WW_CVXR", intercept = T)
    }))
    WW_CVXR_opt_C <- C_list[WW_opt_idx]
    WW_fit <- WW_pri_opt(X = X_train, y = y_train, C = WW_CVXR_opt_C, intercept = F)
    
    mean(predict(WW_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  WW_ksvm_acc <- tryCatch({
    WW_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "WW_ksvm", intercept = T)
    }))
    WW_ksvm_opt_C <- C_list[WW_opt_idx]
    
    WW_fit <-  ksvm(x = X_train, y = as.factor(y_train), type = "kbb-svc", C = WW_ksvm_opt_C, scaled = F)
    mean(predict(WW_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  CS_CVXR_acc <- tryCatch({
    CS_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "CS_CVXR", intercept = F)
    }))
    CS_CVXR_opt_C <- C_list[CS_opt_idx]
    CS_fit <- CS_pri_opt(X = X_train, y = y_train, C = CS_CVXR_opt_C, intercept = F)
    mean(predict(CS_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  CS_ksvm_acc <- tryCatch({
    CS_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "CS_ksvm", intercept = F)
    }))
    CS_ksvm_opt_C <- C_list[CS_opt_idx]
    CS_fit <-  ksvm(x = X_train, y = as.factor(y_train), type = "spoc-svc", C = CS_ksvm_opt_C, scaled = F)
    mean(predict(CS_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  CS_SDCA_acc <- tryCatch({
    CS_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "CS_SDCA", intercept = F)
    }))
    CS_SDCA_opt_C <- C_list[CS_opt_idx]
    CS_fit <- SDCA_kernel_CS(X = X_train, y = y_train, C = CS_SDCA_opt_C, intercept = F)
    mean(predict(CS_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  Duchi_CVXR_acc <- tryCatch({
    Duchi_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "Duchi_CVXR", intercept = F)
    }))
    Duchi_CVXR_opt_C <- C_list[Duchi_opt_idx]
    Duchi_fit <- Duchi_pri_opt(X = X_train, y = y_train, C = Duchi_CVXR_opt_C, intercept = F)
    mean(predict(Duchi_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  Duchi_SDCA_acc <- tryCatch({
    Duchi_SDCA_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "Duchi_SDCA", intercept = F)
    }))
    Duchi_SDCA_opt_C <- C_list[Duchi_SDCA_opt_idx]
    Duchi_fit <- SDCA_kernel_Duchi(X = X_train, y = y_train, C = Duchi_SDCA_opt_C, intercept = F)
    mean(predict(Duchi_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  MDuchi_CVXR_acc <- tryCatch({
    MDuchi_CVXR_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "MDuchi_CVXR", intercept = F)
    }))
    MDuchi_CVXR_opt_C <- C_list[MDuchi_CVXR_opt_idx]
    MDuchi_fit <- MDuchi_pri_opt(X = X_train, y = y_train, C = MDuchi_CVXR_opt_C, intercept = F)
    mean(predict(MDuchi_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  MDuchi_SDCA_acc <- tryCatch({
    MDuchi_SDCA_opt_idx <- which.max(sapply(C_list, function(C) {
      CV(flds, X_train, y_train, cost = C, name = "MDuchi_SDCA", intercept = F)
    }))
    MDuchi_SDCA_opt_C <- C_list[MDuchi_SDCA_opt_idx]
    MDuchi_SDCA_fit <- MDuchi_pri_opt(X = X_train, y = y_train, C = MDuchi_SDCA_opt_C, intercept = F)
    
    mean(predict(MDuchi_SDCA_fit, X_test) == y_test)
  }, error = function(e) {
    return(NA)
  })
  
  
  # res <- list(WW_acc = WW_acc, Duchi_acc = Duchi_acc, MDuchi_acc = MDuchi_acc, CS_acc = CS_acc, OVA_acc = OVA_acc, LLW_acc = LLW_acc, MSVM8_acc = MSVM8_acc, MSVM7_acc = MSVM7_acc)
  acc_res <-  c(WW_CVXR = WW_CVXR_acc, WW_ksvm = WW_ksvm_acc, CS_CVXR = CS_CVXR_acc, CS_ksvm = CS_ksvm_acc, CS_SDCA = CS_SDCA_acc,
                Duchi_CVXR = Duchi_CVXR_acc, Duchi_SDCA = Duchi_SDCA_acc, MDuchi_CVXR = MDuchi_CVXR_acc, MDuchi_SDCA = MDuchi_SDCA_acc)
  # C_res <- c(OVO_opt_C = OVO_opt_C, WW_opt_C = WW_opt_C, CS_opt_C = CS_opt_C, Duchi_opt_C = Duchi_opt_C, MDuchi_opt_C = MDuchi_opt_C,  OVA_opt_C = OVA_opt_C, LLW_opt_C = LLW_opt_C, MSVM8_opt_C = MSVM8_opt_C)
  
  return(c(acc_res))
}, C_list = C_list, X_dat = X_dat, y_dat = y_dat)



stopCluster(cl)


