library(CVXR)
library(MASS)
library(nnet)
library(parallel)
library(CVXfromR)
library(kernlab)
set.seed(12123)
rm(list = ls())
par(mfrow=c(1,1))
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/kernel primary form functions.R")
# source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")

# oracle_data <- three_class_data_generate(100, sep = .99)

# oracle <- multinom(oracle_data$y ~ oracle_data$X, trace = F)
# mean(predict(oracle, oracle_data$X) == oracle_data$y)

# plot(oracle_data$X[,1], oracle_data$X     [,2], col=oracle_data$y, xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "LLW")


############################################
## Plan I Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95
sep_list <- c(1)

#### Plan II Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95
# sep_list <- c(.71, .81, .91, .99, 1.11, 1.23, 1.36, 1.53, 1.73, 2.02)

#### Plan III Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95, .99
# sep_list <- c(.26, .44, .62, .80, 1, 1.2, 1.45, 1.72, 2.1, 2.6, 3.7)

#### Plan IV Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95, .99  
# sep_list <- c(.63, .73, .83, .93, 1.04, 1.17, 1.31, 1.47, 1.68, 2.01, 3.51)

# sep_list <- c(.93, 1.13)
n_list <- c(50, 100)
# n_list <- c(10,15)

C_list <- 2^seq(-15,16)
C_list <- 1/(50*C_list)

# C_list <- 2^seq(16,15)
# C=1

rep_n <- 100


res.array <- array(0, dim = c(length(sep_list), length(n_list), 13, rep_n))
i <- j <- 1

cl <- makeCluster(10)
for(i in 1:length(sep_list)){
  
  # oracle <- multinom(oracle_data$y ~ oracle_data$X, trace = F)
  # oracle_acc <- mean(predict(oracle, oracle_data$X) == oracle_data$y)
  
  for( j in 1:length(n_list)){
    clusterExport(cl=cl, varlist=c("WW_kernel_pri_opt", "Duchi_kernel_pri_opt","MDuchi_kernel_pri_opt", "CS_kernel_pri_opt",
                                   "LLW_kernel_pri_opt", "OVA_kernel_pri_opt", "MSVM8_kernel_pri_opt","MSVM7_kernel_pri_opt","New1_kernel_pri_opt","New3_kernel_pri_opt","ginv",
                                   "solve", "Minimize", "sum_squares", "quad_form", "sum_entries", "vstack", "reshape_expr","polydot","kernelMatrix","quad_form",
                                   "Variable","Problem","max_entries","min_entries","predict.msvm_kernel","three_class_data_generate","data_generate","mvrnorm","Maximize","vec","rbfdot"), envir=environment())
    acc.res <- parSapply(cl, 1:rep_n, function(s, C_list, n, sep){
      # n = n_list[j]; sep = sep_list[i];
      ##############NonLinear Data##########################  ##
      oracle_data <- three_class_data_generate(10^5, sep = sep)
      train_data <- three_class_data_generate(n, sep = sep)
      val_data <- three_class_data_generate(n, sep = sep)
      # sigma <- median(c(apply(train_data$X[train_data$y==1,], MARGIN = 1, function(x,y){sqrt(rowSums((matrix(1,nrow = nrow(y))%*%x-y)^2))}, y = train_data$X[train_data$y==2,]),
      #   apply(train_data$X[train_data$y==1,], MARGIN = 1, function(x,y){sqrt(rowSums((matrix(1,nrow = nrow(y))%*%x-y)^2))}, y = train_data$X[train_data$y==3,]),
      #   apply(train_data$X[train_data$y==1,], MARGIN = 1, function(x,y){sqrt(rowSums((matrix(1,nrow = nrow(y))%*%x-y)^2))}, y = train_data$X[train_data$y==3,])))
      # kernel = rbfdot(1/sigma^2)
      kernel = polydot(degree = 2, scale = 1, offset = 1)
      
      ##############Linear Data############################
       # oracle_data <- data_generate(10^5, sep = sep)
       # train_data <- data_generate(n, sep = sep)
       # val_data <- data_generate(n, sep = sep)

      WW_acc <- tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          res <- tryCatch({
          WW_fit <- WW_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
          mean(predict(WW_fit, val_data$X) == val_data$y)}, error = function(e){
            return(NA)
          })
          return(res)
        }))
        WW_fit <- WW_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx], kernel)
        mean(predict(WW_fit, oracle_data$X) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      CS_acc <- tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          res <- tryCatch({
            CS_fit <- CS_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            mean(predict(CS_fit, val_data$X) == val_data$y)}, error = function(e){
              return(NA)
            })
          return(res)
        }))
        CS_fit <- CS_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx], kernel)
        mean(predict(CS_fit, oracle_data$X) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      Duchi_acc <- tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          res <- tryCatch({
            Duchi_fit <- Duchi_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            mean(predict(Duchi_fit, val_data$X) == val_data$y)}, error = function(e){
              return(NA)
            })
          return(res)
        }))
        Duchi_fit <- Duchi_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx], kernel)
        mean(predict(Duchi_fit, oracle_data$X) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })

      MDuchi_acc <- tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          res <- tryCatch({
            MDuchi_fit <- MDuchi_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            mean(predict(MDuchi_fit, val_data$X) == val_data$y)}, error = function(e){
              return(NA)
            })
          return(res)
        }))
        MDuchi_fit <- MDuchi_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx], kernel)
        mean(predict(MDuchi_fit, oracle_data$X) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      New1_acc <- tryCatch({
        opt_idx <- apply(sapply(C_list, function(C){
          res <- tryCatch({
            New1_fit <- New1_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            c(mean(predict(New1_fit, val_data$X, rule = "dagger_New1") == val_data$y), 
              mean(predict(New1_fit, val_data$X, rule = "0_max") == val_data$y))}, error = function(e){
                return(c(NA,NA))
              })
          return(res)
        }), MARGIN = 1, which.max)
        res1 <- tryCatch({ New1_fit_1 <- New1_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx[1]], kernel)
        mean(predict(New1_fit_1, oracle_data$X, rule = "dagger_New1") == oracle_data$y)}, error=function(e){
          return(NA)
        })
        res2 <- tryCatch({New1_fit_2 <- New1_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx[2]], kernel)
         mean(predict(New1_fit_2, oracle_data$X, rule = "0_max") == oracle_data$y)}, error=function(e){
          return(NA)
        })
        c(res1, res2)
      }, error = function(e) {
        return(c(NA,NA))
      })
      
      New3_acc <- tryCatch({
        opt_idx <- apply(sapply(C_list, function(C){
          res <- tryCatch({
            New3_fit <- New3_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            c(mean(predict(New3_fit, val_data$X, rule = "dagger_New1") == val_data$y), 
              mean(predict(New3_fit, val_data$X, rule = "0_max") == val_data$y))}, error = function(e){
                return(c(NA,NA))
              })
          return(res)
        }), MARGIN = 1, which.max)
        res1 <- tryCatch({ New3_fit_1 <- New3_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx[1]], kernel)
        mean(predict(New3_fit_1, oracle_data$X, rule = "dagger_New1") == oracle_data$y)}, error=function(e){
          return(NA)
        })
        res2 <- tryCatch({New3_fit_2 <- New3_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx[2]], kernel)
        mean(predict(New3_fit_2, oracle_data$X, rule = "0_max") == oracle_data$y)}, error=function(e){
          return(NA)
        })
        c(res1, res2)
      }, error = function(e) {
        return(c(NA,NA))
      })
      
      LLW_acc <-  tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          res <- tryCatch({
            LLW_fit <- LLW_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            mean(predict(LLW_fit, val_data$X) == val_data$y)}, error = function(e){
              return(NA)
            })
          return(res)
        }))
         LLW_fit <- LLW_kernel_pri_opt(train_data$X,train_data$y,C_list[opt_idx], kernel)
          mean(predict(LLW_fit, oracle_data$X) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      OVA_acc <-tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          res <- tryCatch({
            OVA_fit <- OVA_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            mean(predict(OVA_fit, val_data$X) == val_data$y)}, error = function(e){
              return(NA)
            })
          return(res)
        }))
        OVA_fit <- OVA_kernel_pri_opt(train_data$X,train_data$y,C_list[opt_idx], kernel)
        mean(predict(OVA_fit, oracle_data$X) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      MSVM7_acc <- tryCatch({
        opt_idx <- apply(sapply(C_list, function(C){
          res <- tryCatch({
            MSVM7_fit <- MSVM7_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
            c(mean(predict(MSVM7_fit, val_data$X, rule = "dagger_MSVM7") == val_data$y), 
              mean(predict(MSVM7_fit, val_data$X, rule = "0_max") == val_data$y))}, error = function(e){
                return(c(NA,NA))
              })
          return(res)
        }), MARGIN = 1, which.max)
        
        MSVM7_fit_1 <- MSVM7_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx[1]], kernel)
        MSVM7_fit_2 <- MSVM7_kernel_pri_opt(train_data$X,train_data$y, C_list[opt_idx[2]], kernel)
        c(mean(predict(MSVM7_fit_1, oracle_data$X, rule = "dagger_MSVM7") == oracle_data$y),
          mean(predict(MSVM7_fit_2, oracle_data$X, rule = "0_max") == oracle_data$y))
      }, error = function(e) {
        return(c(NA,NA))
      })
      
      MSVM8_acc <-  tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          tryCatch({
          MSVM8_fit <- MSVM8_kernel_pri_opt(train_data$X,train_data$y,C, kernel)
          mean(predict(MSVM8_fit, val_data$X) == val_data$y)},
          error =function(e){
            return(NA)
            })
        }))
        MSVM8_fit <- MSVM8_kernel_pri_opt(train_data$X,train_data$y,C_list[opt_idx], kernel)
        mean(predict(MSVM8_fit, oracle_data$X) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
          
      res <-
        c(
          WW_acc = WW_acc,
          CS_acc = CS_acc,
          Duchi_acc = Duchi_acc,
          MDuchi_acc = MDuchi_acc,
          New1_acc = New1_acc,
          New3_acc = New3_acc,
          OVA_acc = OVA_acc,
          LLW_acc = LLW_acc,
          MSVM7_acc = MSVM7_acc,
          MSVM8_acc = MSVM8_acc
        )
      return(res)
    }, C_list=C_list, n = n_list[j], sep = sep_list[i])
    
    res.array[i,j,,] <- acc.res
  }
}
stopCluster(cl)
# 
# res.array <- array(0, dim = c(length(sep_list), length(n_list), 13, 98))
# res.array[1,1,,] <- do.call(rbind,acc.res[-c(19,65)])

# res.array
  
summary.res.arry <- apply(res.array, c(1,2,3), function(x){median(x,na.rm = T)})
dimnames(summary.res.arry)[[1]] <- c("1")
dimnames(summary.res.arry)[[3]] <-  c("WW", "CS", "Duchi", "MDuchi", "New1,d","New1,0", "New3,d","New3,0", "OVA", "LLW", "MSVM7,d","MSVM7,0", "MSVM8")

round((1-summary.res.arry[,,])*100,2)

round(summary.res.arry[,2,]*100,2)

save(sep_list, rep_n, n_list, res.array, summary.res.arry, file = 
       paste("~/Desktop/Multiclass Classification/MSVM Code/Replication Poly Kernel of Liu and Yuan, Nonlinear Data, sep= 1, rep=",rep_n,".Rdata",sep=""))

##100: OVA, LLW, MSVM8  68.83 67.72 68.11
##100 - c(68.83, 67.72, 68.11)
