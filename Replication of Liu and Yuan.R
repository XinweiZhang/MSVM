library(CVXR)
library(MASS)
library(nnet)
library(parallel)
library(CVXfromR)
set.seed(12123)
rm(list = ls())
par(mfrow=c(1,1))
setwd("~/Desktop/Multiclass Classification/MSVM Code")
source("~/Desktop/Multiclass Classification/MSVM Code/primary form functions.R")
source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")

# oracle_data <- four_class_data_generate(10000, sep = .99)
# oracle <- multinom(oracle_data$y ~ oracle_data$X, trace = F)
# mean(predict(oracle, oracle_data$X) == oracle_data$y)


# plot(oracle_data$X[,1], oracle_data$X[,2], col=oracle_data$y, xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "LLW")


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
n_list <- c(100)
# n_list <- c(10,15)

C_list <- 2^seq(15,-16)

rep_n <- 100


res.array <- array(0, dim = c(length(sep_list), length(n_list), 4, rep_n))
i <- j <- 1

cl <- makeCluster(10)
for(i in 1:length(sep_list)){
  
  # oracle <- multinom(oracle_data$y ~ oracle_data$X, trace = F)
  # oracle_acc <- mean(predict(oracle, oracle_data$X) == oracle_data$y)
  
  for( j in 1:length(n_list)){
    clusterExport(cl=cl, varlist=c("WW_pri_opt", "Duchi_pri_opt","MDuchi_pri_opt", "CS_pri_opt",
                                   "LLW_pri_opt", "OVA_pri_opt", "MSVM8_pri_opt","MSVM7_pri_opt","New1_pri_opt","New3_pri_opt","ginv",
                                   "solve", "Minimize", "sum_squares", "sum_entries", "vstack", "reshape_expr",
                                   "Variable","Problem","max_entries","min_entries","pred","data_generate","mvrnorm","Maximize","vec"), envir=environment())
    acc.res <- parSapply(cl, 1:rep_n, function(s, C_list, n, sep){
      # n = n_list[j]; sep = sep_list[i];
      oracle_data <- data_generate(10^5, sep = sep)
      train_data <- data_generate(n, sep = sep)
      val_data <- data_generate(n, sep = sep)
      WW_acc <- tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          WW_w <- WW_pri_opt(train_data$X,train_data$y,C)
          mean(pred(val_data$X,WW_w) == val_data$y)
        }))
        WW_w <- WW_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
        mean(pred(oracle_data$X,WW_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      LLW_acc <-  tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          LLW_w <- LLW_pri_opt(train_data$X,train_data$y,C)
          mean(pred(val_data$X,LLW_w) == val_data$y)
        }))
        LLW_w <- LLW_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
        mean(pred(oracle_data$X,LLW_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      OVA_acc <-tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          OVA_w <- OVA_pri_opt(train_data$X,train_data$y,C)
          mean(pred(val_data$X,OVA_w) == val_data$y)
        }))
        OVA_w <- OVA_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
        mean(pred(oracle_data$X,OVA_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
     
      MSVM8_acc <-  tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          MSVM8_w <- MSVM8_pri_opt(train_data$X,train_data$y,C)
          mean(pred(val_data$X,MSVM8_w) == val_data$y)
        }))
        MSVM8_w <- MSVM8_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
        mean(pred(oracle_data$X,MSVM8_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      
      # res <- list(WW_acc = WW_acc, Duchi_acc = Duchi_acc, MDuchi_acc = MDuchi_acc, CS_acc = CS_acc, OVA_acc = OVA_acc, LLW_acc = LLW_acc, MSVM8_acc = MSVM8_acc, MSVM7_acc = MSVM7_acc)
      res <- c(WW_acc = WW_acc, OVA_acc = OVA_acc, LLW_acc = LLW_acc, MSVM8_acc = MSVM8_acc)
      
      return(res)
    }, C_list=C_list, n = n_list[j], sep = sep_list[i])
    
    res.array[i,j,,] <- acc.res
  }
}
stopCluster(cl)

# res.array


summary.res.arry <- apply(res.array, c(1,2,3), function(x){mean(x,na.rm = T)})
dimnames(summary.res.arry)[[1]] <- c("1")
dimnames(summary.res.arry)[[3]] <-  c("WW","OVA","LLW","MSVM8")

round(summary.res.arry[,1,]*100,2)
      
save(sep_list, rep_n, n_list, res.array, summary.res.arry, file = 
       paste("~/Desktop/Multiclass Classification/MSVM Code/Raplication of Liu and Yuan, Plan I, sep= 1, n = 50, rep=",rep_n,".Rdata",sep=""))

