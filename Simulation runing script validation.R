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
# source("~/Desktop/Multiclass Classification/MSVM Code/Dual form functions.R")
m <- 3
# oracle_data <- data_generate(10000, sep = 3.9)
# oracle <- multinom(oracle_data$y ~ oracle_data$X, trace = F)
# mean(predict(oracle, oracle_data$X) == oracle_data$y)


# plot(oracle_data$X[1:100,],  col='red', xlim = c(-10,10), ylim=c(-10,10), xlab = "X1", ylab = "X2", main = "LLW")
# 
# points(oracle_data$X[101:200,], col='green')
# 
# points(oracle_data$X[201:300,], col='blue')

############################################
### Plan I Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95, .99
# sep_list <- c(.36, .45, .55, .65, .76, .88, 1.01, 1.17,  1.37, 1.66, 2.30)
  
#### Plan II Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95, .99
sep_list <- c(.59, .68, .77, .86, .97, 1.08, 1.22, 1.38, 1.60, 1.95,3.9)

#### Plan III Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95, .99
# sep_list <- c(.26, .44, .62, .80, 1, 1.2, 1.45, 1.72, 2.1, 2.6, 3.7)

#### Plan IV Separation: .5, .55, .60, .65, .70, .75, .80, .85, 90, .95, .99
# sep_list <- c(.63, .73, .83, .93, 1.04, 1.17, 1.31, 1.47, 1.68, 2.01, 3.51)

# sep_list <- c(.93, 1.13)
n_list <- c(4, 6, 8)
# n_list <- c(4,6)

C_list <- 2^seq(4,-4)

# C_list <- 2^seq(-4,4)
# C_list
rep_n <- 10


# 
# i <- j <- 1
res.array <- array(0, dim = c(length(sep_list), length(n_list), 8, rep_n))

cl <- makeCluster(10)
for(i in 1:length(sep_list)){
  oracle_data <- data_generate(10000, sep = sep_list[i])
  oracle <- multinom(oracle_data$y ~ oracle_data$X, trace = F)
  oracle_acc <- mean(predict(oracle, oracle_data$X) == oracle_data$y)
  
  
  for( j in 1:length(n_list)){
    clusterExport(cl=cl, varlist=c("WW_pri_opt", "Duchi_pri_opt","MDuchi_pri_opt", "CS_pri_opt",
                                   "LLW_pri_opt", "OVA_pri_opt", "MSVM8_pri_opt","MSVM7_pri_opt","Duchi_dual_opt","ginv",
                                   "solve", "Minimize", "sum_squares", "sum_entries", "vstack", "reshape_expr",
                                   "Variable","Problem","max_entries","min_entries","pred","data_generate","mvrnorm","Maximize","vec"), envir=environment())
    acc.res <- parSapply(cl, 1:rep_n, function(s, C_list, oracle_data, n, sep){
      # n = n_list[j]; sep = sep_list[i];
      train_data <- data_generate(n, sep = sep)
      test_data <- data_generate(500, sep = sep)
      WW_acc <- tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          WW_w <- WW_pri_opt(train_data$X,train_data$y,C)
          mean(pred(test_data$X,WW_w) == test_data$y)
        }))
        WW_w <- WW_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
        mean(pred(oracle_data$X,WW_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      CS_acc <-tryCatch({
          opt_idx <- which.max(sapply(C_list, function(C){
            CS_w <- CS_pri_opt(train_data$X,train_data$y,C)
            mean(pred(test_data$X,CS_w) == test_data$y)
          }))
          CS_w <- CS_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
          mean(pred(oracle_data$X,CS_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      

      Duchi_acc <- tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          Duchi_w <- Duchi_dual_opt(train_data$X,train_data$y,C)
          mean(pred(test_data$X,Duchi_w) == test_data$y)
        }))
        Duchi_w <- Duchi_dual_opt(train_data$X,train_data$y,C_list[opt_idx])
        mean(pred(oracle_data$X,Duchi_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })

      MDuchi_acc <-  tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
            MDuchi_w <- MDuchi_pri_opt(train_data$X,train_data$y,C)
            mean(pred(test_data$X,MDuchi_w) == test_data$y)
          }))
          MDuchi_w <- MDuchi_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
          mean(pred(oracle_data$X,MDuchi_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      
      LLW_acc <-  tryCatch({
          opt_idx <- which.max(sapply(C_list, function(C){
            LLW_w <- LLW_pri_opt(train_data$X,train_data$y,C)
            mean(pred(test_data$X,LLW_w) == test_data$y)
          }))
          LLW_w <- LLW_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
          mean(pred(oracle_data$X,LLW_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      
      OVA_acc <-tryCatch({
          opt_idx <- which.max(sapply(C_list, function(C){
            OVA_w <- OVA_pri_opt(train_data$X,train_data$y,C)
            mean(pred(test_data$X,OVA_w) == test_data$y)
          }))
          OVA_w <- OVA_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
          mean(pred(oracle_data$X,OVA_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      MSVM7_acc <- tryCatch({
          opt_idx <- which.max(sapply(C_list, function(C){
            MSVM7_w <- MSVM7_pri_opt(train_data$X,train_data$y,C)
            mean(pred(test_data$X,MSVM7_w,"s") == test_data$y)
          }))
          MSVM7_w <- MSVM7_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
          mean(pred(oracle_data$X,MSVM7_w,"s") == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
      
      MSVM8_acc <-  tryCatch({
        opt_idx <- which.max(sapply(C_list, function(C){
          MSVM8_w <- MSVM8_pri_opt(train_data$X,train_data$y,C)
          mean(pred(test_data$X,MSVM8_w) == test_data$y)
        }))
        MSVM8_w <- MSVM8_pri_opt(train_data$X,train_data$y,C_list[opt_idx])
        mean(pred(oracle_data$X,MSVM8_w) == oracle_data$y)
      }, error = function(e) {
        return(NA)
      })
           
      # res <- list(WW_acc = WW_acc, Duchi_acc = Duchi_acc, MDuchi_acc = MDuchi_acc, CS_acc = CS_acc, OVA_acc = OVA_acc, LLW_acc = LLW_acc, MSVM8_acc = MSVM8_acc, MSVM7_acc = MSVM7_acc)
      res <- c(WW_acc = WW_acc, CS_acc = CS_acc, Duchi_acc = Duchi_acc, MDuchi_acc = MDuchi_acc, OVA_acc = OVA_acc, LLW_acc = LLW_acc, MSVM7_acc = MSVM7_acc, MSVM8_acc = MSVM8_acc)
      
      return(res)
    }, C_list=C_list, oracle_data = oracle_data, n = n_list[j], sep = sep_list[i])
    
    res.array[i,j,,] <- acc.res
  }
}
stopCluster(cl)

res.array


summary.res.arry <- apply(res.array, c(1,2,3), function(x){mean(x,na.rm = T)})
dimnames(summary.res.arry)[[1]] <- c(".5", ".55", ".6", ".65", ".7", ".75", ".8", ".85", ".9", ".95", ".99")
dimnames(summary.res.arry)[[3]] <-  c("WW","CS","Duchi","MDuchi","OVA","LLW","MSVM7","MSVM8")

summary.res.arry[,1,]
summary.res.arry[,2,]
summary.res.arry[,3,]

save(sep_list, rep_n, n_list, res.array, summary.res.arry, file = paste("~/Desktop/Multiclass Classification/MSVM Code/MSVM simulation (Plan II, val = 500), sep=",length(sep_list),", n_list=",length(n_list),", rep=",rep_n,".Rdata",sep=""))

