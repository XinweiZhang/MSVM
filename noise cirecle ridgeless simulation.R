library(CVXR)
library(MASS)
library(nnet)
library(parallel)
library(kernlab)
set.seed(12344)
rm(list = ls())
par(mfrow=c(1,1))
setwd("~/Desktop/Multiclass Classification/MSVM Code")

source("primary form functions.R")

# 
# generate_noise_circle <- function(n, noise_level = 0.95){
#   draw_X <- function(y){
#     t <- sapply(y, function(z){
#       if(z==1){
#         runif(1, min = 0, max = 5)
#       }else if(z==2){
#         runif(1, min = 5, max = 11)
#       }else if(z==3){
#         runif(1, min = 11, max = 20)
#       }
#     })
# 
#     X = t(sapply(t, function(x){c(cos(x*pi/10), sin(x*pi/10))}))
#     return(X)
#   }
# 
#   generate_asy_noise_y <- function(y, noise_level){
#     idx <- sample(length(y), size = length(y)*noise_level, replace = F)
#     y_fake <- sapply(idx, function(x){
#       if(y[x]==1){
#         sample(c(1,2,3), size =1, prob = c(11/27,11/27,5/27))
#       }else if(y[x]==2){
#         sample(c(1,2,3), size =1, prob = c(5/27, 11/27, 11/27))
#       }else if(y[x]==3){
#         sample(c(1,2,3), size =1, prob = c(11/27, 5/27, 11/27))
#       }
#     })
# 
#     y[idx] <- y_fake
#     return(y)
#   }
# 
#   y <- sample(c(1,2,3), size = n, replace = T)
#   X <- draw_X(y)
#   y_noise <- generate_asy_noise_y(y, noise_level=noise_level)
#   return(list(X=X,y=y_noise))
# }



generate_noise_circle <- function(n){

  t <- sample(c(1,2,3), size = n, replace = T)
  dat <- t(sapply(t, function(z){
    # if(z==1){
    #   xz <- runif(1, min = 0, max = 3)
    # }else if(z==2){
    #   xz <- runif(1, min = 3, max = 11)
    # }else if(z==3){
    #   xz <- runif(1, min = 11, max = 20)
    # }
    # x = c(cos(xz*pi/10), sin(xz*pi/10))
    # if(z==1){
    #   xz <- runif(1, min = 0, max = .2)
    # }else if(z==2){
    #   xz <- runif(1, min = 7, max = 7.2)
    # }else if(z==3){
    #   xz <- runif(1, min = 14, max = 14.2)
    # }
    # x = c(cos(xz*pi/10.5), sin(xz*pi/10.5))
    if(z==1){
      xz <- runif(1, min = 0, max = 7)
    }else if(z==2){
      xz <- runif(1, min = 7, max = 14)
    }else if(z==3){
      xz <- runif(1, min = 14, max = 21)
    }
    x = c(cos(xz*pi/10.5), sin(xz*pi/10.5))
    # if(z==1){
    #   y = sample(c(1,2,3), size=1,prob= c(.5, 1/3, 1/6))
    #   # y = sample(c(1,2,3), size=1,prob= c(.9, 0, .1))
    # }else if(z==2){
    #   y = sample(c(1,2,3), size=1,prob= c(1/3, .5, 1/6))
    #   # y = sample(c(1,2,3), size=1,prob= c(.2, .7, .1))
    # }else if(z==3){
    #   y = sample(c(1,2,3), size=1,prob= c(0, 0, 1))
    #   # y = sample(c(1,2,3), size=1,prob= c(1/6, 1/3, .5))
    # }
    # if(z==1){
    #   y = sample(c(1,2,3), size=1,prob= c(.45, .4, .15))
    # }else if(z==2){
    #   y = sample(c(1,2,3), size=1,prob= c(.15, .45, .4))
    # }else if(z==3){
    #   y = sample(c(1,2,3), size=1,prob= c(.4, .15, .45))
    # }

    if(z==1){
      y = sample(c(1,2,3), size=1,prob= c(.7, .2, .1))
    }else if(z==2){
      y = sample(c(1,2,3), size=1,prob= c(.1, .7, .2))
    }else if(z==3){
      y = sample(c(1,2,3), size=1,prob= c(.2, .1, .7))
    }
    return(c(x,y))
  }))
  X <- dat[,1:2]
  y <- dat[,3]
  list(X=X,y=y, oracle_acc = mean(t==y))
}


generate_noise_circle_5class <- function(n){

  t <- sample(c(1,2,3,4), size = n, replace = T)
  dat <- t(sapply(t, function(z){
    if(z==1){
      xz <- runif(1, min = 0, max = 1)
    }else if(z==2){
      xz <- runif(1, min = 1, max = 5)
    }else if(z==3){
      xz <- runif(1, min = 5, max = 14)
    }else if(z==4){
      xz <- runif(1, min = 14, max = 20)
    }
    x = c(cos(xz*pi/10), sin(xz*pi/10))

    if(z==1){
      # y = sample(seq(1,4), size=1, prob= c(.5, .2, .25,.15))
      y = sample(seq(1,4), size=1, prob= c(.8, .05, .15,.05))
    }else if(z==2){
      # y = sample(seq(1,4), size=1, prob= c(.15, .4, .25, .2))
      y = sample(seq(1,4), size=1, prob= c(.05, .8, .15,.05))
    }else if(z==3){
      # y = sample(seq(1,4), size=1, prob= c(.05, .2, .6 , 0.15))
      y = sample(seq(1,4), size=1, prob= c(.15, .05, .8,.05))
    }else if(z==4){
      # y = sample(seq(1,4), size=1, prob= c(.15, .25, .15,  .45 ))
      y = sample(seq(1,4), size=1, prob= c(.15, .05, .05,.8))
    }
    # if(z==1){
    #   y = sample(c(1,2,3), size=1,prob= c(.8, .1, .1))
    # }else if(z==2){
    #   y = sample(c(1,2,3), size=1,prob= c(.3, .4, .3))
    # }else if(z==3){
    #   y = sample(c(1,2,3), size=1,prob= c(.15, .40, .45))
    # }

    # if(z==1){
    #   y = sample(c(1,2,3), size=1,prob= c(.45, .15, .40))
    # }else if(z==2){
    #   y = sample(c(1,2,3), size=1,prob= c(.40, .45, .15))
    # }else if(z==3){
    #   y = sample(c(1,2,3), size=1,prob= c(.15, .40, .45))
    # }
    return(c(x,y))
  }))
  X <- dat[,1:2]
  y <- dat[,3]
  list(X=X,y=y, oracle_acc = mean(t==y))
}




############################################

n_list <- c(500,1000,1500,2000)

rep_n <- 100

# generate_noise_circle_5class(100000)$oracle_acc

res.array <- array(0, dim = c(length(n_list), 13, rep_n))
i <- j <- 1

cl <- makeCluster(10)
for( j in 1:length(n_list)){
  clusterExport(cl=cl, varlist=c("WW_pri_opt", "Duchi_pri_opt","MDuchi_pri_opt", "CS_pri_opt","generate_noise_circle_5class","generate_noise_circle",
                                 "LLW_pri_opt", "OVA_pri_opt", "MSVM8_pri_opt","MSVM7_pri_opt","New1_pri_opt","New3_pri_opt","ginv",
                                 "solve", "Minimize", "sum_squares", "quad_form", "sum_entries", "vstack", "reshape_expr","polydot","vanilladot","kernelMatrix","quad_form",
                                 "Variable","Problem","max_entries","min_entries","predict.msvm","mvrnorm","Maximize","vec","rbfdot"), envir=environment())
  acc.res <- parSapply(cl, 1:rep_n, function(s, n){
    # n = n_list[j];

    #############Linear Data############################
    # oracle_data <- generate_noise_circle_5class(10^5)
    # train_data <- generate_noise_circle_5class(n)
    # # val_data <- generate_noise_circle(n)
    
    oracle_data <- generate_noise_circle(10^5)
    train_data <- generate_noise_circle(n)
    # val_data <- generate_noise_circle(n)

    
    WW_acc <- tryCatch({
      # opt_idx <- which.max(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     WW_fit <- WW_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     mean(predict(WW_fit, val_data$X) == val_data$y)}, error = function(e){
      #       return(NA)
      #     })
      #   return(res)
      # }))
      WW_fit <- WW_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(WW_fit, oracle_data$X) == oracle_data$y)
    }, error = function(e) {
      return(NA)
    })
    
    CS_acc <- tryCatch({
      # opt_idx <- which.max(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     CS_fit <- CS_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     mean(predict(CS_fit, val_data$X) == val_data$y)}, error = function(e){
      #       return(NA)
      #     })
      #   return(res)
      # }))
      CS_fit <- CS_pri_opt(train_data$X,train_data$y,  C=1, lambda=0, intercept = F)
      mean(predict(CS_fit, oracle_data$X) == oracle_data$y)
    }, error = function(e) {
      return(NA)
    })
    
    Duchi_acc <- tryCatch({
      # opt_idx <- which.max(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     Duchi_fit <- Duchi_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     mean(predict(Duchi_fit, val_data$X) == val_data$y)}, error = function(e){
      #       return(NA)
      #     })
      #   return(res)
      # }))
      Duchi_fit <- Duchi_pri_opt(train_data$X,train_data$y,  C=1, lambda=0, intercept = F)
      mean(predict(Duchi_fit, oracle_data$X) == oracle_data$y)
    }, error = function(e) {
      return(NA)
    })
    
    MDuchi_acc <- tryCatch({
      # opt_idx <- which.max(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     MDuchi_fit <- MDuchi_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     mean(predict(MDuchi_fit, val_data$X) == val_data$y)}, error = function(e){
      #       return(NA)
      #     })
      #   return(res)
      # }))
      MDuchi_fit <- MDuchi_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(MDuchi_fit, oracle_data$X) == oracle_data$y)
    }, error = function(e) {
      return(NA)
    })
    
    New1_acc <- tryCatch({
      # opt_idx <- apply(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     New1_fit <- New1_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     c(mean(predict(New1_fit, val_data$X, rule = "dagger") == val_data$y), 
      #       mean(predict(New1_fit, val_data$X, rule = "0_max") == val_data$y))}, error = function(e){
      #         return(c(NA,NA))
      #       })
      #   return(res)
      # }), MARGIN = 1, which.max)
      res1 <- tryCatch({ New1_fit_1 <- New1_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(New1_fit_1, oracle_data$X, rule = "dagger") == oracle_data$y)}, error=function(e){
        return(NA)
      })
      res2 <- tryCatch({New1_fit_2 <- New1_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(New1_fit_2, oracle_data$X, rule = "0_max") == oracle_data$y)}, error=function(e){
        return(NA)
      })
      c(res1, res2)
    }, error = function(e) {
      return(c(NA,NA))
    })
    
    New3_acc <- tryCatch({
      # opt_idx <- apply(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     New3_fit <- New3_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     c(mean(predict(New3_fit, val_data$X, rule = "dagger") == val_data$y), 
      #       mean(predict(New3_fit, val_data$X, rule = "0_max") == val_data$y))}, error = function(e){
      #         return(c(NA,NA))
      #       })
      #   return(res)
      # }), MARGIN = 1, which.max)
      res1 <- tryCatch({ New3_fit_1 <- New3_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(New3_fit_1, oracle_data$X, rule = "dagger") == oracle_data$y)}, error=function(e){
        return(NA)
      })
      res2 <- tryCatch({New3_fit_2 <- New3_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(New3_fit_2, oracle_data$X, rule = "0_max") == oracle_data$y)}, error=function(e){
        return(NA)
      })
      c(res1, res2)
    }, error = function(e) {
      return(c(NA,NA))
    })
    
    LLW_acc <-  tryCatch({
      # opt_idx <- which.max(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     LLW_fit <- LLW_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     mean(predict(LLW_fit, val_data$X) == val_data$y)}, error = function(e){
      #       return(NA)
      #     })
      #   return(res)
      # }))
      LLW_fit <- LLW_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(LLW_fit, oracle_data$X) == oracle_data$y)
    }, error = function(e) {
      return(NA)
    })
    
    OVA_acc <-tryCatch({
      # opt_idx <- which.max(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     OVA_fit <- OVA_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     mean(predict(OVA_fit, val_data$X) == val_data$y)}, error = function(e){
      #       return(NA)
      #     })
      #   return(res)
      # }))
      OVA_fit <- OVA_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(OVA_fit, oracle_data$X) == oracle_data$y)
    }, error = function(e) {
      return(NA)
    })
    
    MSVM7_acc <- tryCatch({
      # opt_idx <- apply(sapply(C_list, function(C){
      #   res <- tryCatch({
      #     MSVM7_fit <- MSVM7_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     c(mean(predict(MSVM7_fit, val_data$X, rule = "dagger") == val_data$y), 
      #       mean(predict(MSVM7_fit, val_data$X, rule = "0_max") == val_data$y))}, error = function(e){
      #         return(c(NA,NA))
      #       })
      #   return(res)
      # }), MARGIN = 1, which.max)
      
      MSVM7_fit_1 <- MSVM7_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      MSVM7_fit_2 <- MSVM7_pri_opt(train_data$X,train_data$y, C=1, lambda=0, intercept = F)
      c(mean(predict(MSVM7_fit_1, oracle_data$X, rule = "dagger") == oracle_data$y),
        mean(predict(MSVM7_fit_2, oracle_data$X, rule = "0_max") == oracle_data$y))
    }, error = function(e) {
      return(c(NA,NA))
    })
    
    MSVM8_acc <-  tryCatch({
      # opt_idx <- which.max(sapply(C_list, function(C){
      #   tryCatch({
      #     MSVM8_fit <- MSVM8_pri_opt(train_data$X,train_data$y,C, intercept = F)
      #     mean(predict(MSVM8_fit, val_data$X) == val_data$y)},
      #     error =function(e){
      #       return(NA)
      #     })
      # }))
      MSVM8_fit <- MSVM8_pri_opt(train_data$X, train_data$y, C=1, lambda=0, intercept = F)
      mean(predict(MSVM8_fit, oracle_data$X) == oracle_data$y)
    }, error = function(e) {
      return(NA)
    })
    
    res <-
      c(WW_acc = WW_acc,
        CS_acc = CS_acc,
        Duchi_acc = Duchi_acc,
        MDuchi_acc = MDuchi_acc,
        New1_acc = New1_acc,
        New3_acc = New3_acc,
        OVA_acc = OVA_acc,
        LLW_acc = LLW_acc,
        MSVM7_acc = MSVM7_acc,
        MSVM8_acc = MSVM8_acc)
    return(res)
  }, n = n_list[j])
  
  res.array[j,,] <- acc.res
}

stopCluster(cl)


summary.res.arry <- apply(res.array, c(1,2), function(x){mean(x,na.rm = T)})

dimnames(summary.res.arry)[[2]] <-  c("WW", "CS", "Duchi", "MDuchi", "New1,d","New1,0", "New3,d","New3,0", "OVA", "LLW", "MSVM7,d","MSVM7,0", "MSVM8")
dimnames(summary.res.arry)[[1]] <- n_list 
round((summary.res.arry[,])*100,2)
# 
# 
# save(rep_n, n_list, res.array, summary.res.arry, file = 
#        paste("~/Desktop/Multiclass Classification/MSVM Code/Noise Circle Ridgeless, rep=",rep_n,".Rdata",sep=""))

##100: OVA, LLW, MSVM8  68.83 67.72 68.11
##100 - c(68.83, 67.72, 68.11)
