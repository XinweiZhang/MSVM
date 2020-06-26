MSVM_Weak_Hard_opt <- function(X,y,type){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})

  X1 = X[y==1,,drop=F]
  X2 = X[y==2,,drop=F]
  X3 = X[y==3,,drop=F]
  p <- ncol(X)
  w1 <- Variable(p)
  w2 <- Variable(p)
  w3 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  b3 <- Variable(1)

  if(type == "max_of_all"){
    objective <- Minimize(max(sum_squares(w2 - w1), sum_squares(w3 - w1),
                              sum_squares(w3 - w2)))
  }else if(type == "max_of_sum"){
    objective <- Minimize(max(sum_squares(w2 - w1) + sum_squares(w3 - w1),
                              sum_squares(w1 - w2) + sum_squares(w3 - w2),
                              sum_squares(w1 - w3) + sum_squares(w2 - w3)))
  }else if(type == "sum_of_max"){
    objective <- Minimize(sum(max(sum_squares(w2 - w1), sum_squares(w3 - w1)),
                              max(sum_squares(w1 - w2), sum_squares(w3 - w2)),
                              max(sum_squares(w1 - w3), sum_squares(w2 - w3))))
  }else if(type == "Duchi"){
    objective <- Minimize(sum(sum_squares(w1), sum_squares(w2), sum_squares(w3)))
  }else if(type == "New1")
  {
    objective <- Minimize(sum(sum_squares(w1-w3), sum_squares(w2-w3)))
  }

  constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                      X1 %*% (w1 - w2) + b1 - b2 >= 1,
                      X1 %*% (w1 - w3) + b1 - b3 >= 1,
                      X2 %*% (w2 - w1) + b2 - b1 >= 1,
                      X2 %*% (w2 - w3) + b2 - b3 >= 1,
                      X3 %*% (w3 - w1) + b3 - b1 >= 1,
                      X3 %*% (w3 - w2) + b3 - b2 >= 1)

  WW <- Problem(objective, constraints)
  CVXR_WW <- solve(WW, solver = "MOSEK")
  CVXR_WW_w1 <- CVXR_WW$getValue(w1)
  CVXR_WW_b1 <- CVXR_WW$getValue(b1)
  CVXR_WW_w2 <- CVXR_WW$getValue(w2)
  CVXR_WW_b2 <- CVXR_WW$getValue(b2)
  CVXR_WW_w3 <- CVXR_WW$getValue(w3)
  CVXR_WW_b3 <- CVXR_WW$getValue(b3)

  beta1 <- c(CVXR_WW_w1,CVXR_WW_b1)
  beta2 <- c(CVXR_WW_w2,CVXR_WW_b2)
  beta3 <- c(CVXR_WW_w3,CVXR_WW_b3)
  CVXR_WW_primary_beta <- rbind(beta1,beta2,beta3)

  return(CVXR_WW_primary_beta)
}
# MSVM_Weak_Hard_opt <- function(X,y,type){
#   class_idx <- sort(unique(y))
#   Y <- sapply(class_idx, function(id){as.numeric(y==id)})
# 
#   X1 = X[y==1,,drop=F]
#   X2 = X[y==2,,drop=F]
#   X3 = X[y==3,,drop=F]
#   p <- ncol(X)
#   w1 <- Variable(p)
#   w2 <- Variable(p)
#   w3 <- Variable(p)
#   b1 <- Variable(1)
#   b2 <- Variable(1)
#   b3 <- Variable(1)
# 
#   if(type == "max_of_all"){
#     objective <- Minimize(max(sum_squares(w2 - w1), sum_squares(w3 - w1),
#                               sum_squares(w3 - w2)))
#   }else if(type == "max_of_sum"){
#     objective <- Minimize(max(sum_squares(w2 - w1) + sum_squares(w3 - w1),
#                               sum_squares(w1 - w2) + sum_squares(w3 - w2),
#                               sum_squares(w1 - w3) + sum_squares(w2 - w3)))
#   }else if(type == "sum_of_max"){
#     objective <- Minimize(sum(max(sum_squares(w2 - w1), sum_squares(w3 - w1)),
#                               max(sum_squares(w1 - w2), sum_squares(w3 - w2)),
#                               max(sum_squares(w1 - w3), sum_squares(w2 - w3))))
#   }else if(type == "Duchi"){
#     objective <- Minimize(sum(sum_squares(w1), sum_squares(w2), sum_squares(w3)))
#   }else if(type == "New1")
#   {
#     objective <- Minimize(sum(sum_squares(w1-w3), sum_squares(w2-w3)))
#   }
# 
#   constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
#                       X1 %*% (w1 - w2) + b1 - b2 >= 1*sum_squares(w1 - w2),
#                       X1 %*% (w1 - w3) + b1 - b3 >= 1*sum_squares(w1 - w3),
#                       X2 %*% (w2 - w1) + b2 - b1 >= 1*sum_squares(w2 - w1),
#                       X2 %*% (w2 - w3) + b2 - b3 >= 1*sum_squares(w2 - w3),
#                       X3 %*% (w3 - w1) + b3 - b1 >= 1*sum_squares(w3 - w1),
#                       X3 %*% (w3 - w2) + b3 - b2 >= 1*sum_squares(w3 - w2))
# 
#   WW <- Problem(objective, constraints)
#   CVXR_WW <- solve(WW, solver = "MOSEK")
#   CVXR_WW_w1 <- CVXR_WW$getValue(w1)
#   CVXR_WW_b1 <- CVXR_WW$getValue(b1)
#   CVXR_WW_w2 <- CVXR_WW$getValue(w2)
#   CVXR_WW_b2 <- CVXR_WW$getValue(b2)
#   CVXR_WW_w3 <- CVXR_WW$getValue(w3)
#   CVXR_WW_b3 <- CVXR_WW$getValue(b3)
# 
#   beta1 <- c(CVXR_WW_w1,CVXR_WW_b1)
#   beta2 <- c(CVXR_WW_w2,CVXR_WW_b2)
#   beta3 <- c(CVXR_WW_w3,CVXR_WW_b3)
#   CVXR_WW_primary_beta <- rbind(beta1,beta2,beta3)
# 
#   return(CVXR_WW_primary_beta)
# }
# 
# Rescale3_beta <- function(X,y,beta, type){
#   class_idx <- sort(unique(y))
#   Y <- sapply(class_idx, function(id){as.numeric(y==id)})
#   
#   X1 = X[y==1,,drop=F]
#   X2 = X[y==2,,drop=F]
#   X3 = X[y==3,,drop=F]
#   p <- ncol(X)
#   w1 <- beta[1,1:2]
#   w2 <- beta[2,1:2]
#   w3 <- beta[3,1:2]
#   b1 <- beta[1,3]
#   b2 <- beta[2,3]
#   b3 <- beta[3,3]
#   gamma <- Variable(1) 
#   
#   if(type == "max_of_all"){
#     objective <- Minimize(gamma^2*max(sum_squares(w2 - w1), sum_squares(w3 - w1), 
#                               sum_squares(w3 - w2))) 
#   }else if(type == "max_of_sum"){
#     objective <- Minimize(gamma^2*max(sum_squares(w2 - w1) + sum_squares(w3 - w1),
#                               sum_squares(w1 - w2) + sum_squares(w3 - w2),
#                               sum_squares(w1 - w3) + sum_squares(w2 - w3)))
#   }else if(type == "sum_of_max"){
#     objective <- Minimize(gamma^2*sum(max(sum_squares(w2 - w1), sum_squares(w3 - w1)),
#                               max(sum_squares(w1 - w2), sum_squares(w3 - w2)),
#                               max(sum_squares(w1 - w3), sum_squares(w2 - w3))))
#   }else if(type == "Duchi"){
#     objective <- Minimize(gamma^2*sum(sum_squares(w1), sum_squares(w2), sum_squares(w3)))
#   }else if(type == "New1")
#   {
#     objective <- Minimize(gamma^2*sum(sum_squares(w1-w2), sum_squares(w2-w3)))
#   }
#   
#   constraints <- list(gamma*(w1+w2+w3) == 0, gamma*(b1+b2+b3) ==0,
#                       gamma*(X1 %*% (w1 - w2) + b1 - b2) >= 1,
#                       gamma*(X1 %*% (w1 - w3) + b1 - b3) >= 1,
#                       gamma*(X2 %*% (w2 - w1) + b2 - b1) >= 1,
#                       gamma*(X2 %*% (w2 - w3) + b2 - b3) >= 1,
#                       gamma*(X3 %*% (w3 - w1) + b3 - b1) >= 1,
#                       gamma*(X3 %*% (w3 - w2) + b3 - b2) >= 1)
#   
#   WW <- Problem(objective, constraints)
#   CVXR_WW <- solve(WW, solver = "MOSEK") 
#   gamma <- CVXR_WW$getValue(gamma)
#   
#   return(list(gamma=gamma,value = CVXR_WW$value))
# }

# 
# MSVM_Weak_Hard_opt4 <- function(X,y,type){
#   class_idx <- sort(unique(y))
#   Y <- sapply(class_idx, function(id){as.numeric(y==id)})
# 
#   X1 = X[y==1,,drop=F]
#   X2 = X[y==2,,drop=F]
#   X3 = X[y==3,,drop=F]
#   X4 = X[y==4,,drop=F]
#   p <- ncol(X)
#   w1 <- Variable(p)
#   w2 <- Variable(p)
#   w3 <- Variable(p)
#   w4 <- Variable(p)
#   b1 <- Variable(1)
#   b2 <- Variable(1)
#   b3 <- Variable(1)
#   b4 <- Variable(1)
#   rho <- Variable(6)
# 
#   if(type == "max_of_all"){
#     objective <- Minimize(max(sum_squares(w2 - w1), sum_squares(w3 - w1), sum_squares(w4 - w1),
#                               sum_squares(w3 - w2), sum_squares(w4 - w2),
#                               sum_squares(w4 - w3)) )
#   }else if(type == "max_of_sum"){
#     objective <- Minimize(max(sum_squares(w2 - w1) + sum_squares(w3 - w1) + sum_squares(w4 - w1),
#                               (sum_squares(w1 - w2) + sum_squares(w3 - w2) + sum_squares(w4 - w2)),
#                               sum_squares(w1 - w3) + sum_squares(w2 - w3) + sum_squares(w4 - w3),
#                               sum_squares(w4 - w1) + sum_squares(w4 - w2) + sum_squares(w4 - w3)))
#   }else if(type == "sum_of_max"){
#     objective <- Minimize(sum(max(sum_squares(w2 - w1), sum_squares(w3 - w1), sum_squares(w4 - w1)),
#                               max(sum_squares(w1 - w2),sum_squares(w3 - w2), sum_squares(w4 - w2)),
#                               max(sum_squares(w1 - w3), sum_squares(w2 - w3), sum_squares(w4 - w3)),
#                               max(sum_squares(w4 - w1), sum_squares(w4 - w2), sum_squares(w4 - w3))))
#   }else if(type == "Duchi"){
#     objective <- Minimize(sum(sum_squares(w1), sum_squares(w2), sum_squares(w3), sum_squares(w4)) -  sum(rho))
#   }else if(type == "New1")
#   {
#     objective <- Minimize(sum(sum_squares(w1-w4), sum_squares(w2-w4), sum_squares(w3-w4)))
#   }
# 
#   constraints <- list(w1+w2+w3+w4== 0, b1+b2+b3+b4 ==0,
#                       X1 %*% (w1 - w2) + b1 - b2 >= rho[1],
#                       X1 %*% (w1 - w3) + b1 - b3 >= rho[2],
#                       X1 %*% (w1 - w4) + b1 - b4 >= rho[3],
#                       X2 %*% (w2 - w1) + b2 - b1 >= rho[1],
#                       X2 %*% (w2 - w3) + b2 - b3 >= rho[4],
#                       X2 %*% (w2 - w4) + b2 - b4 >= rho[5],
#                       X3 %*% (w3 - w1) + b3 - b1 >= rho[2],
#                       X3 %*% (w3 - w2) + b3 - b2 >= rho[4],
#                       X3 %*% (w3 - w4) + b3 - b4 >= rho[6],
#                       X4 %*% (w4 - w1) + b4 - b1 >= rho[3],
#                       X4 %*% (w4 - w2) + b4 - b2 >= rho[5],
#                       X4 %*% (w4 - w3) + b4 - b3 >= rho[6],
#                       rho >= 0
#                       )
# 
#   WW <- Problem(objective, constraints)
#   CVXR_WW <- solve(WW, solver = "MOSEK")
#   CVXR_WW_w1 <- CVXR_WW$getValue(w1)
#   CVXR_WW_b1 <- CVXR_WW$getValue(b1)
#   CVXR_WW_w2 <- CVXR_WW$getValue(w2)
#   CVXR_WW_b2 <- CVXR_WW$getValue(b2)
#   CVXR_WW_w3 <- CVXR_WW$getValue(w3)
#   CVXR_WW_b3 <- CVXR_WW$getValue(b3)
#   CVXR_WW_w4 <- CVXR_WW$getValue(w4)
#   CVXR_WW_b4 <- CVXR_WW$getValue(b4)
# 
#   beta1 <- c(CVXR_WW_w1,CVXR_WW_b1)
#   beta2 <- c(CVXR_WW_w2,CVXR_WW_b2)
#   beta3 <- c(CVXR_WW_w3,CVXR_WW_b3)
#   beta4 <- c(CVXR_WW_w4,CVXR_WW_b4)
# 
#   CVXR_WW_primary_beta <- rbind(beta1,beta2,beta3,beta4)
# 
#   return(CVXR_WW_primary_beta)
# }

MSVM_Weak_Hard_opt4 <- function(X,y,type){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})

  X1 = X[y==1,,drop=F]
  X2 = X[y==2,,drop=F]
  X3 = X[y==3,,drop=F]
  X4 = X[y==4,,drop=F]
  p <- ncol(X)
  w1 <- Variable(p)
  w2 <- Variable(p)
  w3 <- Variable(p)
  w4 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  b3 <- Variable(1)
  b4 <- Variable(1)


  if(type == "max_of_all"){
    objective <- Minimize(max(sum_squares(w2 - w1), sum_squares(w3 - w1), sum_squares(w4 - w1),
                              sum_squares(w3 - w2), sum_squares(w4 - w2),
                              sum_squares(w4 - w3)) )
  }else if(type == "max_of_sum"){
    objective <- Minimize(max(sum_squares(w2 - w1) + sum_squares(w3 - w1) + sum_squares(w4 - w1),
                              (sum_squares(w1 - w2) + sum_squares(w3 - w2) + sum_squares(w4 - w2)),
                              sum_squares(w1 - w3) + sum_squares(w2 - w3) + sum_squares(w4 - w3),
                              sum_squares(w4 - w1) + sum_squares(w4 - w2) + sum_squares(w4 - w3)))
  }else if(type == "sum_of_max"){
    objective <- Minimize(sum(max(sum_squares(w2 - w1), sum_squares(w3 - w1), sum_squares(w4 - w1)),
                              max(sum_squares(w1 - w2),sum_squares(w3 - w2), sum_squares(w4 - w2)),
                              max(sum_squares(w1 - w3), sum_squares(w2 - w3), sum_squares(w4 - w3)),
                              max(sum_squares(w4 - w1), sum_squares(w4 - w2), sum_squares(w4 - w3))))
  }else if(type == "Duchi"){
    objective <- Minimize(sum(sum_squares(w1), sum_squares(w2), sum_squares(w3), sum_squares(w4)))
  }else if(type == "New1")
  {
    objective <- Minimize(sum(sum_squares(w1-w4), sum_squares(w2-w4), sum_squares(w3-w4)))
  }

  constraints <- list(w1+w2+w3+w4== 0, b1+b2+b3+b4 ==0,
                      X1 %*% (w1 - w2) + b1 - b2 >= 1,
                      X1 %*% (w1 - w3) + b1 - b3 >= 1,
                      X1 %*% (w1 - w4) + b1 - b4 >= 1,
                      X2 %*% (w2 - w1) + b2 - b1 >= 1,
                      X2 %*% (w2 - w3) + b2 - b3 >= 1,
                      X2 %*% (w2 - w4) + b2 - b4 >= 1,
                      X3 %*% (w3 - w1) + b3 - b1 >= 1,
                      X3 %*% (w3 - w2) + b3 - b2 >= 1,
                      X3 %*% (w3 - w4) + b3 - b4 >= 1,
                      X4 %*% (w4 - w1) + b4 - b1 >= 1,
                      X4 %*% (w4 - w2) + b4 - b2 >= 1,
                      X4 %*% (w4 - w3) + b4 - b3 >= 1
                      )

  WW <- Problem(objective, constraints)
  CVXR_WW <- solve(WW, solver = "MOSEK")
  CVXR_WW_w1 <- CVXR_WW$getValue(w1)
  CVXR_WW_b1 <- CVXR_WW$getValue(b1)
  CVXR_WW_w2 <- CVXR_WW$getValue(w2)
  CVXR_WW_b2 <- CVXR_WW$getValue(b2)
  CVXR_WW_w3 <- CVXR_WW$getValue(w3)
  CVXR_WW_b3 <- CVXR_WW$getValue(b3)
  CVXR_WW_w4 <- CVXR_WW$getValue(w4)
  CVXR_WW_b4 <- CVXR_WW$getValue(b4)

  beta1 <- c(CVXR_WW_w1,CVXR_WW_b1)
  beta2 <- c(CVXR_WW_w2,CVXR_WW_b2)
  beta3 <- c(CVXR_WW_w3,CVXR_WW_b3)
  beta4 <- c(CVXR_WW_w4,CVXR_WW_b4)

  CVXR_WW_primary_beta <- rbind(beta1,beta2,beta3,beta4)

  return(CVXR_WW_primary_beta)
}




Rescale4_beta <- function(X,y,beta, type){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  
  X1 = X[y==1,,drop=F]
  X2 = X[y==2,,drop=F]
  X3 = X[y==3,,drop=F]
  X4 = X[y==4,,drop=F]
  p <- ncol(X)
  w1 <- beta[1,1:2]
  w2 <- beta[2,1:2]
  w3 <- beta[3,1:2]
  w4 <- beta[4,1:2]
  b1 <- beta[1,3]
  b2 <- beta[2,3]
  b3 <- beta[3,3]
  b4 <- beta[4,3]
  gamma <- Variable(1) 
  
  if(type == "max_of_all"){
    objective <- Minimize(gamma^2*max(sum_squares(w2 - w1), sum_squares(w3 - w1), sum_squares(w4 - w1), 
                                    sum_squares(w3 - w2), sum_squares(w4 - w2),
                                    sum_squares(w4 - w3)))
  }else if(type == "max_of_sum"){
    objective <- Minimize(gamma^2*max(sum_squares(w2 - w1) + sum_squares(w3 - w1) + sum_squares(w4 - w1),
                                    sum_squares(w1 - w2) + sum_squares(w3 - w2) + sum_squares(w4 - w2),
                                    sum_squares(w1 - w3) + sum_squares(w2 - w3) + sum_squares(w4 - w3), 
                                    sum_squares(w4 - w1) + sum_squares(w4 - w2) + sum_squares(w4 - w3)))
  }else if(type == "sum_of_max"){
    objective <- Minimize(gamma^2*sum(max(sum_squares(w2 - w1), sum_squares(w3 - w1), sum_squares(w4 - w1)), 
                                    max(sum_squares(w1 - w2), sum_squares(w3 - w2), sum_squares(w4 - w2)),
                                    max(sum_squares(w1 - w3), sum_squares(w2 - w3), sum_squares(w4 - w3)),
                                    max(sum_squares(w4 - w1), sum_squares(w4 - w2), sum_squares(w4 - w3))))
  }else if(type == "Duchi"){
    objective <- Minimize(gamma^2*sum(sum_squares(w1), sum_squares(w2), sum_squares(w3) + sum_squares(w4)))
  }else if(type == "New1"){
    objective <- Minimize(gamma^2*sum(sum_squares(w1-w4), sum_squares(w2-w4), sum_squares(w3-w4)))
  }
  
  constraints <- list(gamma*(w1+w2+w3+w4)== 0, gamma*(b1+b2+b3+b4) ==0,
                      gamma*(X1 %*% (w1 - w2) + b1 - b2) >= 1,
                      gamma*(X1 %*% (w1 - w3) + b1 - b3) >= 1,
                      gamma*(X1 %*% (w1 - w4) + b1 - b4) >= 1,
                      gamma*(X2 %*% (w2 - w1) + b2 - b1) >= 1,
                      gamma*(X2 %*% (w2 - w3) + b2 - b3) >= 1,
                      gamma*(X2 %*% (w2 - w4) + b2 - b4) >= 1,
                      gamma*(X3 %*% (w3 - w1) + b3 - b1) >= 1,
                      gamma*(X3 %*% (w3 - w2) + b3 - b2) >= 1,
                      gamma*(X3 %*% (w3 - w4) + b3 - b4) >= 1,
                      gamma*(X4 %*% (w4 - w1) + b4 - b1) >= 1,
                      gamma*(X4 %*% (w4 - w2) + b4 - b2) >= 1,
                      gamma*(X4 %*% (w4 - w3) + b4 - b3) >= 1
  )
  
  WW <- Problem(objective, constraints)
  CVXR_WW <- solve(WW, solver = "MOSEK") 
  gamma <- CVXR_WW$getValue(gamma)
  
  return(list(gamma=gamma,value = CVXR_WW$value))
}

plot_decision_boundary <- function(X, y, beta, title = NULL, np_resolution = 500, xlim =NULL, ylim = NULL, dagger_rule_w = F, dagger_rule_s = F){
  X_dat <- as.data.frame(X)
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$V1,X_dat$V2))-2, ceiling(max(X_dat$V1,X_dat$V2))+2)
  }
  if(is.null(ylim)){
    ylim = c( floor(min(X_dat$V2,X_dat$V1))-2, ceiling(max(X_dat$V2,X_dat$V1))+2)
  }
  
  nd.x = seq(from = xlim[1], to =  xlim[2], length.out = np_resolution)
  nd.y = seq(from = ylim[1], to =  ylim[2], length.out = np_resolution)
  nd = expand.grid(Var1 = nd.x, Var2 = nd.y)
  if(nrow(beta) == length(unique(y))){
    prd = apply(as.matrix(cbind(nd,1))%*%t(beta), MARGIN = 1, FUN = which.max)
    prd_p = as.matrix(cbind(nd,1))%*%t(beta)
    idx = apply(apply(prd_p,MARGIN = 1, FUN = max)%*%matrix(1,ncol=nrow(beta)) - prd_p, MARGIN = 1, FUN = function(x){ sum(x>=1) == (nrow(beta)-1)})
    band = prd
    band[!idx]=0
  }else{
    if(dagger_rule_w==T && dagger_rule_s == F){
      prd_p = as.matrix(cbind(nd,1))%*%t(beta)
      prd_p =  cbind(prd_p, 0 - apply(as.matrix(cbind(prd_p[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(prd_p[,2]+1,0)), MARGIN = 1, FUN = max))
      prd = apply(prd_p, MARGIN = 1, FUN = which.max)
      idx = apply(apply(prd_p,MARGIN = 1, FUN = max)%*%matrix(1,ncol=nrow(beta)+1) - prd_p, MARGIN = 1, FUN = function(x){ sum(x>=1) == (nrow(beta))})
      band = prd
      band[!idx]=0
    }else if(dagger_rule_w==F && dagger_rule_s == T){
      prd_p = as.matrix(cbind(nd,1))%*%t(beta)
      prd_p =  cbind(prd_p, 1 - apply(as.matrix(cbind(prd_p[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(prd_p[,2]+1,0)), MARGIN = 1, FUN = max))
      prd = apply(prd_p, MARGIN = 1, FUN = which.max)
    
    }else{
      prd_p = as.matrix(cbind(nd,1))%*%t(rbind(beta,0))
      prd = apply(prd_p, MARGIN = 1, FUN = which.max)
      idx = apply(apply(prd_p,MARGIN = 1, FUN = max)%*%matrix(1,ncol=nrow(beta)) - prd_p, MARGIN = 1, FUN = function(x){ sum(x>=1) == (nrow(beta))})
      band = prd
      band[!idx]=0
    }
  }
  
  # op <-  par(mfrow = c(1,1), mar=c(5.1, 4.1, 4.1, 7), xpd=TRUE)
  plot(X_dat$V1, X_dat$V2, col = as.factor(y), ylim=ylim, xlim=xlim, xlab ="X", ylab = "Y", main = title)
  
  contour(x = nd.x, y = nd.y, z = matrix(prd, nrow = np_resolution, ncol = np_resolution), 
          levels = unique(y), add = TRUE, drawlabels = FALSE)
  
  contour(x = nd.x, y = nd.y, z = matrix(band, nrow = np_resolution, ncol = np_resolution), 
          levels = unique(y), add = TRUE, drawlabels = FALSE, lty = 2)
  
  
  # legend("topright", inset=c(-0.3,0),legend = sapply(unique(y), function(x){paste("Class ",x)}), col= unique(y), pch = 1)
  
  # legend("bottomright", legend = sapply(unique(y), function(x){paste("Class ",x)}), col= unique(y), pch = 1)
  # par(op)
}



plot_nearest_neighbour_decision_boundary <- function(X, y, dis = "L2", title = NULL, np_resolution = 500, xlim =NULL, ylim = NULL, dagger_rule = T){
  X_dat <- as.data.frame(X)
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$V1,X_dat$V2))-2, ceiling(max(X_dat$V1,X_dat$V2))+2)
  }
  if(is.null(ylim)){
    ylim = c( floor(min(X_dat$V2,X_dat$V1))-2, ceiling(max(X_dat$V2,X_dat$V1))+2)
  }
  
  
  nd.x = seq(from = xlim[1], to =  xlim[2], length.out = np_resolution)
  nd.y = seq(from = ylim[1], to =  ylim[2], length.out = np_resolution)
  nd = expand.grid(Var1 = nd.x, Var2 = nd.y)
  
  if(dis == "L2"){
    prd <- apply(nd, MARGIN = 1, function(s,X){which.min(apply(X,MARGIN=1,FUN = function(t){sqrt(sum((s-t)^2))}))},X)
  }else{
    prd <- apply(nd, MARGIN = 1, function(s,X){which.min(apply(X,MARGIN=1,FUN = function(t){sum(abs(s-t))}))},X)
  }
  # op <-  par(mfrow = c(1,1), mar=c(5.1, 4.1, 4.1, 7), xpd=TRUE)
  plot(X_dat$V1, X_dat$V2, col = as.factor(y), ylim=ylim, xlim=xlim, xlab ="X", ylab = "Y", main = title)
  
  contour(x = nd.x, y = nd.y, z = matrix(prd, nrow = np_resolution, ncol = np_resolution), 
          levels = unique(y), add = TRUE, drawlabels = FALSE)
  
  
  legend("bottomright", legend = sapply(unique(y), function(x){paste("Class ",x)}), col= unique(y), pch = 1)
  # par(op)
}
