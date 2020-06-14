WW_pri_opt <- function(X,y,C){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  objective <- Minimize(sum_squares(w)/2 +  C*sum_entries((1-Y)*slack))
  
  constraints <- list( sum_entries(b) == 0,
                       sum_entries(w, axis = 2) ==0,
                       ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >=1-slack,
                       slack >=0)
  
  WW <- Problem(objective, constraints)
  CVXR_WW <- solve(WW, solver = "MOSEK")
  
  return(cbind(CVXR_WW$getValue(w),CVXR_WW$getValue(b)))
}


CS_pri_opt <- function(X,y,C){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(max_entries((1-Y)*slack, axis = 1)))
  
  constraints <- list( sum_entries(b) == 0,
                       sum_entries(w, axis = 2) ==0,
                       ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-slack),
                       slack >=0)
  
  CS <- Problem(objective, constraints)
  CVXR_CS <- solve(CS, solver = "MOSEK")
  
  return(cbind(CVXR_CS$getValue(w),CVXR_CS$getValue(b)))
}



Duchi_pri_opt <- function(X,y,C){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  objective <- Minimize(sum_squares(w)/2 +  C*sum_entries((1-Y)*slack))
  
  
  slack1 <- Variable(rows = nrow(X1), cols = 3)
  slack2 <- Variable(rows = nrow(X2), cols = 3)
  slack3 <- Variable(rows = nrow(X3), cols = 3)
  xi1 <- Variable(rows = nrow(X1))
  xi2 <- Variable(rows = nrow(X2))
  xi3 <- Variable(rows = nrow(X3))

  objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum(xi1)+sum(xi2)+sum(xi3)))
  constraints <- list(b1+b2+b3 ==0,
                      slack1[,1] == 1,
                      X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
                      X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,3],
                      slack2[,2] == 1,
                      X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                      X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,3],
                      slack3[,3] == 1,
                      X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
                      X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
                      xi1 >= max_entries(slack1, axis = 1)-1,
                      xi2 >= max_entries(slack2, axis = 1)-1,
                      xi3 >= max_entries(slack3, axis = 1)-1,
                      xi1 >= (sum_entries(slack1, axis = 1) - min_entries(slack1, axis = 1))/2-1/2,
                      xi2 >= (sum_entries(slack2, axis = 1) - min_entries(slack2, axis = 1))/2-1/2,
                      xi3 >= (sum_entries(slack3, axis = 1) - min_entries(slack3, axis = 1))/2-1/2,
                      xi1 >= sum_entries(slack1, axis = 1)/3 - 1/3,
                      xi2 >= sum_entries(slack2, axis = 1)/3 - 1/3,
                      xi3 >= sum_entries(slack3, axis = 1)/3 - 1/3,
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)
  
  
  Duchi <- Problem(objective, constraints)
  CVXR_Duchi <- solve(Duchi, solver = "MOSEK")
  
  
  
  CVXR_Duchi_b1 <- CVXR_Duchi$getValue(b1)
  CVXR_Duchi_b2 <- CVXR_Duchi$getValue(b2)
  CVXR_Duchi_b3 <- CVXR_Duchi$getValue(b3)


  beta1 <- c(w1,CVXR_Duchi_b1)
  beta2 <- c(w2,CVXR_Duchi_b2)
  beta3 <- c(w3,CVXR_Duchi_b3)

  CVXR_Duchi_primary_beta <- rbind(beta1,beta2,beta3)
  return(CVXR_Duchi_primary_beta)
}

MDuchi_pri_opt <- function(X,y,C){
  class_idx <- unique(y)
  X1 <- X[y==class_idx[1],]
  X2 <- X[y==class_idx[2],]
  X3 <- X[y==class_idx[3],]
  p <- ncol(X)
  
  w1 <- Variable(p)
  w2 <- Variable(p)
  w3 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  b3 <- Variable(1)
  
  
  slack1 <- Variable(rows = nrow(X1), cols = 2)
  slack2 <- Variable(rows = nrow(X2), cols = 2)
  slack3 <- Variable(rows = nrow(X3), cols = 2)
  
  xi1 <- Variable(rows = nrow(X1))
  xi2 <- Variable(rows = nrow(X2))
  xi3 <- Variable(rows = nrow(X3))
  
  objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum(xi1)+sum(xi2)+sum(xi3)))
  constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                      X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,1],
                      X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,2],
                      X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                      X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,2],
                      X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
                      X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
                      xi1 >= max_entries(slack1, axis = 1)/2,
                      xi2 >= max_entries(slack2, axis = 1)/2,
                      xi3 >= max_entries(slack3, axis = 1)/2,
                      xi1 >= sum_entries(slack1, axis = 1)/3,
                      xi2 >= sum_entries(slack2, axis = 1)/3,
                      xi3 >= sum_entries(slack3, axis = 1)/3,
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)
  
  M_Duchi <- Problem(objective, constraints)
  M_Duchi_primary_time <- system.time(CVXR_M_Duchi <- solve(M_Duchi, solver = "MOSEK"))
  
  CVXR_M_Duchi_w1 <- CVXR_M_Duchi$getValue(w1)
  CVXR_M_Duchi_b1 <- CVXR_M_Duchi$getValue(b1)
  CVXR_M_Duchi_w2 <- CVXR_M_Duchi$getValue(w2)
  CVXR_M_Duchi_b2 <- CVXR_M_Duchi$getValue(b2)
  CVXR_M_Duchi_w3 <- CVXR_M_Duchi$getValue(w3)
  CVXR_M_Duchi_b3 <- CVXR_M_Duchi$getValue(b3)
  
  beta1 <- c(CVXR_M_Duchi_w1,CVXR_M_Duchi_b1)
  beta2 <- c(CVXR_M_Duchi_w2,CVXR_M_Duchi_b2)
  beta3 <- c(CVXR_M_Duchi_w3,CVXR_M_Duchi_b3)
  
  CVXR_M_Duchi_primary_beta <- rbind(beta1,beta2,beta3)
  return(CVXR_M_Duchi_primary_beta)
}


OVA_pri_opt <- function(X,y,C){
  class_idx <- unique(y)
  X1 <- X[y==class_idx[1],]
  X2 <- X[y==class_idx[2],]
  X3 <- X[y==class_idx[3],]
  p <- ncol(X)
  
  
  w1 <- Variable(p)
  w2 <- Variable(p)
  w3 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  b3 <- Variable(1)
  
  
  slack1 <- Variable(rows = nrow(X1), cols = 3)
  slack2 <- Variable(rows = nrow(X2), cols = 3)
  slack3 <- Variable(rows = nrow(X3), cols = 3)

  objective <- Minimize(1/2*sum_squares(vstack(w1,w2,w3)) +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))
  constraints <- list(X1 %*% w1  + b1  >= 1 - slack1[,1],
                      X1 %*% w2  + b2 <= -1 + slack1[,2],
                      X1 %*% w3  + b3 <= -1 + slack1[,3],
                      X2 %*% w1  + b1 <= -1 + slack2[,1],
                      X2 %*% w2  + b2  >= 1 - slack2[,2],
                      X2 %*% w3  + b3 <= -1 + slack2[,3],
                      X3 %*% w1  + b1 <= -1 + slack3[,1],
                      X3 %*% w2  + b2 <= -1 +  slack3[,2],
                      X3 %*% w3  + b3  >= 1- slack3[,3], 
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)
  
  OVA <- Problem(objective, constraints)
  CVXR_OVA <- solve(OVA, solver="MOSEK")
  
  CVXR_OVA_w1 <- CVXR_OVA$getValue(w1)
  CVXR_OVA_b1 <- CVXR_OVA$getValue(b1)
  CVXR_OVA_w2 <- CVXR_OVA$getValue(w2)
  CVXR_OVA_b2 <- CVXR_OVA$getValue(b2)
  CVXR_OVA_w3 <- CVXR_OVA$getValue(w3)
  CVXR_OVA_b3 <- CVXR_OVA$getValue(b3)
  
  OVA_beta1 <- c(CVXR_OVA_w1,CVXR_OVA_b1)
  OVA_beta2 <- c(CVXR_OVA_w2,CVXR_OVA_b2)
  OVA_beta3 <- c(CVXR_OVA_w3,CVXR_OVA_b3)
  

  CVXR_OVA_primary_beta <- rbind(OVA_beta1,OVA_beta2,OVA_beta3)
  return(CVXR_OVA_primary_beta)
}


MSVM8_pri_opt <- function(X,y,C){
  class_idx <- unique(y)
  X1 <- X[y==class_idx[1],]
  X2 <- X[y==class_idx[2],]
  X3 <- X[y==class_idx[3],]
  p <- ncol(X)
  
  w1 <- Variable(p)
  w2 <- Variable(p)
  w3 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  b3 <- Variable(1)
  
  slack1 <- Variable(rows = nrow(X1), cols = 3)
  slack2 <- Variable(rows = nrow(X2), cols = 3)
  slack3 <- Variable(rows = nrow(X3), cols = 3) 

  objective <- Minimize(1/2*sum_squares(vstack(w1,w2,w3)) +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))
  constraints <- list(X1 %*% w1  + b1  >= 1 - sum_entries(slack1, axis=1),
                      X1 %*% w2  + b2 <= -1 + slack1[,2],
                      X1 %*% w3  + b3 <= -1 + slack1[,3],
                      X2 %*% w1  + b1 <= -1 + slack2[,1],
                      X2 %*% w2  + b2  >= 1-sum_entries(slack2, axis=1),
                      X2 %*% w3  + b3 <= -1 + slack2[,3],
                      X3 %*% w1  + b1 <= -1 + slack3[,1],
                      X3 %*% w2  + b2 <= -1 +  slack3[,2],
                      X3 %*% w3  + b3  >= 1-sum_entries(slack3, axis=1), 
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)
  
  MSVM8 <- Problem(objective, constraints)
  
  CVXR_MSVM8 <- solve(MSVM8, solver = "MOSEK")
  
  CVXR_MSVM8_w1 <- CVXR_MSVM8$getValue(w1)
  CVXR_MSVM8_b1 <- CVXR_MSVM8$getValue(b1)
  CVXR_MSVM8_w2 <- CVXR_MSVM8$getValue(w2)
  CVXR_MSVM8_b2 <- CVXR_MSVM8$getValue(b2)
  CVXR_MSVM8_w3 <- CVXR_MSVM8$getValue(w3)
  CVXR_MSVM8_b3 <- CVXR_MSVM8$getValue(b3)
  
  MSVM8_primary_beta <- rbind(c(CVXR_MSVM8_w1,CVXR_MSVM8_b1),
                              c(CVXR_MSVM8_w2,CVXR_MSVM8_b2),
                              c(CVXR_MSVM8_w3,CVXR_MSVM8_b3))
  return(MSVM8_primary_beta)
}




LLW_pri_opt <- function(X,y,C){
  class_idx <- unique(y)
  X1 <- X[y==class_idx[1],]
  X2 <- X[y==class_idx[2],]
  X3 <- X[y==class_idx[3],]
  p <- ncol(X)
  
  w1 <- Variable(p)
  w2 <- Variable(p)
  w3 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  b3 <- Variable(1)
  
  slack1 <- Variable(rows = nrow(X1), cols = 2)
  slack2 <- Variable(rows = nrow(X2), cols = 2)
  slack3 <- Variable(rows = nrow(X3), cols = 2)
  
  objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))
  constraints <- list(w1+w2+w3 == 0, b1+b2+b3 ==0,
                      X1 %*% w2  + b2 <= -1 + slack1[,1],
                      X1 %*% w3  + b3 <= -1 + slack1[,2],
                      X2 %*% w1  + b1 <= -1 + slack2[,1],
                      X2 %*% w3  + b3 <= -1 + slack2[,2],
                      X3 %*% w1  + b1 <= -1 + slack3[,1],
                      X3 %*% w2  + b2 <= -1 +  slack3[,2],
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)
  
  LLW <- Problem(objective, constraints)
  
  
  CVXR_LLW <- solve(LLW, solver = "MOSEK")
  
  CVXR_LLW_w1 <- CVXR_LLW$getValue(w1)
  CVXR_LLW_b1 <- CVXR_LLW$getValue(b1)
  CVXR_LLW_w2 <- CVXR_LLW$getValue(w2)
  CVXR_LLW_b2 <- CVXR_LLW$getValue(b2)
  CVXR_LLW_w3 <- CVXR_LLW$getValue(w3)
  CVXR_LLW_b3 <- CVXR_LLW$getValue(b3)
  
  LLW_beta1 <- c(CVXR_LLW_w1,CVXR_LLW_b1)
  LLW_beta2 <- c(CVXR_LLW_w2,CVXR_LLW_b2)
  LLW_beta3 <- c(CVXR_LLW_w3,CVXR_LLW_b3)
  
  LLW_primary_beta <- rbind(LLW_beta1, LLW_beta2, LLW_beta3)
  return(LLW_primary_beta)
}


MSVM7_pri_opt <- function(X,y,C){
  class_idx <- unique(y)
  X1 <- X[y==class_idx[1],]
  X2 <- X[y==class_idx[2],]
  X3 <- X[y==class_idx[3],]
  p <- ncol(X)

  w1 <- Variable(p)
  w2 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  
  slack1 <- Variable(rows = nrow(X1), cols = 2)
  slack2 <- Variable(rows = nrow(X2), cols = 2)
  slack3 <- Variable(rows = nrow(X3), cols = 2)
  
  objective <- Minimize(sum_squares(vstack(w1,w2))/2 +  C*(sum_entries(slack1) + sum_entries(slack2) + sum_entries(slack3)))
  constraints <- list(X1 %*% w1  + b1  >= 2 - sum_entries(slack1, axis=1),
                      X1 %*% w2  + b2 <= -1 + slack1[,2],
                      X2 %*% w1  + b1 <= -1 + slack2[,1],
                      X2 %*% w2  + b2  >= 2 - sum_entries(slack2, axis=1),
                      X3 %*% w1  + b1 <= -1 + slack3[,1],
                      X3 %*% w2  + b2 <= -1 +  slack3[,2],
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)
  
  MSVM7 <- Problem(objective, constraints)
  
  CVXR_MSVM7 <- solve(MSVM7, solver = "MOSEK")
  
  CVXR_MSVM7_w1 <- CVXR_MSVM7$getValue(w1)
  CVXR_MSVM7_b1 <- CVXR_MSVM7$getValue(b1)
  CVXR_MSVM7_w2 <- CVXR_MSVM7$getValue(w2)
  CVXR_MSVM7_b2 <- CVXR_MSVM7$getValue(b2)
  
  MSVM7_beta1 <- c(CVXR_MSVM7_w1,CVXR_MSVM7_b1)
  MSVM7_beta2 <- c(CVXR_MSVM7_w2,CVXR_MSVM7_b2)
  MSVM7_primary_beta <- rbind(MSVM7_beta1, MSVM7_beta2)
  return(MSVM7_primary_beta)
}


pred <- function(X_test,w, rule = "simple_max"){
  lpred <- cbind(X_test,1)%*%t(w)
  if(rule == "simple_max"){
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else
  {
    lpred <- cbind(lpred, 1 - apply(as.matrix(cbind(lpred[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(lpred[,2]+1,0)), MARGIN = 1, FUN = max))
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }
  return(y_pred)
}

data_generate <- function(n, sep =  1,  v = 1.5^2){
  # X1 <- mvrnorm(n, sep*c(0,2), diag(v,nrow=2))
  # X2 <- mvrnorm(n, sep*c(sqrt(3),-1), diag(v,nrow=2))
  # X3 <- mvrnorm(n, sep*c(-sqrt(3),-1), diag(v,nrow=2))
  ########  Plan II     #############
  X1 <- mvrnorm(n, sep*c(0,2), diag(c(8,1),nrow=2))
  X2 <- mvrnorm(n, sep*c(sqrt(3),-1),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2))
  X3 <- mvrnorm(n, sep*c(-sqrt(3),-1),  matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2))

  ########  Plan III     #############
  # X1 <- mvrnorm(n, sep*c(0,2), diag(c(1,8),nrow=2))
  # X2 <- mvrnorm(n, sep*c(sqrt(3),-1),   matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2))
  # X3 <- mvrnorm(n, sep*c(-sqrt(3),-1), matrix(c(4.5,3.5,3.5,4.5), nrow= 2))
  # 
  ########  Plan IV     #############
  # X1 <- mvrnorm(n, sep*c(0,2), diag(c(10,2),nrow=2))
  # X2 <- mvrnorm(n, sep*c(sqrt(3),-1),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2))
  # a <- matrix(c(cos(pi/3),-sin(pi/3), sin(pi/3), cos(pi/3)), nrow = 2)
  # Sigma3 <- a%*%diag(c(8,1))%*%t(a)
  # X3 <- mvrnorm(n, sep*c(-sqrt(3),-1), Sigma3)

  X <- rbind(X1,X2,X3)
  y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
  return(list(X = X,y = y))
}



Duchi_dual_opt <- function(X,y,C){
  n <- nrow(X)
  m <- length(unique(y))
  p <- ncol(X)
  Y <- sapply(unique(y), function(id){as.numeric(y==id)})
  
  alpha <- Variable(rows = n, cols = m)
  beta <- Variable(rows = n, cols = m)
  tau <- Variable(rows = n*m, cols = m)
  
  objective <- Maximize(-sum_squares(t(sum_entries(alpha, axis=1)%*%matrix(1,nrow=1,ncol = m)*Y - alpha)%*%X)/2 + 
                          sum_entries(alpha) - sum_entries(beta*sapply(1:m, function(j){rep(1/j,n)})))
  
  constraints <- list(sum_entries(Y*(sum_entries(alpha, axis = 1)%*%matrix(1,nrow=1,ncol = m)), axis = 2) == sum_entries(alpha, axis = 2),
                      C - sum_entries(beta, axis =1) ==0, 
                      vec(beta*sapply(1:m, function(j){rep(1/j,n)}))%*%matrix(1,nrow=1,ncol=m) - tau >=0, 
                      beta - reshape_expr(sum_entries(tau, axis =1), c(n,m))==0,
                      alpha>=0,
                      beta>=0,
                      tau>=0) 
  
  constraints <- c(constraints,lapply(0:(m-1), function(id){
    sum_entries(reshape_expr(tau, c(n, m*m))[,(id*m+1):((id+1)*m)], axis=1) - alpha[,id+1] >=0}))
  
  
  Duchi_dual <- Problem(objective, constraints)
  Duchi_dual_time <- system.time(Duchi_dual <- solve(Duchi_dual, solver = "MOSEK"))
  
  alpha_sol <- Duchi_dual$getValue(alpha)
  tau_sol <- Duchi_dual$getValue(tau)
  
  solve_w_and_b <- function(X, Y, alpha_sol, tau_sol, zero_tol = 4){
    m = ncol(Y)
    w <- t(rowSums(alpha_sol)%*%matrix(1,nrow=1,ncol = m)*Y - alpha_sol)%*%X
    lpred_dif <- X%*%t(w) - (Y*X%*%t(w))%*%matrix(1,nrow = m, ncol=m)
    tau_sum <- do.call(cbind,lapply(0:(m-1), function(id,tau){
      rowSums(tau[,(id*m+1):((id+1)*m)])}, tau = matrix(c(tau_sol), n, m*m, byrow = F))) 
    b_system_mat <- rep(1,m)
    b_system_v <- 0
    for(i in 1:m){
      for( j in (i+1):m)
      {
        if(i==m){
          break
        }
        lpred_dif_b <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(tau_sum - alpha_sol,zero_tol) >0) * (y==i)%*%matrix(1,nrow=1, ncol = m)
        lpred_dif_b[lpred_dif_b == 0] <- NA
        b_tmp <- colMeans(lpred_dif_b, na.rm = T)[j]
        
        lpred_dif_b_inv <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(tau_sum - alpha_sol,zero_tol) >0) * (y==j)%*%matrix(1,nrow=1, ncol = m)
        lpred_dif_b_inv[lpred_dif_b_inv == 0] <- NA
        b_tmp <- c(b_tmp ,-colMeans(lpred_dif_b_inv, na.rm = T)[i])
        
        if(is.na(mean(b_tmp,na.rm = T))!=T){
          tmp <- rep(0,m)
          tmp[i] <- 1
          tmp[j] <- -1
          b_system_mat <- rbind(b_system_mat, tmp)
          b_system_v <- c(b_system_v, mean(b_tmp,na.rm = T))
        }
      }
    }
    if(sum(svd(b_system_mat)$d>0)==m){
      b <- ginv(b_system_mat)%*%b_system_v
    }else{
      class_idx <- unique(y)
      X1 <- X[y==class_idx[1],]
      X2 <- X[y==class_idx[2],]
      X3 <- X[y==class_idx[3],]
      
      w1 <- w[1,]
      w2 <- w[2,]
      w3 <- w[3,]
      b1 <- Variable(1)
      b2 <- Variable(1)
      b3 <- Variable(1)
      
      slack1 <- Variable(rows = nrow(X1), cols = 3)
      slack2 <- Variable(rows = nrow(X2), cols = 3)
      slack3 <- Variable(rows = nrow(X3), cols = 3)
      xi1 <- Variable(rows = nrow(X1))
      xi2 <- Variable(rows = nrow(X2))
      xi3 <- Variable(rows = nrow(X3))
      
      objective <- Minimize(sum_squares(vstack(w1,w2,w3))/2 +  C*(sum(xi1)+sum(xi2)+sum(xi3)))
      constraints <- list(b1+b2+b3 ==0,
                          slack1[,1] == 1,
                          X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
                          X1 %*% (w1 - w3) + b1 - b3 >= 1-slack1[,3],
                          slack2[,2] == 1,
                          X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                          X2 %*% (w2 - w3) + b2 - b3 >= 1-slack2[,3],
                          slack3[,3] == 1,
                          X3 %*% (w3 - w1) + b3 - b1 >= 1-slack3[,1],
                          X3 %*% (w3 - w2) + b3 - b2 >= 1-slack3[,2],
                          xi1 >= max_entries(slack1, axis = 1)-1,
                          xi2 >= max_entries(slack2, axis = 1)-1,
                          xi3 >= max_entries(slack3, axis = 1)-1,
                          xi1 >= (sum_entries(slack1, axis = 1) - min_entries(slack1, axis = 1))/2-1/2,
                          xi2 >= (sum_entries(slack2, axis = 1) - min_entries(slack2, axis = 1))/2-1/2,
                          xi3 >= (sum_entries(slack3, axis = 1) - min_entries(slack3, axis = 1))/2-1/2,
                          xi1 >= sum_entries(slack1, axis = 1)/3 - 1/3,
                          xi2 >= sum_entries(slack2, axis = 1)/3 - 1/3,
                          xi3 >= sum_entries(slack3, axis = 1)/3 - 1/3,
                          slack1 >=0,
                          slack2 >=0,
                          slack3 >=0)
      
      
      Duchi <- Problem(objective, constraints)
      CVXR_Duchi <- solve(Duchi, solver = "MOSEK")
      
      b <- c(CVXR_Duchi$getValue(b1), CVXR_Duchi$getValue(b2), CVXR_Duchi$getValue(b3))
    }
    return(cbind(w,b))
  }
  
  return(solve_w_and_b(X,Y,alpha_sol, tau_sol))
  
}

