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



Duchi_pri_opt <- function(X,y,C, intercept_only = F, w = NULL){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  if(intercept_only == F){
    w <- Variable(rows = m, cols = p)
    b <- Variable(m)
    slack <- Variable(rows = n, cols = m)
    epsilon <- Variable(n)
    t <- Variable(n*m)
    u <- Variable(rows = n*m, cols = m)
    
    objective <- Minimize(sum_squares(w)/2 +  C*sum(epsilon))
    constraints <- list( sum_entries(b) == 0,
                         sum_entries(w, axis = 2) ==0,
                         ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                         vec(t(epsilon%*%matrix(1,ncol=m))) >= t + sum_entries(rep(1/seq(1,m),n)%*%matrix(1,ncol = m)*u,axis=1) - rep(1/seq(1,m),n),
                         t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m))%*%slack,
                         u>=0,
                         slack >=0)
    
    Duchi <- Problem(objective, constraints)
    CVXR_Duchi <- solve(Duchi, solver = "MOSEK")
    
    return( cbind(CVXR_Duchi$getValue(w), CVXR_Duchi$getValue(b)))
  }else{
    b <- Variable(m)
    slack <- Variable(rows = n, cols = m)
    epsilon <- Variable(n)
    t <- Variable(n*m)
    u <- Variable(rows = n*m, cols = m)
    
    objective <- Minimize(C*sum(epsilon))
    constraints <- list(sum_squares(w)/2 + sum_entries(b) == 0,
                         ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                         vec(t(epsilon%*%matrix(1,ncol=m))) >= t + sum_entries(rep(1/seq(1,m),n)%*%matrix(1,ncol = m)*u,axis=1) - rep(1/seq(1,m),n),
                         t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m))%*%slack,
                         u>=0,
                         slack >=0)
    
    Duchi <- Problem(objective, constraints)
    CVXR_Duchi <- solve(Duchi, solver = "MOSEK")
    return(CVXR_Duchi$getValue(b))
  }
}

MDuchi_pri_opt <- function(X,y,C, intercept_only = F, w = NULL){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  if(intercept_only == F){
    w <- Variable(rows = m, cols = p)
    b <- Variable(m)
    slack <- Variable(rows = n, cols = m)
    epsilon <- Variable(n)
    t <- Variable(n*(m-1))
    u <- Variable(rows = n*(m-1), cols = m)
  
    objective <- Minimize(sum_squares(w)/2 +  C*sum(epsilon))
    constraints <- list( sum_entries(b) == 0,
                         sum_entries(w, axis = 2) ==0,
                         ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-Y)*(1-slack),
                         vec(t(epsilon%*%matrix(1,ncol=m-1))) >=  rep(seq(1,m-1),n)/rep(seq(2,m),n)*t + 1/rep(seq(2,m),n)*sum_entries(u,axis=1),
                         t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m-1))%*%slack,
                         u>=0,
                         slack >=0)
    
    M_Duchi <- Problem(objective, constraints)
      CVXR_M_Duchi <- solve(M_Duchi, solver = "MOSEK")
  return( cbind(CVXR_M_Duchi$getValue(w), CVXR_M_Duchi$getValue(b)))
  }else{
    b <- Variable(m)
    slack <- Variable(rows = n, cols = m)
    epsilon <- Variable(n)
    t <- Variable(n*(m-1))
    u <- Variable(rows = n*(m-1), cols = m)
  
    objective <- Minimize(sum_squares(w)/2 +  C*sum(epsilon))
    constraints <- list( sum_entries(b) == 0,
                         ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-Y)*(1-slack),
                         vec(t(epsilon%*%matrix(1,ncol=m-1))) >=  rep(seq(1,m-1),n)/rep(seq(2,m),n)*t + 1/rep(seq(2,m),n)*sum_entries(u,axis=1),
                         t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m-1))%*%slack,
                         u>=0,
                         slack >=0)
    
    M_Duchi <- Problem(objective, constraints)
    CVXR_M_Duchi <- solve(M_Duchi, solver = "MOSEK")
    
    return(CVXR_M_Duchi$getValue(b))
  }
}


OVA_pri_opt <- function(X,y,C){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  b <- Variable(rows = m)
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(slack))
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (2*Y-1)) >=1- slack,
                      slack >=0)
  
  OVA <- Problem(objective, constraints)
  CVXR_OVA <- solve(OVA, solver="MOSEK")

  return(cbind(CVXR_OVA$getValue(w),CVXR_OVA$getValue(b)))
}


LLW_pri_opt <- function(X,y,C){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  b <- Variable(rows = m)
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(slack))
  
  constraints <- list(sum_entries(b) == 0,
                      sum_entries(w, axis = 2) ==0,
                      ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1- slack) * (1-Y),
                      slack >=0)
  
  LLW <- Problem(objective, constraints)
  CVXR_LLW <- solve(LLW)
  
  return(cbind(CVXR_LLW$getValue(w),CVXR_LLW$getValue(b)))
}


MSVM7_pri_opt <- function(X,y,C, base_class = NULL){
  class_idx <- sort(unique(y))
  m <- length(class_idx)
  if(is.null(base_class)){
    Y <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y==id)})
  }else{
    Y <- sapply(class_idx[class_idx!=base_class], function(id){as.numeric(y==id)})
  }
  
  n <- nrow(X)
  p <- ncol(X)
  w <- Variable(rows = m-1, cols = p)
  b <- Variable(rows = m-1)
  slack <- Variable(rows = n, cols = m-1)
  
  objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(slack))
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(1- slack),
                      sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >=(1- sum_entries(slack, axis=1))*sum_entries(Y,axis=1),
                      slack >=0)
  
  MSVM7 <- Problem(objective, constraints)
  CVXR_MSVM7 <- solve(MSVM7, solver = "MOSEK")
  
  return(cbind(CVXR_MSVM7$getValue(w),CVXR_MSVM7$getValue(b)))
}



MSVM8_pri_opt <- function(X,y,C){
  
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  b <- Variable(rows = m)
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(sum_squares(w)/2 +  C*sum_entries(slack))
  
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(1- slack),
                      sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >=1- sum_entries(slack, axis=1),
                      slack >=0)
  
  MSVM8 <- Problem(objective, constraints)
  CVXR_MSVM8 <- solve(MSVM8, solver = "MOSEK")
  
  return(cbind(CVXR_MSVM8$getValue(w),CVXR_MSVM8$getValue(b)))
}


New1_pri_opt <- function(X,y,C){
  
  X1 <- X[y==1,,drop=F]
  X2 <- X[y==2,,drop=F]
  X3 <- X[y==3,,drop=F]
  p <- ncol(X)
  w1 <- Variable(p)
  w2 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  
  slack1 <- Variable(rows = nrow(X1), cols = 2)
  slack2 <- Variable(rows = nrow(X2), cols = 2)
  slack3 <- Variable(rows = nrow(X3), cols = 2)
  xi1 <- Variable(rows = nrow(X1))
  xi2 <- Variable(rows = nrow(X2))
  
  objective <- Minimize(sum_squares(vstack(w1,w2))/2 +  C*(sum(xi1)+sum(xi2)+sum_entries(slack3)))
  constraints <- list(X1 %*% w1 + b1 >= 1-xi1,
                      slack1[,1] == 1,
                      X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
                      X2 %*% w2 + b2 >= 1-xi2,
                      slack2[,2] == 1,
                      X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
                      X3 %*% w1 + b1 <= -1 + slack3[,1],
                      X3 %*% w2 + b2 <= -1 + slack3[,2],
                      xi1 >= max_entries(slack1, axis = 1)-1,
                      xi2 >= max_entries(slack2, axis = 1)-1,
                      xi1 >= sum_entries(slack1, axis = 1)/2 - 1/2,
                      xi2 >= sum_entries(slack2, axis = 1)/2 - 1/2,
                      slack1 >=0,
                      slack2 >=0,
                      slack3 >=0)
  New1 <- Problem(objective, constraints)
  CVXR_New1 <- solve(New1, solver = "MOSEK")
  
  beta1 <- c(CVXR_New1$getValue(w1),CVXR_New1$getValue(b1))
  beta2 <- c(CVXR_New1$getValue(w2),CVXR_New1$getValue(b2))
  New1_primary_beta <- rbind(beta1,beta2)
  
  return(New1_primary_beta)
}


New3_pri_opt <- function(X,y,C){
  
  X1 <- X[y==1,,drop=F]
  X2 <- X[y==2,,drop=F]
  X3 <- X[y==3,,drop=F]
  p <- ncol(X)
  w1 <- Variable(p)
  w2 <- Variable(p)
  b1 <- Variable(1)
  b2 <- Variable(1)
  
  slack1 <- Variable(rows = nrow(X1), cols = 2)
  slack2 <- Variable(rows = nrow(X2), cols = 2)
  slack3 <- Variable(rows = nrow(X3), cols = 2)
  xi1 <- Variable(rows = nrow(X1))
  xi2 <- Variable(rows = nrow(X2))
  
  objective <- Minimize(sum_squares(vstack(w1,w2))/2 +  C*(sum(xi1)+sum(xi2)+sum_entries(slack3)))
  constraints <- list(
    X1 %*% (w1 - w2) + b1 - b2 >= 1-slack1[,2],
    X1 %*% w1 + b1 >= 1 - xi1,
    X2 %*% (w2 - w1) + b2 - b1 >= 1-slack2[,1],
    X2 %*% w2 + b2 >= 1 - xi2,
    X3 %*% w1 + b1 <= -1 + slack3[,1],
    X3 %*% w2 + b2 <= -1 + slack3[,2],
    xi1 >= max_entries(slack1, axis = 1)/2,
    xi2 >= max_entries(slack2, axis = 1)/2,
    slack1 >=0,
    slack2 >=0,
    slack3 >=0)
  
  New3 <- Problem(objective, constraints)
  CVXR_New3 <- solve(New3, solver = "MOSEK")
  
  beta1 <- c(CVXR_New3$getValue(w1),CVXR_New3$getValue(b1))
  beta2 <- c(CVXR_New3$getValue(w2),CVXR_New3$getValue(b2))
  New3_primary_beta <- rbind(beta1,beta2)
  
  return(New3_primary_beta)
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
  y <- sort(sample(c(1,2,3), size = n-6, replace = T, prob =  1/rep(3,3)))
  n1 <- sum(y==1)+1
  n2 <- sum(y==2)+1
  n3 <- sum(y==3)+1
  X1 <- matrix(mvrnorm(n1, sep*c(0,2), diag(v,nrow=2)), nrow=n1)
  X2 <- matrix(mvrnorm(n2, sep*c(sqrt(3),-1), diag(v,nrow=2)), nrow=n2)
  X3 <- matrix(mvrnorm(n3, sep*c(-sqrt(3),-1), diag(v,nrow=2)), nrow=n3)
  ########  Plan II     #############
  # X1 <- matrix(mvrnorm(n1, sep*c(-sqrt(3),-1),  matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2)), nrow=n1)
  # X2 <- matrix(mvrnorm(n2, sep*c(sqrt(3),-1),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2)), nrow=n2)
  # X3 <- matrix(mvrnorm(n3, sep*c(0,2), diag(c(8,1),nrow=2)), nrow=n3)

  ########  Plan III     #############
  # X1 <- mvrnorm(n1, sep*c(0,2), diag(c(1,8),nrow=2))
  # X2 <- mvrnorm(n2, sep*c(sqrt(3),-1),   matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2))
  # X3 <- mvrnorm(n3, sep*c(-sqrt(3),-1), matrix(c(4.5,3.5,3.5,4.5), nrow= 2))
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



