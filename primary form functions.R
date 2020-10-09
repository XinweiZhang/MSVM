library(ggplot2)
WW_pri_opt <- function(X,y,C = 1, lambda = 1, intercept = T){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  if(intercept == T){
    b <- Variable(m)
  }else{
    b <- rep(0,m)
  }
  slack <- Variable(rows = n, cols = m)
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum_entries((1-Y)*slack))
  
  constraints <- list( sum_entries(b) == 0,
                       sum_entries(w, axis = 2) ==0,
                       ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >=1-slack,
                       slack >=0)
  
  WW <- Problem(objective, constraints)
  CVXR_WW <- solve(WW, solver = "MOSEK")
  
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_WW$getValue(w),CVXR_WW$getValue(b)), X = X, y = y, value = CVXR_WW$value, status = CVXR_WW$status)
    class(fit) <- "msvm"
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_WW$getValue(w),0), X = X, y = y, value = CVXR_WW$value, status = CVXR_WW$status)
    class(fit) <- "msvm"
    return(fit)
  }
}

CS_pri_opt <- function(X,y,C = 1, lambda = 1, intercept = T){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  if(intercept == T){
    b <- Variable(m)
  }else{
    b <- rep(0,m)
  }
  slack <- Variable(rows = n, cols = m)
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum_entries(max_entries((1-Y)*slack, axis = 1)))
  
  constraints <- list( sum_entries(b) == 0,
                       sum_entries(w, axis = 2) ==0,
                       ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-slack),
                       slack >=0)
  
  CS <- Problem(objective, constraints)
  CVXR_CS <- solve(CS, solver = "MOSEK")
  
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_CS$getValue(w),CVXR_CS$getValue(b)), X = X, y = y, value = CVXR_CS$value, status = CVXR_CS$status)
    class(fit) <- "msvm"
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_CS$getValue(w),0), X = X, y = y, value = CVXR_CS$value, status = CVXR_CS$status)
    class(fit) <- "msvm"
    return(fit)
  }
}



Duchi_pri_opt <- function(X,y,C = 1, lambda = 1, intercept = T, intercept_only = F, w = NULL){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  if(intercept_only == F){
    w <- Variable(rows = m, cols = p)
    if(intercept == T){
      b <- Variable(m)
    }else{
      b <- rep(0,m)
    }
    slack <- Variable(rows = n, cols = m)
    epsilon <- Variable(n)
    t <- Variable(n*m)
    u <- Variable(rows = n*m, cols = m)
    
    objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum(epsilon))
    constraints <- list( sum_entries(b) == 0,
                         sum_entries(w, axis = 2) ==0,
                         ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                         vec(t(epsilon%*%matrix(1,ncol=m))) >= t + sum_entries(rep(1/seq(1,m),n)%*%matrix(1,ncol = m)*u,axis=1) - rep(1/seq(1,m),n),
                         t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m))%*%slack,
                         u>=0,
                         slack >=0)
    
    Duchi <- Problem(objective, constraints)
    CVXR_Duchi <- solve(Duchi, solver = "MOSEK")
    
    if(intercept == T){
      fit <- list(beta = cbind(CVXR_Duchi$getValue(w),CVXR_Duchi$getValue(b)), X = X, y = y, value = CVXR_Duchi$value, status = CVXR_Duchi$status)
      class(fit) <- "msvm"
      return(fit)
    }else{  
      fit <- list(beta = cbind(CVXR_Duchi$getValue(w),0), X = X, y = y, value = CVXR_Duchi$value, status = CVXR_Duchi$status)
      class(fit) <- "msvm"
      return(fit)
    }
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

MDuchi_pri_opt <- function(X, y, C = 1, lambda = 1, intercept = T, intercept_only = F, w = NULL){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  if(intercept_only == F){
    w <- Variable(rows = m, cols = p)
    if(intercept == T){
      b <- Variable(m)
    }else{
      b <- rep(0,m)
    }
    slack <- Variable(rows = n, cols = m)
    epsilon <- Variable(n)
    t <- Variable(n*(m-1))
    u <- Variable(rows = n*(m-1), cols = m)
  
    objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum(epsilon))
    constraints <- list( sum_entries(b) == 0,
                         sum_entries(w, axis = 2) ==0,
                         ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-Y)*(1-slack),
                         vec(t(epsilon%*%matrix(1,ncol=m-1))) >=  rep(seq(1,m-1),n)/rep(seq(2,m),n)*t + 1/rep(seq(2,m),n)*sum_entries(u,axis=1),
                         t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m-1))%*%slack,
                         u>=0,
                         slack >=0)
    
    M_Duchi <- Problem(objective, constraints)
    CVXR_M_Duchi <- solve(M_Duchi, solver = "MOSEK")
    if(intercept == T){
      fit <- list(beta = cbind(CVXR_M_Duchi$getValue(w),CVXR_M_Duchi$getValue(b)), X = X, y = y, value = CVXR_M_Duchi$value, status = CVXR_M_Duchi$status)
      class(fit) <- "msvm"
      return(fit)
    }else{  
      fit <- list(beta = cbind(CVXR_M_Duchi$getValue(w),0), X = X, y = y, value = CVXR_M_Duchi$value, status = CVXR_M_Duchi$status)
      class(fit) <- "msvm"
      return(fit)
    }
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


OVA_pri_opt <- function(X, y, C = 1, lambda = 1, intercept = T){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  if(intercept == T){
    b <- Variable(m)
  }else{
    b <- rep(0,m)
  }
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum_entries(slack))
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (2*Y-1)) >=1- slack,
                      slack >=0)
  
  OVA <- Problem(objective, constraints)
  CVXR_OVA <- solve(OVA, solver="MOSEK")
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_OVA$getValue(w),CVXR_OVA$getValue(b)), X = X, y = y, value = CVXR_OVA$value, status = CVXR_OVA$status)
    class(fit) <- "msvm"
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_OVA$getValue(w),0), X = X, y = y, value = CVXR_OVA$value, status = CVXR_OVA$status)
    class(fit) <- "msvm"
    return(fit)
  }
}


LLW_pri_opt <- function(X, y, C = 1, lambda = 1, intercept = T){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  if(intercept == T){
    b <- Variable(m)
  }else{
    b <- rep(0,m)
  }
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum_entries(slack))
  
  constraints <- list(sum_entries(b) == 0,
                      sum_entries(w, axis = 2) ==0,
                      ((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1- slack) * (1-Y),
                      slack >=0)
  
  LLW <- Problem(objective, constraints)
  CVXR_LLW <- solve(LLW)
  
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_LLW$getValue(w),CVXR_LLW$getValue(b)), X = X, y = y, value = CVXR_LLW$value, status = CVXR_LLW$status)
    class(fit) <- "msvm"
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_LLW$getValue(w),0), X = X, y = y, value = CVXR_LLW$value, status = CVXR_LLW$status)
    class(fit) <- "msvm"
    return(fit)
  }
}


MSVM7_pri_opt <- function(X, y, C = 1, lambda = 1, intercept = T,  base_class = NULL){
  class_idx <- sort(unique(y))
  m <- length(class_idx)
  if(is.null(base_class)){
    base_class <- m
    Y <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y==id)})
  }else{
    Y <- sapply(class_idx[class_idx!=base_class], function(id){as.numeric(y==id)})
  }
  
  n <- nrow(X)
  p <- ncol(X)
  w <- Variable(rows = m-1, cols = p)
  if(intercept == T){
    b <- Variable(m-1)
  }else{
    b <- rep(0,m-1)
  }
  slack <- Variable(rows = n, cols = m-1)
  
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum_entries(slack))
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(0 - slack),
                      sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >= (1 - sum_entries(slack, axis=1))*sum_entries(Y,axis=1),
                      slack >=0)
  
  MSVM7 <- Problem(objective, constraints)
  CVXR_MSVM7 <- solve(MSVM7, solver = "MOSEK")
  
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_MSVM7$getValue(w),CVXR_MSVM7$getValue(b)), X = X, y = y, value = CVXR_MSVM7$value, status = CVXR_MSVM7$status)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_MSVM7$getValue(w),0), X = X, y = y, value = CVXR_MSVM7$value, status = CVXR_MSVM7$status)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }
}



  MSVM8_pri_opt <- function(X, y, C = 1, lambda = 1, intercept = T){
  
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  w <- Variable(rows = m, cols = p)
  if(intercept == T){
    b <- Variable(m)
  }else{
    b <- rep(0,m)
  }
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*sum_entries(slack))
  
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(1 - slack),
                      sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >= 1 - sum_entries(slack, axis=1),
                      slack >=0)
  
  MSVM8 <- Problem(objective, constraints)
  CVXR_MSVM8 <- solve(MSVM8, solver = "MOSEK")
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_MSVM8$getValue(w),CVXR_MSVM8$getValue(b)), X = X, y = y, value = CVXR_MSVM8$value, status = CVXR_MSVM8$status)
    class(fit) <- "msvm"
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_MSVM8$getValue(w),0), X = X, y = y, value = CVXR_MSVM8$value, status = CVXR_MSVM8$status)
    class(fit) <- "msvm"
    return(fit)
  }
}


New1_pri_opt <- function(X, y, C = 1, lambda = 1, intercept = T, base_class = NULL){
  
  class_idx <- sort(unique(y))
  m = length(class_idx)
  if(is.null(base_class)){
    base_class <- m
  }
  rel_class <- class_idx[class_idx!=base_class]
  n <- nrow(X)
  X_p1 <- X[y!=base_class,,drop = F]
  X_p2 <- X[y==base_class,,drop = F]
  n_p1 <- nrow(X_p1)
  n_p2 <- nrow(X_p2)
  y_p1 <- y[y!=base_class]
  Y_p1 <- sapply(rel_class, function(id){y_p1==id})
  
  p <- ncol(X)
  w <- Variable(rows = m-1, cols = p)
  if(intercept == T){
    b <- Variable(m-1)
  }else{
    b <- rep(0,m-1)
  }
  
  slack_p1 <- Variable(rows = n_p1, cols = m-1)
  slack_p2 <- Variable(rows = n_p2, cols = m-1)
  epsilon_p1 <- Variable(n_p1)
  t <- Variable(n_p1*(m-1))
  u <- Variable(rows = n_p1*(m-1), cols = m-1)
  
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*(sum(epsilon_p1)+sum_entries(slack_p2)))
  constraints <- list(((X_p1%*%t(w) + matrix(1, nrow = n_p1,ncol = 1)%*%t(b))*Y_p1)%*%matrix(1, nrow = m-1, ncol = m-1) - (X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b)) >= 1-slack_p1,
                      vec(t(epsilon_p1%*%matrix(1, ncol = m-1))) >= t + sum_entries(rep(1/seq(1, m-1), n_p1)%*%matrix(1, ncol = m-1)*u, axis=1) - rep(1/seq(1,m-1),n_p1),
                      t%*%matrix(1,nrow=1,ncol=m-1) + u >=  (diag(1,n_p1)%x%matrix(1,nrow=m-1))%*%slack_p1,
                      u>=0,
                      sum_entries((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1, axis = 1) >= 1-epsilon_p1,
                      X_p2%*%t(w) + matrix(1,nrow=n_p2,ncol =1)%*%t(b) <= 0 + slack_p2,
                      slack_p1 >=0,
                      slack_p2 >=0)
  
  New1 <- Problem(objective, constraints)
  CVXR_New1 <- solve(New1, solver = "MOSEK")
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_New1$getValue(w),CVXR_New1$getValue(b)), X = X, y = y, value = CVXR_New1$value, status = CVXR_New1$status)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_New1$getValue(w),0), X = X, y = y, value = CVXR_New1$value, status = CVXR_New1$status)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }
}


New3_pri_opt <- function(X, y, C = 1, lambda = 1, intercept = T, base_class = NULL){
  
  
  class_idx <- sort(unique(y))
  m = length(class_idx)
  if(is.null(base_class)){
    base_class <- m
  }
  rel_class <- class_idx[class_idx!=base_class]
  n <- nrow(X)
  X_p1 <- X[y!=base_class,,drop = F]
  X_p2 <- X[y==base_class,,drop = F]
  n_p1 <- nrow(X_p1)
  n_p2 <- nrow(X_p2)
  y_p1 <- y[y!=base_class]
  Y_p1 <- sapply(rel_class, function(id){y_p1==id})
  
  p <- ncol(X)
  w <- Variable(rows = m-1, cols = p)
  if(intercept == T){
    b <- Variable(m-1)
  }else{
    b <- rep(0,m-1)
  }
  
  slack_p1 <- Variable(rows = n_p1, cols = m-1)
  slack_p2 <- Variable(rows = n_p2, cols = m-1)
  epsilon_p1 <- Variable(n_p1)
  t <- Variable(n_p1*(m-2))
  u <- Variable(rows = n_p1*(m-2), cols = m-1)
  
  
  objective <- Minimize(lambda*sum_squares(w)/2 +  C*(sum(epsilon_p1)+sum_entries(slack_p2)))
  constraints <- list(((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1)%*%matrix(1, nrow = m-1, ncol = m-1) - (X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b)) >= (1-Y_p1)*(1-slack_p1),
                      vec(t(epsilon_p1%*%matrix(1,ncol=m-2))) >=  rep(seq(1,m-2),n_p1)/rep(seq(2,m-1),n_p1)*t + 1/rep(seq(2,m-1),n_p1)*sum_entries(u,axis=1),
                      t%*%matrix(1,nrow=1,ncol=m-1) + u >=  (diag(1,n_p1)%x%matrix(1,nrow=m-2))%*%slack_p1,
                      u>=0,
                      sum_entries((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1, axis = 1) >= 1-epsilon_p1,
                      X_p2%*%t(w) + matrix(1,nrow=n_p2,ncol =1)%*%t(b) <= -0 + slack_p2,
                      slack_p1 >=0,
                      slack_p2 >=0)
  
  
  New3 <- Problem(objective, constraints)
  CVXR_New3 <- solve(New3, solver = "MOSEK")
  
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_New3$getValue(w),CVXR_New3$getValue(b)), X = X, y = y, value = CVXR_New3$value, status = CVXR_New3$status)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_New3$getValue(w),0), X = X, y = y, value = CVXR_New3$value, status = CVXR_New3$status)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }
}


predict.msvm <- function(model, X_test, rule = "simple_max"){
  if(rule == "simple_max"){
    lpred <- round(cbind(X_test,1)%*%t(model$beta),7)
    # y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
    y_pred <-  apply(lpred, MARGIN = 1, FUN= function(t){
      ifelse(sum(abs(max(t) - t)<=1e-3)==1, which.max(t), 0)
    })
  }else if(rule == "dagger")
  {
    lpred <- matrix(0, nrow = nrow(X_test), ncol = length(unique(model$y)))
    base_class <- attr(model, "base_class") 
    rel_class <- setdiff(unique(model$y), base_class)
    lpred[, -c(base_class)] <-cbind(X_test,1)%*%t(model$beta)
    lpred[, base_class] <- 1-rowSums(sapply(rel_class, function(idx){apply(as.matrix(cbind(lpred[,idx],0)), MARGIN = 1, FUN = max)}))
    # y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
    y_pred <-  apply(round(lpred,7), MARGIN = 1, FUN= function(t){
      ifelse(sum(abs(max(t) - t)<=1e-3)==1, which.max(t), 0)
    })
  }else if(rule == "dagger2")
  {
    lpred <- matrix(0, nrow = nrow(X_test), ncol = length(unique(model$y)))
    base_class <- attr(model, "base_class") 
    rel_class <- setdiff(unique(model$y), base_class)
    lpred[, -c(base_class)] <-cbind(X_test,1)%*%t(model$beta)
    lpred[, base_class] <- 1/2-rowSums(sapply(rel_class, function(idx){apply(as.matrix(cbind(lpred[,idx]+1/2,0)), MARGIN = 1, FUN = max)}))
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "dagger_MSVM7")
  {
    lpred <- cbind(lpred, 1 - apply(as.matrix(cbind(lpred[,1]+1,0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(lpred[,2]+1,0)), MARGIN = 1, FUN = max))
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }
  else if(rule == "0_max"){
    lpred <- matrix(0, nrow = nrow(X_test), ncol = length(unique(model$y)))
    base_class <- attr(model, "base_class") 
    rel_idx <- setdiff(unique(model$y), base_class)
    lpred[, -c(base_class)] <-cbind(X_test,1)%*%t(model$beta)
    lpred[, base_class] <- 0
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }
  return(y_pred)
}

plot_1D_decision_boundary <- function(model, X = NULL, y = NULL, title = NULL, np_resolution = 500, xlim =NULL, ylim = NULL, rule = "simple_max", print = F){
  
  if(is.null(X)|| is.null(y))
  {
    X_dat <- as.data.frame(cbind(model$X,0))
    colnames(X_dat) <- c("x","y")
    X_dat$class <- as.factor(model$y)  
  }
  else{
    X_dat <- as.data.frame(X)
    colnames(X_dat) <- c("x","y")
    X_dat$class <- as.factor(y)  
  }
  
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$x,X_dat$y))-.2, ceiling(max(X_dat$x,X_dat$y))+.2)
  }
  if(is.null(ylim)){
    ylim = c(-1,1)
  }
  
  nd.x = seq(from = xlim[1], to =  xlim[2], length.out = np_resolution)
  nd.y = seq(from = ylim[1], to =  ylim[2], length.out = np_resolution)
  nd = expand.grid(x = nd.x, y = nd.y)
  model$beta <- cbind(model$beta[1,], 0, model$beta[2,])
  nd$class <- as.factor(predict(model, as.matrix(nd), rule = rule))
  
  colorfun <- function(n,l=65,c=100) { hues = seq(15, 375, length=n+1); hcl(h=hues, l=l, c=c)[1:n] } # default ggplot2 colours
  colors <- colorfun(length(unique(model$y)))
  colorslight <- colorfun(length(unique(model$y)),l=90,c=50)
  names(colorslight) <- unique(model$y) 
  plt <- ggplot(nd, aes(x=x, y=y)) +
    geom_raster(data=nd, aes(x = x, y = y, fill = factor(class)),alpha=0.7,show.legend=FALSE) +
    geom_contour(data=nd, aes(x = x, y = y, z= as.numeric(nd$class)), colour="red2", alpha=0.5, breaks=c(1.5,2.5)) +
    geom_point(data = X_dat, size = 2, aes(pch = class,  colour = class)) +
    scale_x_continuous(limits = xlim, expand=c(0,0)) +
    scale_y_continuous(limits = ylim, expand=c(0,0)) +
    scale_fill_manual(values=colorslight,guide=F) + 
    ggtitle(title)
  if(print == T){
    suppressWarnings(print(plt))
  }
  return(plt)
}


data_generate <- function(n, sep =  1,  v = 1.5^2){
  y <- sort(sample(c(1,2,3), size = n-3, replace = T, prob =  1/rep(3,3)))
  n1 <- sum(y==1)+1
  n2 <- sum(y==2)+1
  n3 <- sum(y==3)+1
  
  X1 <- matrix(mvrnorm(n1, sep*c(0,2), diag(v,nrow=2)), nrow=n1)
  X2 <- matrix(mvrnorm(n2, sep*c(sqrt(3),-1), diag(v,nrow=2)), nrow=n2)
  X3 <- matrix(mvrnorm(n3, sep*c(-sqrt(3),-1), diag(v,nrow=2)), nrow=n3)
  # ########  Plan II     #############
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

plot_decision_boundary <- function(model, X = NULL, y = NULL, title = NULL, np_resolution = 500, xlim =NULL, ylim = NULL, rule = "simple_max", print = F){
  
  if(is.null(X)|| is.null(y))
  {
    X_dat <- as.data.frame(model$X)
    colnames(X_dat) <- c("x","y")
    X_dat$class <- as.factor(model$y)  
  }
  else{
    X_dat <- as.data.frame(X)
    colnames(X_dat) <- c("x","y")
    X_dat$class <- as.factor(y)  
  }
  
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$x,X_dat$y))-.2, ceiling(max(X_dat$x,X_dat$y))+.2)
  }
  if(is.null(ylim)){
    ylim = c( floor(min(X_dat$x,X_dat$y))-.2, ceiling(max(X_dat$x,X_dat$y))+.2)
  }
  
  nd.x = seq(from = xlim[1], to =  xlim[2], length.out = np_resolution)
  nd.y = seq(from = ylim[1], to =  ylim[2], length.out = np_resolution)
  nd = expand.grid(x = nd.x, y = nd.y)
  
  nd$class <- as.factor(predict(model, as.matrix(nd), rule = rule))
  
  colorfun <- function(n,l=65,c=100) { hues = seq(15, 375, length=n+1); hcl(h=hues, l=l, c=c)[1:n] } # default ggplot2 colours
  colors <- colorfun(length(unique(model$y)))
  colorslight <- colorfun(length(unique(model$y)),l=90,c=50)
  names(colorslight) <- sort(unique(model$y)) 
  plt <- ggplot(nd, aes(x=x, y=y)) +
    geom_raster(data=nd, aes(x = x, y = y, fill = factor(class)),alpha=0.7,show.legend=FALSE) +
    geom_contour(data=nd, aes(x = x, y = y, z= as.numeric(nd$class)), colour="red2", alpha=0.5, breaks=c(1.5,2.5)) +
    geom_point(data = X_dat, size = 2, aes(pch = class,  colour = class)) +
    scale_x_continuous(limits = xlim, expand=c(0,0)) +
    scale_y_continuous(limits = ylim, expand=c(0,0)) +
    scale_fill_manual(values=colorslight,guide=F) + 
    ggtitle(title)
  if(print == T){
    suppressWarnings(print(plt))
  }
  return(plt)
}



plot_nearest_neighbour_decision_boundary <- function(X, y, dis = "L2", title = NULL, np_resolution = 500, xlim =NULL, ylim = NULL, dagger_rule = T){
  X_dat <- as.data.frame(X)
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$V1))-1, ceiling(max(X_dat$V1))+1)
  }
  if(is.null(ylim)){
    ylim = c( floor(min(X_dat$V2))-1, ceiling(max(X_dat$V2))+1)
  }
  
  nd.x = seq(from = xlim[1], to =  xlim[2], length.out = np_resolution)
  nd.y = seq(from = ylim[1], to =  ylim[2], length.out = np_resolution)
  nd = expand.grid(Var1 = nd.x, Var2 = nd.y)
  
  if(dis == "L2"){
    prd <- apply(nd, MARGIN = 1, function(s,X){which.min(apply(X,MARGIN=1,FUN = function(t){sqrt(sum((s-t)^2))}))},X)
  }else{
    prd <- apply(nd, MARGIN = 1, function(s,X){which.min(apply(X,MARGIN=1,FUN = function(t){sum(abs(s-t))}))},X)
  }
  op <-  par(mfrow = c(1,1), mar=c(5.1, 4.1, 4.1, 7), xpd=TRUE)
  plot(X_dat$V1, X_dat$V2, col = as.factor(y), ylim=ylim, xlim=xlim, xlab ="X", ylab = "Y", main = title)
  
  contour(x = nd.x, y = nd.y, z = matrix(prd, nrow = np_resolution, ncol = np_resolution), 
          levels = unique(y), add = TRUE, drawlabels = FALSE)
  
  
  legend("topright", inset=c(-0.3,0),legend = sapply(unique(y), function(x){paste("Class ",x)}), col= unique(y), pch = 1)
  par(op)
}

four_class_data_generate <- function(n, sep =  1,  v = 1.5^2){
  y <- sort(sample(c(1,2,3,4), size = n-4, replace = T, prob =  1/rep(4,4)))
  n1 <- sum(y==1)+1
  n2 <- sum(y==2)+1
  n3 <- sum(y==3)+1
  n4 <- sum(y==4)+1
  # 
  X1 <- matrix(mvrnorm(n1, sep*c(2,2), diag(v,nrow=2)), nrow=n1)
  X2 <- matrix(mvrnorm(n2, sep*c(2,-2), diag(v,nrow=2)), nrow=n2)
  X3 <- matrix(mvrnorm(n3, sep*c(-2,2), diag(v,nrow=2)), nrow=n3)
  X4 <- matrix(mvrnorm(n4, sep*c(-2,-2), diag(v,nrow=2)), nrow=n4)
  ########  Plan II     #############
  # X1 <- matrix(mvrnorm(n1, sep*c(-2,-2),  matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2)), nrow=n1)
  # X2 <- matrix(mvrnorm(n2, sep*c(2,-2),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2)), nrow=n2)
  # X3 <- matrix(mvrnorm(n3, sep*c(2,2), matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2)), nrow=n3)
  # X4 <- matrix(mvrnorm(n4, sep*c(-2,2),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2)), nrow=n4)

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
  
  X <- rbind(X1,X2,X3,X4)
  y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)),rep(4,nrow(X4)))
  return(list(X = X,y = y))
}

