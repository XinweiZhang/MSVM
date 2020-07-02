WW_kernel_pri_opt <- function(X,y,C, kernel){
  
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  v <- Variable(n,m)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  
  # objective <- Minimize(sum_squares(chol(K+diag(1e-5,nrow=nrow(K)))%*%v)/2  +  C*sum_entries(slack))
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2 +  C*sum_entries(slack))
  
  constraints <- list( ((K%*%v + matrix(1, nrow=n, ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (K%*%v  + matrix(1, nrow=n, ncol =1)%*%t(b)) >= (1-Y)*(1-slack),
                       sum_entries(b)*matrix(1,nrow=n) + K %*% sum_entries(v, axis = 1) == 0,
                       slack >=0)
  
  
  WW <- Problem(objective, constraints)
  CVXR_WW <- solve(WW)
  CVXR_WW_v <- CVXR_WW$getValue(v)
  CVXR_WW_b <- CVXR_WW$getValue(b)
  
  model <- list(v = CVXR_WW_v, b = as.numeric(CVXR_WW_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}



CS_kernel_pri_opt <- function(X,y,C, kernel){
  
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  v <- Variable(n,m)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2  + C*sum_entries(max_entries((1-Y)*slack, axis = 1)))
  constraints <- list( ((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (K%*%v  + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                       sum_entries(b)*matrix(1,nrow=n) + K %*% sum_entries(v, axis = 1) == 0,
                       slack >=0)
  
  
  CS <- Problem(objective, constraints)
  CVXR_CS <- solve(CS, solver = "MOSEK")
  CVXR_CS_v <- CVXR_CS$getValue(v)
  CVXR_CS_b <- CVXR_CS$getValue(b)
  
  model <- list(v = CVXR_CS_v, b = as.numeric(CVXR_CS_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}




Duchi_kernel_pri_opt <- function(X,y,C, kernel){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  v <- Variable(n,m)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  epsilon <- Variable(n)
  t <- Variable(n*m)
  u <- Variable(rows = n*m, cols = m)
    
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2+  C*sum(epsilon))
    constraints <- list(((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (K%*%v  + matrix(1,nrow=n,ncol =1)%*%t(b)) >= 1-slack,
                         sum_entries(b)*matrix(1,nrow=n) + K %*% sum_entries(v, axis = 1) == 0,
                         vec(t(epsilon%*%matrix(1,ncol=m))) >= t + sum_entries(rep(1/seq(1,m),n)%*%matrix(1,ncol = m)*u,axis=1) - rep(1/seq(1,m),n),
                         t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m))%*%slack,
                         u>=0,
                         slack >=0)
    
    Duchi <- Problem(objective, constraints)
    CVXR_Duchi <- solve(Duchi, solver = "MOSEK")
    CVXR_Duchi_v <- CVXR_Duchi$getValue(v)
    CVXR_Duchi_b <- CVXR_Duchi$getValue(b)
    
    model <- list(v = CVXR_Duchi_v, b = as.numeric(CVXR_Duchi_b), kernel = kernel, X = X, y= y)
    class(model) <- "msvm_kernel"
    return(model)
}



MDuchi_kernel_pri_opt <- function(X,y,C, kernel){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  v <- Variable(n,m)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  epsilon <- Variable(n)
  t <- Variable(n*(m-1))
  u <- Variable(rows = n*(m-1), cols = m)
  
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2 +  C*sum(epsilon))
  constraints <- list(((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (K%*%v  + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-Y)*(1-slack),
                       sum_entries(b)*matrix(1,nrow=n) + K %*% sum_entries(v, axis = 1) == 0,
                       vec(t(epsilon%*%matrix(1,ncol=m-1))) >=  rep(seq(1,m-1),n)/rep(seq(2,m),n)*t + 1/rep(seq(2,m),n)*sum_entries(u,axis=1),
                       t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m-1))%*%slack,
                       u>=0,
                       slack >=0)
  
  MDuchi <- Problem(objective, constraints)
  CVXR_MDuchi <- solve(MDuchi, solver = "MOSEK")
  CVXR_MDuchi_v <- CVXR_MDuchi$getValue(v)
  CVXR_MDuchi_b <- CVXR_MDuchi$getValue(b)
  
  model <- list(v = CVXR_MDuchi_v, b = as.numeric(CVXR_MDuchi_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}



MDuchi_kernel_pri_opt <- function(X,y,C, kernel){
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  v <- Variable(n,m)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  epsilon <- Variable(n)
  t <- Variable(n*(m-1))
  u <- Variable(rows = n*(m-1), cols = m)
  
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2 +  C*sum(epsilon))
  constraints <- list(((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*Y)%*%matrix(1, nrow = m, ncol = m) - (K%*%v  + matrix(1,nrow=n,ncol =1)%*%t(b)) >= (1-Y)*(1-slack),
                      sum_entries(b)*matrix(1,nrow=n) + K %*% sum_entries(v, axis = 1) == 0,
                      vec(t(epsilon%*%matrix(1,ncol=m-1))) >=  rep(seq(1,m-1),n)/rep(seq(2,m),n)*t + 1/rep(seq(2,m),n)*sum_entries(u,axis=1),
                      t%*%matrix(1,nrow=1,ncol=m) + u >=  (diag(1,n)%x%matrix(1,nrow=m-1))%*%slack,
                      u>=0,
                      slack >=0)
  
  MDuchi <- Problem(objective, constraints)
  CVXR_MDuchi <- solve(MDuchi, solver = "MOSEK")
  CVXR_MDuchi_v <- CVXR_MDuchi$getValue(v)
  CVXR_MDuchi_b <- CVXR_MDuchi$getValue(b)
  
  model <- list(v = CVXR_MDuchi_v, b = as.numeric(CVXR_MDuchi_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}




New1_kernel_pri_opt <- function(X,y,C, kernel){
  class_idx <- sort(unique(y))
  m = length(class_idx)
  n <- nrow(X)
  
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  K_p1 <- K[y!=class_idx[m],,drop = F]
  K_p2 <- K[y==class_idx[m],,drop = F]
  n_p1 <- nrow(K_p1)
  n_p2 <- nrow(K_p2)
  Y_p1 <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y[1:n_p1]==id)})
  p <- ncol(X)
  v <- Variable(n,m-1)
  b <- Variable(m-1)
  
  slack_p1 <- Variable(rows = n_p1, cols = m-1)
  slack_p2 <- Variable(rows = n_p2, cols = m-1)
  epsilon_p1 <- Variable(n_p1)
  t <- Variable(n_p1*(m-1))
  u <- Variable(rows = n_p1*(m-1), cols = m-1)
  
  objective <- Minimize(sum(do.call(rbind,sapply(1:(m-1), FUN = function(k){quad_form(v[,k],K)})))/2 +  C*(sum(epsilon_p1)+sum_entries(slack_p2)))
  constraints <- list(((K_p1%*%v + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1)%*%matrix(1, nrow = m-1, ncol = m-1) - (K_p1%*%v + matrix(1,nrow=n_p1,ncol =1)%*%t(b)) >= 1 - slack_p1,
                      vec(t(epsilon_p1%*%matrix(1,ncol=m-1))) >= t + sum_entries(rep(1/seq(1,m-1),n_p1)%*%matrix(1,ncol = m-1)*u,axis=1) - rep(1/seq(1,m-1),n_p1),
                      t%*%matrix(1,nrow=1,ncol=m-1) + u >=  (diag(1,n_p1)%x%matrix(1,nrow=m-1))%*%slack_p1,
                      u>=0,
                      sum_entries((K_p1%*%v + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1, axis = 1) >= 1 - epsilon_p1,
                      K_p2%*%v + matrix(1,nrow=n_p2,ncol =1)%*%t(b) <= 0 + slack_p2,
                      slack_p1 >=0,
                      slack_p2 >=0)
  
  New1 <- Problem(objective, constraints)
  CVXR_New1 <- solve(New1, solver = "MOSEK")
  
  CVXR_New1_v <- CVXR_New1$getValue(v)
  CVXR_New1_b <- CVXR_New1$getValue(b)
  
  model <- list(v = CVXR_New1_v, b = as.numeric(CVXR_New1_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}


New3_kernel_pri_opt <- function(X,y,C, kernel){
  
  class_idx <- sort(unique(y))
  m = length(class_idx)
  n <- nrow(X)
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  K_p1 <- K[y!=class_idx[m],,drop = F]
  K_p2 <- K[y==class_idx[m],,drop = F]
  n_p1 <- nrow(K_p1)
  n_p2 <- nrow(K_p2)
  Y_p1 <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y[1:n_p1]==id)})
  p <- ncol(X)
  v <- Variable(n,m-1)
  b <- Variable(m-1)
  
  
  slack_p1 <- Variable(rows = n_p1, cols = m-1)
  slack_p2 <- Variable(rows = n_p2, cols = m-1)
  epsilon_p1 <- Variable(n_p1)
  t <- Variable(n_p1*(m-2))
  u <- Variable(rows = n_p1*(m-2), cols = m-1)
  
  
  objective <- Minimize(sum(do.call(rbind,sapply(1:(m-1), FUN = function(k){quad_form(v[,k],K)})))/2 +  C*(sum(epsilon_p1)+sum_entries(slack_p2)))
  constraints <- list(((K_p1%*%v + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1)%*%matrix(1, nrow = m-1, ncol = m-1) - (K_p1%*%v + matrix(1,nrow=n_p1,ncol =1)%*%t(b)) >= (1-Y_p1)*(1 - slack_p1),
                      vec(t(epsilon_p1%*%matrix(1,ncol=m-2))) >=  rep(seq(1,m-2),n_p1)/rep(seq(2,m-1),n_p1)*t + 1/rep(seq(2,m-1),n_p1)*sum_entries(u,axis=1),
                      t%*%matrix(1,nrow=1,ncol=m-1) + u >=  (diag(1,n_p1)%x%matrix(1,nrow=m-2))%*%slack_p1,
                      u>=0,
                      sum_entries((K_p1%*%v + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1, axis = 1) >= 1 - epsilon_p1,
                      K_p2%*%v + matrix(1,nrow=n_p2,ncol =1)%*%t(b) <= 0 + slack_p2,
                      slack_p1 >=0,
                      slack_p2 >=0)
  
  
  New3 <- Problem(objective, constraints)
  CVXR_New3 <- solve(New3, solver = "MOSEK")
  
  CVXR_New3_v <- CVXR_New3$getValue(v)
  CVXR_New3_b <- CVXR_New3$getValue(b)
  
  model <- list(v = CVXR_New3_v, b = as.numeric(CVXR_New3_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}



LLW_kernel_pri_opt <- function(X,y,C, kernel, base_class = NULL){
  
  class_idx <- sort(unique(y))
  m <- length(class_idx)
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  n <- nrow(X)
  m <- length(unique(y))
  v <- Variable(n,m)
  b <- Variable(m)
  slack <- Variable(rows = n, cols = m)
  
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2  + C*sum_entries(slack))
  constraints <- list((K%*%v +  matrix(1,nrow=n,ncol =1)%*%t(b))*(Y-1) >= (1-Y)*(1 - slack),
                      slack >= 0,
                      sum_entries(b)*matrix(1,nrow=n) + K %*% sum_entries(v, axis = 1) == 0)
  
  LLW <- Problem(objective, constraints)
  CVXR_LLW <- solve(LLW, solver = "MOSEK")
  
  CVXR_LLW_v <- CVXR_LLW$getValue(v)
  CVXR_LLW_b <- CVXR_LLW$getValue(b)
  
  model <- list(v = CVXR_LLW_v, b = as.numeric(CVXR_LLW_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}

MSVM7_kernel_pri_opt <- function(X,y,C, kernel, base_class = NULL){
  
  class_idx <- sort(unique(y))
  m <- length(class_idx)
  if(is.null(base_class)){
    Y <- sapply(class_idx[1:(m-1)], function(id){as.numeric(y==id)})
  }else{
    Y <- sapply(class_idx[class_idx!=base_class], function(id){as.numeric(y==id)})
  }
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")

  n <- nrow(X)
  v <- Variable(n,m-1)
  b <- Variable(m-1)
  slack <- Variable(rows = n, cols = m-1)
  
  objective <- Minimize(sum(do.call(rbind,sapply(1:(m-1), FUN = function(k){quad_form(v[,k],K)})))/2  + C*sum_entries(slack))
  constraints <- list(sum_entries((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*Y, axis = 1) >= (1 - sum_entries(slack, axis = 1))*sum_entries(Y,axis=1),
                      (K%*%v +  matrix(1,nrow=n,ncol =1)%*%t(b))*(Y-1) >= (1-Y)*(0 - slack),
                      slack >= 0)
  
  MSVM7 <- Problem(objective, constraints)
  CVXR_MSVM7 <- solve(MSVM7, solver = "MOSEK")
  
  CVXR_MSVM7_v <- CVXR_MSVM7$getValue(v)
  CVXR_MSVM7_b <- CVXR_MSVM7$getValue(b)
  
  
  model <- list(v = CVXR_MSVM7_v, b = as.numeric(CVXR_MSVM7_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
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

MSVM8_kernel_pri_opt <- function(X,y,C, kernel){
  
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  v <- Variable(n,m)
  b <- Variable(m)
  
  slack <- Variable(rows = n, cols = m)
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2  + C*sum_entries(slack)   )
  constraints <- list(sum_entries((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*Y,1) >= 1 - sum_entries(slack, axis = 1),
                      (K%*%v +  matrix(1,nrow=n,ncol =1)%*%t(b))*(Y-1) >= (1-Y)*(1 - slack),
                      slack >= 0)
  
  MSVM8 <- Problem(objective, constraints)
  
  CVXR_MSVM8 <- solve(MSVM8, solver = "MOSEK")
  
  CVXR_MSVM8_v <- CVXR_MSVM8$getValue(v)
  CVXR_MSVM8_b <- CVXR_MSVM8$getValue(b)
  
  
  model <- list(v = CVXR_MSVM8_v, b = as.numeric(CVXR_MSVM8_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}



OVA_kernel_pri_opt <- function(X,y,C, kernel){
  
  class_idx <- sort(unique(y))
  Y <- sapply(class_idx, function(id){as.numeric(y==id)})
  n <- nrow(X)
  p <- ncol(X)
  m <- length(class_idx)
  K <- kernelMatrix(kernel, X)
  suppressWarnings(class(K) <- "matrix")
  
  v <- Variable(n,m)
  b <- Variable(m)
  
  slack <- Variable(rows = n, cols = m)
  objective <- Minimize(sum(do.call(rbind,sapply(1:m, FUN = function(k){quad_form(v[,k],K)})))/2  + C*sum_entries(slack)   )
  constraints <- list((K%*%v + matrix(1,nrow=n,ncol =1)%*%t(b))*(2*Y-1) >= 1 - slack,
                      slack >= 0)
  
  OVA <- Problem(objective, constraints)
  
  CVXR_OVA <- solve(OVA, solver = "MOSEK")
  
  CVXR_OVA_v <- CVXR_OVA$getValue(v)
  CVXR_OVA_b <- CVXR_OVA$getValue(b)
  
  model <- list(v = CVXR_OVA_v, b = as.numeric(CVXR_OVA_b), kernel = kernel, X = X, y= y)
  class(model) <- "msvm_kernel"
  return(model)
}


predict.msvm_kernel <- function(model, X_test,  rule = "simple_max"){
  K_test <- kernelMatrix(model$kernel, X_test, model$X)
  lpred <- K_test%*%model$v + matrix(1, nrow=nrow(K_test))%*%t(model$b)
  if(rule == "simple_max"){
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "dagger_New1"){
    lpred <- cbind(lpred, 1 - apply(as.matrix(cbind(lpred[,1],0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(lpred[,2],0)), MARGIN = 1, FUN = max))
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "dagger_MSVM7"){
    lpred <- cbind(lpred, 1 - apply(as.matrix(cbind(lpred[,1],0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(lpred[,2],0)), MARGIN = 1, FUN = max))
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "0_max"){
    lpred <- cbind(lpred, 0)
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }
  return(y_pred)
}



plot_decision_kernel_boundary <- function(model, title = NULL, np_resolution = 500, xlim =NULL, ylim = NULL, rule = "simple_max"){
  X_dat <- as.data.frame(model$X)
  colnames(X_dat) <- c("x","y")
  X_dat$class <- as.factor(model$y)  
  if(is.null(xlim)){
    xlim = c( floor(min(X_dat$x,X_dat$y))-2, ceiling(max(X_dat$x,X_dat$y))+2)
  }
  if(is.null(ylim)){
    ylim = c( floor(min(X_dat$x,X_dat$y))-2, ceiling(max(X_dat$x,X_dat$y))+2)
  }
  
  nd.x = seq(from = xlim[1], to =  xlim[2], length.out = np_resolution)
  nd.y = seq(from = ylim[1], to =  ylim[2], length.out = np_resolution)
  nd = expand.grid(x = nd.x, y = nd.y)
  
  nd$class <- predict(model, as.matrix(nd), rule = rule)
  
  colorfun <- function(n,l=65,c=100) { hues = seq(15, 375, length=n+1); hcl(h=hues, l=l, c=c)[1:n] } # default ggplot2 colours
  colors <- colorfun(length(unique(model$y)))
  colorslight <- colorfun(length(unique(model$y)),l=90,c=50)
  plt <- ggplot(nd, aes(x=x, y=y)) +
    geom_raster(data=nd, aes(x = x, y = y, fill = factor(class)),alpha=0.7,show.legend=FALSE) +
    geom_contour(data=nd, aes(x = x, y = y, z= as.numeric(nd$class)), colour="red2", alpha=0.5, breaks=c(1.5,2.5)) +
    geom_point(data = X_dat, size = 2, aes(pch = class,  colour = class)) +
    scale_x_continuous(limits = xlim, expand=c(0,0)) +
    scale_y_continuous(limits = ylim, expand=c(0,0)) +
    scale_fill_manual(values=colorslight,guide=F) + 
      ggtitle(title)
  suppressWarnings(print(plt))
} 

three_class_data_generate <- function(n, sep =  1){
  y <- sort(sample(c(1,2,3), size = n-3, replace = T, prob =  1/rep(3,3)))
  n1 <- sum(y==1)+1
  n2 <- sum(y==2)+1
  n3 <- sum(y==3)+1
  
  X1 <- matrix(mvrnorm(n1, sep*c(2,0), diag(1.5^2,nrow=2)), nrow=n1)
  c2 <- rbinom(n2, size = 1, .5)
  X2_1 <- matrix(mvrnorm(n2, sep*c(0,2), diag(1.5^2,nrow=2)), nrow=n2)
  X2_2 <- matrix(mvrnorm(n2, sep*c(0,-2), diag(1.5^2,nrow=2)), nrow=n2)
  X2 <- rbind(X2_1[c2==1,], X2_2[c2 == 0,])
  X3 <- matrix(mvrnorm(n3, sep*c(-2,0), diag(1.5^2,nrow=2)), nrow=n3)
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


four_class_data_generate <- function(n, sep =  1){
  y <- sort(sample(c(1,2,3,4), size = n-4, replace = T, prob =  1/rep(4,4)))
  n1 <- sum(y==1)+1
  n2 <- sum(y==2)+1
  n3 <- sum(y==3)+1
  n4 <- sum(y==4)+1
  # 
  # X1 <- matrix(mvrnorm(n1, sep*c(2,2), diag(1.5^2,nrow=2)), nrow=n1)
  # X2 <- matrix(mvrnorm(n2, sep*c(2,-2), diag(1.5^2,nrow=2)), nrow=n2)
  # X3 <- matrix(mvrnorm(n3, sep*c(-2,2), diag(1.5^2,nrow=2)), nrow=n3)
  # X4 <- matrix(mvrnorm(n4, sep*c(-2,-2), diag(1.5^2,nrow=2)), nrow=n4)
  ########  Plan II     #############
  X1 <- matrix(mvrnorm(n1, sep*c(-2,-2),  matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2)), nrow=n1)
  X2 <- matrix(mvrnorm(n2, sep*c(2,-2),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2)), nrow=n2)
  X3 <- matrix(mvrnorm(n3, sep*c(2,2), matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2)), nrow=n3)
  X4 <- matrix(mvrnorm(n4, sep*c(-2,2),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2)), nrow=n4)

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

five_class_data_generate <- function(n, sep =  1){
  y <- sort(sample(c(1,2,3,4,5), size = n-5, replace = T, prob =  1/rep(5,5)))
  n1 <- sum(y==1)+1
  n2 <- sum(y==2)+1
  n3 <- sum(y==3)+1
  n4 <- sum(y==4)+1
  n5 <- sum(y==5)+1
  # 
  # X1 <- matrix(mvrnorm(n1, sep*c(2,2), diag(1.5^2,nrow=2)), nrow=n1)
  # X2 <- matrix(mvrnorm(n2, sep*c(2,-2), diag(1.5^2,nrow=2)), nrow=n2)
  # X3 <- matrix(mvrnorm(n3, sep*c(-2,2), diag(1.5^2,nrow=2)), nrow=n3)
  # X4 <- matrix(mvrnorm(n4, sep*c(-2,-2), diag(1.5^2,nrow=2)), nrow=n4)
  ########  Plan II     #############
  X1 <- matrix(mvrnorm(n1, sep*c(-2,-2),  matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2)), nrow=n1)
  X2 <- matrix(mvrnorm(n2, sep*c(2,-2),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2)), nrow=n2)
  X3 <- matrix(mvrnorm(n3, sep*c(2,2), matrix(c(4.5,-3.5,-3.5,4.5), nrow= 2)), nrow=n3)
  X4 <- matrix(mvrnorm(n4, sep*c(-2,2),  matrix(c(4.5,3.5,3.5,4.5), nrow= 2)), nrow=n4)
  X5 <- matrix(mvrnorm(n5, sep*c(0,0),  matrix(c(4,0,0,4), nrow= 2)), nrow=n5)
  
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
  
  X <- rbind(X1,X2,X3,X4,X5)
  y <- c(rep(1,nrow(X1)), rep(2,nrow(X2)), rep(3,nrow(X3)), rep(4,nrow(X4)), rep(5,nrow(X5)))
  return(list(X = X,y = y))
}
  