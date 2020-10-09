
MSVM8_pri_opt1 <- function(X, y, C = 1, lambda = 1, intercept = T){
  
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
  
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(0 - slack),
                      sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >=1 - sum_entries(slack, axis=1),
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

MSVM8_pri_opt2 <- function(X, y, C = 1, lambda = 1, intercept = T){
  
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
  
  constraints <- list(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * (Y-1)) >=(1-Y)*(1/m - slack),
                      sum_entries(((X%*%t(w) + matrix(1,nrow=n,ncol =1)%*%t(b)) * Y), axis=1) >= (m-1)/m - sum_entries(slack, axis=1),
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

draw_X <- function(y){
  t <- sapply(y, function(z){
    if(z==1){
      runif(1, min = 0, max = 5)
    }else if(z==2){
      runif(1, min = 5, max = 11)
    }else if(z==3){
      runif(1, min = 11, max = 20)
    }
  })
  
  X = t(sapply(t, function(x){c(cos(x*pi/10), sin(x*pi/10))}))
  return(X)
}

generate_noise_y <- function(y, level){
  idx <- sample(length(y), size = length(y)*level, replace = F)
  y_fake <- sample(c(1,2,3), size = length(idx), replace = T)
  y[idx] <- y_fake
  return(y)
}

#
n = 300
y <- sample(c(1,2,3), size = n, replace = T)
X <- draw_X(y)

y_noise <- generate_noise_y(y, level=0.90)


# MSVM8_pri_opt(X,y_noise,C = 1, lambda = 0, intercept = T)$beta

set.seed(123)
p <- 2
m <- 3
v <- 5
X1 <- mvrnorm(10, c(-2,3), diag(v,p))
X2 <- mvrnorm(10, c(3,-2), diag(v,p))
X3 <- mvrnorm(10, c(-3,-3), diag(v,p))
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
X <- rbind(X1,X2,X3)
Y <- sapply(unique(y), function(id){as.numeric(y==id)})
  
fit1 <- MSVM8_pri_opt1(X, y, C = 2, lambda = 1, intercept = T)

fit2 <- MSVM8_pri_opt2(X, y, C = 2, lambda = 1, intercept = T)

fit1$beta
fit2$beta

plot_decision_boundary(fit1, X,y, title = "canonical")
plot_decision_boundary(fit2, X,y, title = "reparametrize")




MSVM8 <- function(gamma){
  return(c(max(1-gamma[1],  max(gamma[2],0) + max(gamma[3],0)),
           max(1-gamma[2],  max(gamma[1],0) + max(gamma[3],0)),
           max(1-gamma[3],  max(gamma[1],0) + max(gamma[2],0))))
}

lpred <- cbind(X,1)%*%t(fit1$beta)

sum(t(apply(lpred, MARGIN=1, MSVM8))*Y)
fit1$value

lpred2 <- cbind(X,1)%*%t(fit2$beta)
sum(t(apply(lpred2, MARGIN=1, MSVM8))*Y)
fit1$value



New1_pri_opt1 <- function(X, y, C = 1, lambda = 1, intercept = T, base_class = NULL){
  
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
    fit <- list(beta = cbind(CVXR_New1$getValue(w),CVXR_New1$getValue(b)), X = X, y = y)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_New1$getValue(w),0), X = X, y = y)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }
}



New1_pri_opt2 <- function(X, y, C = 1, lambda = 1, intercept = T, base_class = NULL){
  
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
                      sum_entries((X_p1%*%t(w) + matrix(1,nrow=n_p1,ncol =1)%*%t(b))*Y_p1, axis = 1) >= 1/2-epsilon_p1,
                      X_p2%*%t(w) + matrix(1,nrow=n_p2,ncol =1)%*%t(b) <= -1/2 + slack_p2,
                      slack_p1 >=0,
                      slack_p2 >=0)
  
  New1 <- Problem(objective, constraints)
  CVXR_New1 <- solve(New1, solver = "MOSEK")
  if(intercept == T){
    fit <- list(beta = cbind(CVXR_New1$getValue(w),CVXR_New1$getValue(b)), X = X, y = y)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }else{  
    fit <- list(beta = cbind(CVXR_New1$getValue(w),0), X = X, y = y)
    class(fit) <- "msvm"
    attr(fit,"base_class") <- base_class
    return(fit)
  }
}


set.seed(123)
p <- 2
m <- 3
v <- 5
X1 <- mvrnorm(10, c(-2,3), diag(v,p))
X2 <- mvrnorm(10, c(3,-2), diag(v,p))
X3 <- mvrnorm(10, c(-3,-3), diag(v,p))
y <- c(rep(1,nrow(X1)),rep(2,nrow(X2)),rep(3,nrow(X3)))
X <- rbind(X1,X2,X3)
Y <- sapply(unique(y), function(id){as.numeric(y==id)})

fit1 <- New1_pri_opt1(X, y, C = 2, lambda = 1, intercept = T)

fit2 <- New1_pri_opt2(X, y, C = 2, lambda = 1, intercept = T)

fit1$beta
fit2$beta

plot_decision_boundary(fit1,rule="dagger", print = T)
# debugonce(plot_decision_boundary)
plot_decision_boundary(fit2,rule="dagger2", print = T)
