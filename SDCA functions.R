L_MDuchi <- function(alpha, Weights){
  j <- alpha[1]
  alpha <- alpha[-1]
  alpha_sort <- sort(alpha[-j], decreasing = T)
  max(0, Weights%*%(alpha_sort - alpha[j] + 1))
}

PL_kernel_MDuchi_Loss <- function(v, b, K, y, lambda, gamma, Weights){
  if(gamma != 0){
    s <- v%*%K + replicate(length(y), b)
  }else{
    s <- v%*%K
  }
  PL_MDuchi_v <- sum(apply(cbind(y,t(s)), MARGIN = 1, FUN = L_MDuchi, Weights = Weights))/length(y)
  PL_MDuchi_v <- PL_MDuchi_v + lambda/2*sum(diag(v%*%K%*%t(v))) + gamma/2*sum(b^2)
  return(PL_MDuchi_v)
}


DL_kernel_MDuchi_Loss <- function(A, K, Y, lambda){
  DL_MDuchi_v <- -sum(colSums(A*Y))/ncol(K) - 1/(2*lambda*ncol(K)^2)*sum(diag(A%*%K%*%t(A)))
  return(DL_MDuchi_v)
}



SDCA_kernel_MDuchi <- function(X, y, intercept = T, kernel = NULL, lambda = NaN, C = NaN, gamma = 0, 
                               gap_cvg_tol = 1e-2, pri_cvg_tol = 1e-5, max_iter_num = 1000, 
                               trace = F, gam_lam_ratio = 100, type = "naive", init = "zero", K = NULL){
  
  X <- t(X)
  Y <- t(sapply(sort(unique(y)), function(id){as.numeric(y==id)}))
  m <- length(unique(y))
  n <- ncol(X)
  p <- nrow(X)
  
  if(!is.nan(C))
  {
    lambda = 1/(n*C)
  }
  if(intercept == T){
    if(gamma == 0){
      gamma = lambda/gam_lam_ratio
    }
    if(is.null(kernel)){
      K <- t(X)%*%X  + lambda/gamma
    }else{
      if(is.null(K)){
        require(kernlab)
        K <- kernelMatrix(kernel, t(X))  + lambda/gamma
      }
    }
  }else{
    if(is.null(kernel)){
      K <- t(X)%*%X  
    }else{
      if(is.null(K)){
        require(kernlab)
        K <- kernelMatrix(kernel, t(X))  
      }
    }
  }
  
  if(init == "zero"){
    v <- A <- matrix(0, m, n)
    b <- rep(0,m)
  }else{
    v <- A <- matrix(rnorm(m*n), m, n)
    b <- rnorm(m)
  }
  
  
  Weights <- replicate(m, 1/seq(1,m))
  Weights[upper.tri(Weights)] <- 0
  
  bt <- max_iter_num
  P_obj <- rep(0, max_iter_num)
  D_obj <- rep(0, length(P_obj))
  
  P_obj[1] <- (m-1)/m
  D_obj[1] <- 0
  
  for(t in 2:max_iter_num){
 
    # if(type == "naive"){
    ###naive
    perm_idx <- sample(n, replace = F)
    # }
    
    for(j in perm_idx){
      
      ##Option 6
      q <- c(A%*%K[,j] - K[j,j]*A[,j])
      q2 <- Y[,j] - (lambda*n*Y[,j] + q)/(K[j,j])
      
      od.q2 <- order(q2[-y[j]], decreasing = T)
      q2tmp <- c(q2[y[j]], q2[-y[j]][od.q2])
      values.cand <- c(Weights%*%q2tmp - 1/seq(1,m))
      
      a <- rep(1/m,m)
      for(r in 1:(m-1)){
        theta <- max(values.cand[r:m])
        tmp <- q2tmp[(r+1):m] - theta
        tmp[tmp<0] <- 0
        a_max <- (1 - sum(tmp))/r
        lam <- rep(a_max,r) - q2tmp[1:r] + theta
        if(a_max >= tmp[1] && lam[1] >=-1e-3 && all(lam[-1] <=1e-3))
        {
          a <- c(rep(a_max,r),tmp)
          break
        }
      }
      
      x_res <- a
      
      A[y[j],j] <- x_res[1] - 1
      A[-y[j],j][od.q2] <- x_res[-1]
    }
    
    v <- -A/lambda/n
    if(intercept == T){
      b <- -rowSums(A)/gamma/n
      P_obj[t] <-  PL_kernel_MDuchi_Loss(v, b, K - lambda/gamma, y, lambda = lambda, gamma = gamma, Weights = Weights[2:m,2:m])
    }else{
      P_obj[t] <- PL_kernel_MDuchi_Loss(v, b, K, y, lambda = lambda, gamma = gamma, Weights = Weights[2:m,2:m])
    }
    D_obj[t] <- DL_kernel_MDuchi_Loss(A, K, Y, lambda)
    
    if(trace == T){
      cat("P:", P_obj[t],"D:",  D_obj[t], "Gap:",  P_obj[t] - D_obj[t], "\n")
    }
    
    if(P_obj[t]-D_obj[t] < gap_cvg_tol || abs(P_obj[t]-P_obj[t-1]) < pri_cvg_tol ){
      bt <- t
      break
    }
    
  }
  
  sv_idx <- apply(v, MARGIN = 2, function(x){sum(x^2)})!=0
  fit <- list(v = t(v[,sv_idx]), b = b, kernel = kernel, SV = t(X)[sv_idx, ], sv_idx = sv_idx, nSV = sum(sv_idx), v_f = t(v), 
              X = t(X), y = y, P_obj = P_obj[1:bt], D_obj = D_obj[1:bt])
  class(fit) <- "msvm_kernel"
  return(fit)
}  





PL_kernel_CS_Loss <- function(v, b, K, Y, lambda, gamma){
  if(gamma != 0){
    s <- v%*%K + replicate(ncol(Y), b)
  }else{
    s <- v%*%K
  }
  s_y <- colSums(s*Y)
  PL_CS_V <- s - t(replicate(nrow(Y),s_y)) + (1-Y)
  PL_CS_v <- mean(apply(PL_CS_V, MARGIN = 2, FUN = max))
  PL_CS_v <- PL_CS_v + lambda/2*sum(diag(v%*%K%*%t(v))) + gamma/2*sum(b^2)
  return(PL_CS_v)
}


DL_kernel_CS_Loss <- function(A, K, Y, lambda){
  DL_CS_v <- -sum(A*Y)/ncol(K) - 1/(2*lambda*ncol(K)^2)*sum(diag(A%*%K%*%t(A)))
  return(DL_CS_v)
}

SDCA_kernel_CS <- function(X, y, intercept = T, kernel = NULL, lambda = NaN, C = NaN, gamma = 0, 
                           gap_cvg_tol = 1e-2, pri_cvg_tol = 1e-5, max_iter_num = 1000, 
                           trace = F, gam_lam_ratio = 100, type = "naive", init = "zero", K = NULL){
 
  X <- t(X)
  m <- length(unique(y))
  n <- ncol(X)
  p <- nrow(X)
  Y <- t(sapply(sort(unique(y)), function(id){as.numeric(y==id)}))
  
  if(!is.nan(C))
  {
    lambda = 1/(n*C)
  }
  if(intercept == T){
    if(gamma == 0){
      gamma = lambda/gam_lam_ratio
    }
    if(is.null(kernel)){
      K <- t(X)%*%X  + lambda/gamma
    }else{
      if(is.null(K)){
        require(kernlab)
        K <- kernelMatrix(kernel, t(X))  + lambda/gamma
      }
    }
  }else{
    if(is.null(kernel)){
      K <- t(X)%*%X  
    }else{
      if(is.null(K)){
        require(kernlab)
        K <- kernelMatrix(kernel, t(X))  
      }
    }
  }
  if(init == "zero"){
    v <- A <- matrix(0, m, n)
    b <- rep(0,m)
  }else{
    v <- A <- matrix(rnorm(m*n), m, n)
    b <- rnorm(m)
  }


  Weights <- replicate(m, 1/seq(1,m))
  Weights[upper.tri(Weights)] <- 0
  
  bt <- max_iter_num
  P_obj <- rep(0, max_iter_num)
  D_obj <- rep(0, length(P_obj))
  
  P_obj[1] <- (m-1)/m
  D_obj[1] <- 0

  for(t in 2:max_iter_num){
 
    ###greedy
    if(type == "naive"){
      ###naive
      perm_idx <- sample(n, replace = F)
    }else if(type == "greedy"){
      Fmat <- -lambda*Y - A%*%K/n
      Ftmp <- matrix(Inf, nrow = m, ncol = n)
      Ftmp[A + Y>0] <- Fmat[A + Y>0]
      psi <- apply(Fmat, MARGIN = 2, max) - apply(Ftmp, MARGIN = 2, min)
      psi_idx <- which(psi > 0)
      perm_idx <- psi_idx[sample(length(psi_idx), replace = F)]
    }else if(type == "greedy2"){
      Fmat <- -lambda*Y - A%*%K/n
      Ftmp <- matrix(Inf, nrow = m, ncol = n)
      Ftmp[A + Y>0] <- Fmat[A + Y>0]
      psi <- apply(Fmat, MARGIN = 2, max) - apply(Ftmp, MARGIN = 2, min)
      psi_idx <- which(psi > 0)
      perm_idx <- sort(psi_idx, decreasing = T)[1:round(length(psi_idx)/4)]
    }
  
   
    for(j in perm_idx){
      
      ##Option 6
      q <- c(A%*%K[,j] - K[j,j]*A[,j])
      q2 <- Y[,j] - (lambda*n*Y[,j] + q)/K[j,j]
      
      values.cand <- c(Weights%*%sort(q2, decreasing = T) - 1/seq(1,m))
      x_res <- q2 - max(values.cand)
      x_res[x_res <= 0] <- 0
      
      A[y[j],j] <- x_res[y[j]] - 1
      A[-y[j],j] <- x_res[-y[j]]
    }
    
    v <- -A/lambda/n
    if(intercept == T){
      b <- -rowSums(A)/gamma/n
      P_obj[t] <- PL_kernel_CS_Loss(v, b, K -  lambda/gamma, Y, lambda, gamma)
    }else{
      P_obj[t] <- PL_kernel_CS_Loss(v, b, K, Y, lambda, gamma)
    }
   
    D_obj[t] <- DL_kernel_CS_Loss(A, K, Y, lambda)
    if(trace == T){
      cat("P:", P_obj[t],"D:",  D_obj[t],"Gap:",  P_obj[t] - D_obj[t], "Active Set:", length(perm_idx), "\n")
    }
   
  
    if(P_obj[t]-D_obj[t] < gap_cvg_tol || abs(P_obj[t]-P_obj[t-1]) < pri_cvg_tol ){
      bt <- t
      break
    }
  }
  
  sv_idx <- apply(v, MARGIN = 2, function(x){sum(x^2)})!=0
  fit <- list(v = t(v[,sv_idx]), b = b, kernel = kernel, SV = t(X)[sv_idx, ], sv_idx = sv_idx, nSV = sum(sv_idx), v_f = t(v), 
              X = t(X), y = y, P_obj = P_obj[1:bt], D_obj = D_obj[1:bt])
  class(fit) <- "msvm_kernel"
  return(fit)
}  


PL_kernel_Duchi_Loss <- function(v, b, K, Y, lambda, gamma, Weights){
  if(gamma != 0){
    s <- v%*%K + replicate(ncol(Y), b)
  }else{
    s <- v%*%K
  }
  PL_Duchi_v <-  1 - (sum(s*Y) - sum(apply(Weights%*%apply(s, MARGIN = 2, function(x) sort(x, decreasing = T)) - replicate(ncol(Y), 1/seq(1,nrow(Y))), MARGIN = 2, max)))/ncol(Y)
  PL_Duchi_v <- PL_Duchi_v + lambda/2*sum(diag(v%*%K%*%t(v))) + gamma/2*sum(b^2)
  return(PL_Duchi_v)
}


DL_kernel_Duchi_Loss <- function(A, K, Y, lambda){
  DL_Duchi_v <- -sum(colSums(A*Y))/ncol(K) - 1/(2*lambda*ncol(K)^2)*sum(diag(A%*%K%*%t(A)))
  return(DL_Duchi_v)
}



SDCA_kernel_Duchi <- function(X, y, intercept = T, kernel = NULL, lambda = NaN, C = NaN, gamma = 0, 
                              gap_cvg_tol = 1e-2, pri_cvg_tol = 1e-5, max_iter_num = 1000, 
                              trace = F, gam_lam_ratio = 100, type = "naive", init = "zero", K = NULL){
  
  X <- t(X)
  Y <- t(sapply(sort(unique(y)), function(id){as.numeric(y==id)}))
  m <- length(unique(y))
  n <- ncol(X)
  p <- nrow(X)
  
  if(!is.nan(C))
  {
    lambda = 1/(n*C)
  }
  if(intercept == T){
    if(gamma == 0){
      gamma = lambda/gam_lam_ratio
    }
    if(is.null(kernel)){
      K <- t(X)%*%X  + lambda/gamma
    }else{
      if(is.null(K)){
        require(kernlab)
        K <- kernelMatrix(kernel, t(X))  + lambda/gamma
      }
    }
  }else{
    if(is.null(kernel)){
      K <- t(X)%*%X  
    }else{
      if(is.null(K)){
        require(kernlab)
        K <- kernelMatrix(kernel, t(X))
      }
    }
  }
  
  if(init == "zero"){
    v <- A <- matrix(0, m, n)
    b <- rep(0,m)
  }else{
    v <- A <- matrix(rnorm(m*n), m, n)
    b <- rnorm(m)
  }
  
  
  Weights <- replicate(m, 1/seq(1,m))
  Weights[upper.tri(Weights)] <- 0
  
  bt <- max_iter_num
  P_obj <- rep(0, max_iter_num)
  D_obj <- rep(0, length(P_obj))
  
  P_obj[1] <- (m-1)/m
  D_obj[1] <- 0
  
  for(t in 2:max_iter_num){
    
    # if(type == "naive"){
    ###naive
    perm_idx <- sample(n, replace = F)
    # }
    
    for(j in perm_idx){
      
      ##Option 6
      q <- c(A%*%K[,j] - K[j,j]*A[,j])
      q2 <- Y[,j] - q/K[j,j]
      rho <- lambda*n/K[j,j]
      q2.od <- order(q2, decreasing = T)
      q2tmp <- q2[q2.od]
      values.cand <- c(Weights%*%q2tmp - (1+rho)/seq(1,m))
      
      a <- rep(1/m,m)
      for(r in 1:(m-1)){
        theta <- max(values.cand[r:m])
        tmp <- q2tmp[(r+1):m] - theta
        tmp[tmp<0] <- 0
        a_max <- (1 - sum(tmp))/r
        lam <- rep(a_max,r) - q2tmp[1:r] + theta
        if(a_max >= tmp[1] &&  all(lam <=1e-3))
        {
          a <- c(rep(a_max,r),tmp)
          break
        }
      }
      x_res <- rep(0,m)
      x_res[q2.od] <- a
      x_res <-  x_res- Y[,j]
      A[,j] <- x_res
    }
    
    v <- -A/lambda/n
    if(intercept == T){
      b <- -rowSums(A)/gamma/n
      P_obj[t] <-  PL_kernel_Duchi_Loss(v, b, K - lambda/gamma, Y, lambda = lambda, gamma = gamma, Weights = Weights)
    }else{
      P_obj[t] <- PL_kernel_Duchi_Loss(v, b, K, Y, lambda = lambda, gamma = gamma, Weights = Weights)
    }
    D_obj[t] <- DL_kernel_Duchi_Loss(A, K, Y, lambda)
    
    if(trace == T){
      cat("P:", P_obj[t],"D:",  D_obj[t], "Gap:",  P_obj[t] - D_obj[t], "\n")
    }
    
    if(P_obj[t]-D_obj[t] < gap_cvg_tol || abs(P_obj[t]-P_obj[t-1]) < pri_cvg_tol ){
      bt <- t
      break
    }
    
  }
  sv_idx <- apply(v, MARGIN = 2, function(x){sum(x^2)})!=0
  fit <- list(v = t(v[,sv_idx]), b = b, kernel = kernel, SV = t(X)[sv_idx, ], sv_idx = sv_idx, nSV = sum(sv_idx), v_f = t(v), 
                X = t(X), y = y, P_obj = P_obj[1:bt], D_obj = D_obj[1:bt])
  class(fit) <- "msvm_kernel"
  return(fit)
}  



predict.msvm_kernel <- function(model, X_test,  rule = "simple_max"){
  if(is.null(model$kernel)){
    K_test <- X_test%*%t(model$SV)
  }else{
    K_test <- kernelMatrix(model$kernel, X_test, model$SV)
  }
  
  lpred <- K_test%*%model$v + t(replicate(nrow(X_test), model$b))
  if(rule == "simple_max"){
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "dagger"){
    lpred <- cbind(lpred, 1 - apply(as.matrix(cbind(lpred[,1],0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(lpred[,2],0)), MARGIN = 1, FUN = max))
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "0_max"){
    lpred <- cbind(lpred, 0)
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }
  return(y_pred)
}

predict_msvm_kernel <- function(model, X_test,  rule = "simple_max"){
  if(is.null(model$kernel)){
    K_test <- X_test%*%t(model$X)
  }else{
    K_test <- kernelMatrix(model$kernel, X_test, model$X)
  }
  
  lpred <- K_test%*%model$v_f + t(replicate(nrow(X_test), model$b))
  if(rule == "simple_max"){
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "dagger"){
    lpred <- cbind(lpred, 1 - apply(as.matrix(cbind(lpred[,1],0)), MARGIN = 1, FUN = max) - apply(as.matrix(cbind(lpred[,2],0)), MARGIN = 1, FUN = max))
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }else if(rule == "0_max"){
    lpred <- cbind(lpred, 0)
    y_pred <- apply(lpred, MARGIN = 1, FUN= which.max)
  }
  return(y_pred)
}

CV_kernel <- function(flds, X_train, y_train, cost, gamma, name, intercept = T){
  flds.acc <- sapply(1:length(flds), function(i){
    X_train_cv <- X_train[flds[[i]], , drop = F]
    y_train_cv <- y_train[flds[[i]]]
    X_test_cv <- X_train[-flds[[i]], , drop = F]
    y_test_cv <- y_train[-flds[[i]]]
    X_train_cv <- scale(X_train_cv, center = T, scale = T)
    X_test_cv <-  scale(X_test_cv, center = attr(X_train_cv,"scaled:center"), scale = attr(X_train_cv,"scaled:scale"))
    
    acc <- NA
    try({
      if(name == "OVO"){
        train_data_cv <- as.data.frame(cbind(as.factor(y_train_cv), X_train_cv))
        test_data_cv <- as.data.frame(cbind(as.factor(y_test_cv), X_test_cv))
        colnames(train_data_cv)[1] <- "y"
        colnames(test_data_cv)[1] <- "y"
        OVO_fit <- svm(as.factor(y)~., data = train_data_cv, kernel="radial", cost = cost, gamma = gamma, scale = F)
        acc <- mean(predict(OVO_fit,  newdata = test_data_cv) == y_test_cv)
        return(acc)
      }
    })
    try({
      if(name == "OVA"){
        OVA_pred <- libsvm_kernel_OVA(X = X_train_cv, y = y_train_cv, C = cost, sigma = gamma, X_test_cv, y_test_cv)
        acc <-mean(OVA_pred == y_test_cv)
        return(acc)
      }
    })
   
    try({
      if(name == "WW"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), kernel = rbfdot(sigma = gamma), type = "spoc-svc", C = cost, scaled = F)
        # fit <- WW_kernel_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept, kernel = rbfdot(sigma = gamma))
      }else if(name == "CS"){
        fit <- ksvm(x = X_train_cv, y = as.factor(y_train_cv), kernel = rbfdot(sigma = gamma), type = "kbb-svc", C = cost, scaled = F)
        # fit <- SDCA_kernel_CS(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept, kernel = rbfdot(sigma = gamma), type = "greedy")
      }else if(name == "Duchi"){
        fit <- SDCA_kernel_Duchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept, kernel = rbfdot(sigma = gamma))
      }else if(name == "MDuchi"){
        fit <- SDCA_kernel_MDuchi(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept, kernel = rbfdot(sigma = gamma))
      }else if(name == "LLW"){
        fit <- LLW_kernel_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept, kernel = rbfdot(sigma = gamma))
      }else if(name == "MSVM8"){
        fit <- MSVM8_kernel_pri_opt(X = X_train_cv, y = y_train_cv, C = cost, intercept = intercept, kernel = rbfdot(sigma = gamma))
      }
      
      pred <- predict(fit,  X_test_cv)
      acc <- mean(pred == y_test_cv)
    })
    return(acc)
  })
  return(mean(flds.acc, na.rm = T))
}



libsvm_kernel_OVA <- function(X, y, C, sigma, X_test, y_test){
  
  # v <- matrix(0, nrow = length(unique(y)), ncol = nrow(X))
  # W <- matrix(0, nrow = length(unique(y)), ncol = ncol(X))
  # b <- rep(0, length(unique(y)))
  prob <-  sapply(sort(unique(y)), function(i){
    
    X_sub1 <- X[y==i,]
    X_sub2 <- X[y!=i,]
    data <- as.data.frame(cbind(class = c(rep(1,nrow(X_sub1)), rep(-1, nrow(X_sub2))), rbind(X_sub1,X_sub2)))
    class <- rep(0,length(y_test)) 
    class[y_test==i] <- 1
    class[y_test!=i] <- -1
    data_test <- as.data.frame(cbind(class = as.factor(class), X_test))
    
    svm_fit <- svm(class~.-1, data = data, type = "C-classification",  kernel="radial", cost = C, gamma = sigma, scale = F)
    return(attr(predict(svm_fit, newdata = data_test,decision.values = TRUE),"decision.values"))
    # w = t(svm_fit$coefs) %*% kernelMatrix(rbfdot(sigma = sigma), svm_fit$SV, X)
    # a==w %*% kernelMatrix(rbfdot(sigma = sigma), X) - svm_fit$rho
    # 
    # v[i, ] <- w
    # b[i] <- -svm_fit$rho
    
  })
  
  pred <- apply(prob, MARGIN = 1, which.max)
  return(pred)
}

MSVMPack_interface <- function(type, X_train, y_train, X_test = NULL, y_test = NULL,  C, kernel, sigma, a = 0.98, tmp.dir = ""){
  
  now <- Sys.time()
  pid <- Sys.getpid()
  train_file_name <- paste0("msvm",pid,format(now, "%m%d%H%M%S"),".train")
  train_fileConn<-file(train_file_name)
  writeLines(paste(nrow(X_train),"\n",paste(ncol(X_train)), sep=""), train_fileConn)
  close(train_fileConn)
  write.table(cbind(X_train, y_train), file = train_file_name, append = T, quote = F, sep = " ", col.names = F, row.names = F)
  
  test_file_name <- paste0("msvm", pid,format(now, "%m%d%H%M%S"),".test")
  fileConn<-file(test_file_name)
  writeLines(paste(nrow(X_test),"\n",paste(ncol(X_test)), sep=""), fileConn)
  close(fileConn)
  write.table(cbind(X_test,y_test), file = test_file_name, append = T, quote = F, sep = " ", col.names = F, row.names = F)
  
  pred_file_name <- paste0("msvm", pid, format(now, "%m%d%H%M%S"),".outputs")
  
  train_file_name <- paste0(tmp.dir, train_file_name)
  test_file_name <- paste0(tmp.dir, test_file_name)
  pred_file_name <- paste0(tmp.dir, pred_file_name)
  
  
  # if(type == "WW"){
  #   svm_options = paste("-m WW -u -t 1","-c", C, "-k", kernel, "-p", sigma, "-q")
  # }else if(type == "CS"){
  #   svm_options = paste("-m CS -u -t 1","-c", C, "-k", kernel, "-p", sigma, "-q")
  # }else if(type == "LLW"){
  #   svm_options = paste("-m LLW -u -t 1","-c", C, "-k", kernel, "-p", sigma, "-q")
  # }
  if(type == "WW"){
    svm_options = paste("-m WW","-c", C, "-k", kernel, "-p", sigma, "-q -t 1 -a 0.95")
  }else if(type == "CS"){
    svm_options = paste("-m CS","-c", C, "-k", kernel, "-p", sigma, "-q -t 1 -a 0.95")
  }else if(type == "LLW"){
    svm_options = paste("-m LLW","-c", C, "-k", kernel, "-p", sigma, "-q -t 1 -a 0.95")
  }
  system(paste("trainmsvm", train_file_name, test_file_name, pred_file_name, svm_options))
  pred <- read.table(pred_file_name, header = FALSE, sep = " ")
  pred_l <- pred[,ncol(pred)]
  
  file.remove(c(train_file_name, test_file_name, pred_file_name))
  
  return(pred_l)
}


