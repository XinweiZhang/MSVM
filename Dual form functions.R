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
      b <- ginv(b_system_mat)%*%b_system_v
      return(cbind(w,b))
  } 
  return(solve_w_and_b(X,Y,alpha_sol, tau_sol))

}


MDuchi_dual_opt <- function(X,y,C){
  n <- nrow(X)
  m <- length(unique(y))
  p <- ncol(X)
  Y <- sapply(unique(y), function(id){as.numeric(y==id)})
  
  Delta <- matrix(1, nrow  = n, ncol = m) - Y
  alpha <- Variable(rows = n, cols = m)
  beta <- Variable(rows = n, cols = m-1)
  tau <- Variable(rows = n*(m-1), cols = m)
  
  objective <- Maximize(-sum_squares(t(sum_entries(alpha, axis=1)%*%matrix(1,nrow=1,ncol = m)*Y - alpha)%*%X)/2 +  sum_entries(alpha))
  constraints <- list(sum_entries(Y*(sum_entries(alpha, axis = 1)%*%matrix(1,nrow=1,ncol = m)), axis = 2) == sum_entries(alpha, axis = 2),
                      sum_entries(alpha*Y, axis = 1) ==0,
                      C - sum_entries(beta, axis =1) ==0, 
                      vec(beta*sapply(1:(m-1), function(j){rep(1/(j+1),n)}))%*%matrix(1,nrow=1,ncol=m) - tau >=0, 
                      beta*sapply(1:(m-1), function(j){rep(j/(j+1),n)}) - reshape_expr(sum_entries(do.call(rbind,lapply(1:(m-1),
                                       function(x){1-Y}))*tau, axis =1), c(n,m-1))==0,
                 
                      sum_entries(do.call(rbind,lapply(1:(m-1), function(x){Y}))*tau, axis =1) == 0,
                      alpha>=0,
                      beta>=0,
                      tau>=0) 
  constraints <- c(constraints,lapply(0:(m-1), function(id){sum_entries(reshape_expr(tau, c(n, (m-1)*m))[,(id*(m-1)+1):((id+1)*(m-1))], axis=1) - alpha[,id+1] >=0}))
  
  
  M_Duchi_dual <- Problem(objective, constraints)
  
  M_Duchi_dual <- solve(M_Duchi_dual, solver = "MOSEK")
  
  alpha_sol <- M_Duchi_dual$getValue(alpha)
  tau_sol <- M_Duchi_dual$getValue(tau)
  
  solve_w_and_b <- function(X, Y, alpha_sol, tau_sol, zero_tol=4){

    w <- t(rowSums(alpha_sol)%*%matrix(1,nrow=1,ncol = m)*Y - alpha_sol)%*%X
    lpred_dif <- X%*%t(w) - (Y*X%*%t(w))%*%matrix(1,nrow = m, ncol=m)
    m = ncol(Y)
    tau_sum <- do.call(cbind,lapply(0:(m-1), function(id,tau){
      rowSums(tau[,(id*(m-1)+1):((id+1)*(m-1))])}, tau = matrix(c(tau_sol), n, (m-1)*m, byrow = F))) 
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
    
    b <- ginv(b_system_mat)%*%b_system_v
    return(cbind(w,b))
  }
  
  return(solve_w_and_b(X,Y,alpha_sol, tau_sol))
  
}





WW_dual_opt <- function(X,y,C){
  n <- nrow(X)
  m <- length(unique(y))
  p <- ncol(X)
  
  Y <- sapply(unique(y), function(id){as.numeric(y==id)})
  
  Delta <- matrix(1, nrow  = n, ncol = m) - Y
  alpha <- Variable(rows = n, cols = m)
  beta <- Variable(rows = n, cols = m)
  
  objective <- Maximize(-sum_squares(t(sum_entries(alpha, axis=1)%*%matrix(1,nrow=1,ncol = m)*Y - alpha)%*%X)/2 +  sum_entries(alpha))
  
  
  constraints <- list(sum_entries(Y*(sum_entries(Delta*alpha, axis = 1)%*%matrix(1,nrow=1,ncol = m)), axis = 2) == sum_entries(Delta*alpha, axis = 2),
                      sum_entries(alpha*Y, axis = 1) ==0,
                      alpha>=0, 
                      alpha<=C) 
  
  
  WW_dual <- Problem(objective, constraints)
  
  WW_dual_time <- system.time(WW_dual <- solve(WW_dual, solver = "MOSEK"))
  
  alpha_sol <- WW_dual$getValue(alpha)
  
  
  solve_w_and_b <- function(X, Y, alpha_sol, C, zero_tol =4){
    w <- t(rowSums(alpha_sol)%*%matrix(1,nrow=1,ncol = m)*Y - alpha_sol)%*%X
    lpred_dif <- X%*%t(w) - (Y*X%*%t(w))%*%matrix(1,nrow = m, ncol=m)
    
    b_system_mat <- rep(1,m)
    b_system_v <- 0
    m = ncol(Y)
    for(i in 1:m){
      for( j in (i+1):m)
      {
        if(i==m){
          break
        }
        lpred_dif_b <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(alpha_sol,zero_tol) < C) * (y==i)%*%matrix(1,nrow=1, ncol = m)
        lpred_dif_b[lpred_dif_b == 0] <- NA
        b_tmp <- colMeans(lpred_dif_b, na.rm = T)[j]
        
        lpred_dif_b_inv <- (1+lpred_dif)*(round(alpha_sol,zero_tol)>0 & round(alpha_sol,zero_tol) < C) * (y==j)%*%matrix(1,nrow=1, ncol = m)
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
    
    b <- ginv(b_system_mat)%*%b_system_v
    return(cbind(w,b))
  }
  
  return(solve_w_and_b(X,Y,alpha_sol, C))
  
}





