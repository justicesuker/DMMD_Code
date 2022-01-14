# Update function for R
updateR <- function(X, P_rall, P_c, dimR){
  n = nrow(X)
  temp_init = (diag(n) - P_c) %*% X %*% P_rall
  R = svd(temp_init)$u[, 1:dimR, drop = F]
  return(R)
}

# Update function for S
updateS <- function(X, P_call, P_r, dimS){
  p = ncol(X)
  temp_init = P_call %*% X %*% (diag(p) - P_r)
  S = svd(temp_init)$v[, 1:dimS, drop = F]
  return(S)
}

# Update function for M
updateM <- function(X1, X2, PR1, PR2, P_rall1, P_rall2, rc){
  n = nrow(X1)
  M1 = (diag(n) - PR1) %*% X1 %*% P_rall1
  M2 = (diag(n) - PR2) %*% X2 %*% P_rall2
  M = svd(cbind(M1, M2))$u[, 1:rc, drop = F]
  return(M)
}

# Update function for N
updateN <- function(X1, X2, PS1, PS2, P_call1, P_call2, rr){
  p = ncol(X1)
  N1 =  P_call1 %*% X1 %*% (diag(p) - PS1)
  N2 =  P_call2 %*% X2 %*% (diag(p) - PS2)
  N = svd(rbind(N1, N2))$v[, 1:rr, drop = F]
  return(N)
}

#' Main function of iterative DMMD algorithm
#'
#' @param X1 The first matrix.
#' @param x2 The second matrix.
#' @param eps Tolerence, default is the square root of machine precision.
#' @param r1 The total rank X1. Default is NULL, which means unknown.
#' @param r2 The total rank X2. Default is NULL, which means unknown.
#' @param rc The joint column rank. Default is NULL, which means unknown.
#' @param rr The joint row rank. Default is NULL, which means unknown.
#' @param kmax The maximum iterations. Default is 1000.
#' @param verbose Do you want to see the progress of the function? Default is F, which means there is no progress shown.
DMMD_i <- function(X1, X2, r1 = NULL, r2 = NULL, rc = NULL, rr = NULL, eps = .Machine$double.eps^0.5, kmax = 1000, verbose = FALSE){
  n = nrow(X1)
  p = ncol(X1)
  # Check if the specified ranks are legal
  if (!is.null(r1) | !is.null(r2)){
    if (max(r1,r2) > min(n, p)){
      stop("The specified rank is not legal, please check.")
    }
  }
  # Check if the specified joint rank is legal
  if (!is.null(rc) | !is.null(rr)){
    if (max(rc, rr) > min(r1, r2)){
      stop("The specified joint rank is not legal, please check.")
    }
  }
  # Save the svd result of the original matrices
  svd_x1 = svd(X1)
  svd_x2 = svd(X2)

  # Get the estimated total rank of X1 and X2. Store it as r1 and r2.
  if (is.null(r1)){
    r1 = ProfileLikCluster(svd_x1$d, variance = 'equal')$index
  }
  if (is.null(r2)){
    r2 = ProfileLikCluster(svd_x2$d, variance = 'equal')$index
  }
  # Get the estimated column/row space of X1 and X2
  X1_est_c = as.matrix(svd_x1$u[,1:r1])
  X2_est_c = as.matrix(svd_x2$u[,1:r2])
  X1_est_r = as.matrix(svd_x1$v[,1:r1])
  X2_est_r = as.matrix(svd_x2$v[,1:r2])
  
  # Calculate joint column space
  # Get the principal angles
  angle_result_c = angle_cal(X1_est_c, X2_est_c, tol = tol)
  principal_angle_c = angle_result_c$angle
  # Get the principal vectors
  pv1_c = angle_result_c$principal_vector1
  pv2_c = angle_result_c$principal_vector2
  # If the specified joint column rank is NULL. Calculate it using the PL method specified.
  if (is.null(rc)){
    rc = joint_angle_cluster(principal_angle_c, variance = 'equal')$joint_rank
  }
  # Calculate joint row space
  angle_result_r = angle_cal(X1_est_r, X2_est_r, tol = tol)
  # Get the principal angles
  principal_angle_r = angle_result_r$angle
  # Get the principal vectors
  pv1_r = angle_result_r$principal_vector1
  pv2_r = angle_result_r$principal_vector2
  # If the specified joint row rank is NULL. Calculate it using the PL method specified.
  if (is.null(rr)){
    rr = joint_angle_cluster(principal_angle_r,variance = 'equal')$joint_rank
  }
  # Get initial estimates for M and N by averaging
  # Calculate joint column space projection matrix (this is MM')
  if (rc == 0){
    joint_space_c = matrix(rep(0, n), nrow = n)
    P_c = projection(joint_space_c, ortho = TRUE)
  }
  else{
    joint_space_c = (pv1_c[,1:rc] + pv2_c[,1:rc])/2
    P_c = projection(joint_space_c)
  }
  # Calculate joint row space projection matrix (this is NN')
  if (rr == 0){
    joint_space_r = matrix(rep(0, p), nrow = p)
    P_r = projection(joint_space_r, ortho = TRUE)
  }
  else{
    joint_space_r = (pv1_r[,1:rr] + pv2_r[,1:rr])/2
    P_r = projection(joint_space_r)
  }
  
  # Initialize R1, R2 and complete full column space
  if (r1 == rc){
    R1 = updateR(X1, P_rall = diag(p), P_c, dimR = r1 - rc)
    P_call1 = projection(R1, ortho = TRUE)
  }
  else{
    R1 = updateR(X1, P_rall = diag(p), P_c, dimR = r1 - rc)
    P_call1 = P_c + projection(R1, ortho = TRUE)
  }
  if (r2 == rc){
    R2 = updateR(X2, P_rall = diag(p), P_c, dimR = r2 - rc)
    P_call2 = projection(R2, ortho = TRUE)
  }
  else{
    R2 = updateR(X2, P_rall = diag(p), P_c, dimR = r2 - rc)
    P_call2 = P_c + projection(R2, ortho = TRUE)
  }
  
  # Initialize S1, S2 and complete full column space
  if (r1 == rr){
    S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
    P_rall1 = projection(S1, ortho = TRUE)
  }
  else{
    S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
    P_rall1 = P_r + projection(S1, ortho = TRUE)
  }
  if (r2 == rr){
    S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
    P_rall2 = projection(S2, ortho = TRUE)
  }
  else{
    S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
    P_rall2 = P_r + projection(S2, ortho = TRUE)
  }
  
  # Calculate estimated signals and objective value
  A1 = P_call1 %*% X1 %*% P_rall1
  A2 = P_call2 %*% X2 %*% P_rall2
  obj_old = sum((X1-A1)^2) + sum((X2-A2)^2)
  
  obj_vec = c(obj_old) 
  k = 0
  error = 1000
  
  # while loop
  while((error > eps)&(k < kmax)){
    # Not MONOTONE right now -> double check what is going on
    # Count iterations
    k = k + 1
    if (verbose){
      print(obj_old)
      print("Update M")
    }
    
    # update M
    if (rc == 0){
      M = matrix(rep(0, n), nrow = n)
      P_c = projection(M, ortho = TRUE)
      P_call1 = projection(R1, ortho = TRUE)
      P_call2 = projection(R2, ortho = TRUE)
    }
    else{
      M = updateM(X1, X2, PR1 = projection(R1, ortho = TRUE), PR2 = projection(R2, ortho = TRUE), P_rall1, P_rall2, rc)
      P_c = projection(M, ortho = TRUE)
      P_call1 = P_c + projection(R1, ortho = TRUE)
      P_call2 = P_c + projection(R2, ortho = TRUE)
    }
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    if (verbose){
      print(obj_new)
      print("Update N")
    }
    # update N
    if (rr == 0){
      N = matrix(rep(0, p), nrow = p)
      P_r = projection(N, ortho = TRUE)
      P_rall1 = projection(S1, ortho = TRUE)
      P_rall2 = projection(S2, ortho = TRUE)
    }
    else{
      N = updateN(X1, X2, PS1 = projection(S1, ortho = TRUE), PS2 = projection(S2, ortho = TRUE), P_call1, P_call2, rr)
      P_r = projection(N, ortho = TRUE)
      P_rall1 = P_r + projection(S1, ortho = TRUE)
      P_rall2 = P_r + projection(S2, ortho = TRUE)
    }
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    
    error2 = 1000
    while(error2 > 1e-2){
      objMN = obj_new
      if (verbose){
        print(obj_new)
        print("Update R1, R2")
      }
      # update R and full column space
      if (r1 == rc){
        R1 = matrix(rep(0, n), nrow = n)
        P_call1 = P_c
      }
      else{
        R1 = updateR(X1, P_rall = P_rall1, P_c, dimR = r1 - rc)
        P_call1 = P_c + projection(R1, ortho = TRUE)
      }
      if (r2 == rc){
        R2 = matrix(rep(0, n), nrow = n)
        P_call2 = P_c
      }
      else{
        R2 = updateR(X2, P_rall = P_rall2, P_c, dimR = r2 - rc)
        P_call2 = P_c + projection(R2, ortho = TRUE)
      }
      A1 = P_call1 %*% X1 %*% P_rall1
      A2 = P_call2 %*% X2 %*% P_rall2
      obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
      
      if (verbose){
        print(obj_new)
        print("Update S1, S2")
      }
      # update S and full row space
      if (r1 == rr){
        S1 = matrix(rep(0, p), nrow = p)
        P_rall1 = P_r 
      }
      else{
        S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
        P_rall1 = P_r + projection(S1, ortho = TRUE)
      }
      if (r2 == rr){
        S2 = matrix(rep(0, p), nrow = p)
        P_rall2 = P_r
      }
      else{
        S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
        P_rall2 = P_r + projection(S2, ortho = TRUE)
      }
      # calculate objective function difference
      A1 = P_call1 %*% X1 %*% P_rall1
      A2 = P_call2 %*% X2 %*% P_rall2
      # obj_new = sum((X1-A1)^2)/X1F2 + sum((X2-A2)^2)/X2F2
      obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
      error2 = abs(obj_new - objMN)
    }
    error = abs(obj_new - obj_old)
    if (verbose){
      print(obj_new)
    }
    obj_old = obj_new
    obj_vec = c(obj_vec,obj_new)
  }
  
  return(list(A1 = A1, A2 = A2, M = M, N = N, r1 = r1, r2 = r2, rc = rc, rr = rr, obj_vec = obj_vec))
}

# Another implementation of iterative DMMD when given ranks. It gives the same results as DMMD_i. 
DMMD_All <- function(X1, X2, r1, r2, rc, rr, eps = .Machine$double.eps^0.5, kmax = 1000){
  n = nrow(X1)
  p = ncol(X1)
  X1F2 = sum(X1^2)
  X2F2 = sum(X2^2)
  
  # Save the svd result of the original matrices
  svd_x1 = svd(X1)
  svd_x2 = svd(X2)
  
  # Get initial estimates for M and N by averaging
  # Calculate joint column space projection matrix (this is MM')
  if (rc == 0){
    joint_space_c = matrix(rep(0, n), nrow = n)
    P_c = projection(joint_space_c, ortho = TRUE)
  }
  else{
    angle_result_c = angle_cal(svd_x1$u[,1:r1], svd_x2$u[,1:r2], tol = tol)
    pv1_c = angle_result_c$principal_vector1
    pv2_c = angle_result_c$principal_vector2
    joint_space_c = (pv1_c[,1:rc] + pv2_c[,1:rc])/2
    P_c = projection(joint_space_c, ortho = TRUE)
  }
  # Calculate joint row space projection matrix (this is NN')
  if (rr == 0){
    joint_space_r = matrix(rep(0, p), nrow = p)
    P_r = projection(joint_space_r, ortho = TRUE)
  }
  else{
    angle_result_r = angle_cal(svd_x1$v[,1:r1], svd_x2$v[,1:r2], tol = tol)
    pv1_r = angle_result_r$principal_vector1
    pv2_r = angle_result_r$principal_vector2
    joint_space_r = (pv1_r[,1:rr] + pv2_r[,1:rr])/2
    P_r = projection(joint_space_r, ortho = TRUE)
  }
  
  # Initialize R1, R2 and complete full column space
  if (r1 == rc){
    R1 = updateR(X1, P_rall = diag(p), P_c, dimR = r1 - rc)
    P_call1 = projection(R1, ortho = TRUE)
  }
  else{
    R1 = updateR(X1, P_rall = diag(p), P_c, dimR = r1 - rc)
    P_call1 = P_c + projection(R1, ortho = TRUE)
  }
  if (r2 == rc){
    R2 = updateR(X2, P_rall = diag(p), P_c, dimR = r2 - rc)
    P_call2 = projection(R2, ortho = TRUE)
  }
  else{
    R2 = updateR(X2, P_rall = diag(p), P_c, dimR = r2 - rc)
    P_call2 = P_c + projection(R2, ortho = TRUE)
  }
  
  # Initialize S1, S2 and complete full column space
  if (r1 == rr){
    S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
    P_rall1 = projection(S1, ortho = TRUE)
  }
  else{
    S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
    P_rall1 = P_r + projection(S1, ortho = TRUE)
  }
  if (r2 == rr){
    S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
    P_rall2 = projection(S2, ortho = TRUE)
  }
  else{
    S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
    P_rall2 = P_r + projection(S2, ortho = TRUE)
  }
  
  # Calculate estimated signals and objective value
  A1 = P_call1 %*% X1 %*% P_rall1
  A2 = P_call2 %*% X2 %*% P_rall2
  obj_old = sum((X1-A1)^2) + sum((X2-A2)^2)
  
  k = 0
  error = 1000
  obj_vec = c(obj_old)
  # while loop
  while((error > eps)&(k < kmax)){
    # Count iterations
    k = k + 1
    print("Update M")
    # update M
    if (rc == 0){
      M = matrix(rep(0, n), nrow = n)
      P_c = projection(M, ortho = TRUE)
      P_call1 = projection(R1, ortho = TRUE)
      P_call2 = projection(R2, ortho = TRUE)
    }
    else{
      M = updateM(X1, X2, PR1 = projection(R1, ortho = TRUE), PR2 = projection(R2, ortho = TRUE), P_rall1, P_rall2, rc)
      P_c = projection(M, ortho = TRUE)
      P_call1 = P_c + projection(R1, ortho = TRUE)
      P_call2 = P_c + projection(R2, ortho = TRUE)
    }
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    print(obj_new)
    
    print("Update N")
    # update N
    if (rr == 0){
      N = matrix(rep(0, p), nrow = p)
      P_r = projection(N, ortho = TRUE)
      P_rall1 = projection(S1, ortho = TRUE)
      P_rall2 = projection(S2, ortho = TRUE)
    }
    else{
      N = updateN(X1, X2, PS1 = projection(S1, ortho = TRUE), PS2 = projection(S2, ortho = TRUE), P_call1, P_call2, rr)
      P_r = projection(N, ortho = TRUE)
      P_rall1 = P_r + projection(S1, ortho = TRUE)
      P_rall2 = P_r + projection(S2, ortho = TRUE)
    }
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    print(obj_new)
    
    print("Update R1, R2")
    # update R and full column space
    if (r1 == rc){
      R1 = matrix(rep(0, n), nrow = n)
      P_call1 = P_c
    }
    else{
      R1 = updateR(X1, P_rall = P_rall1, P_c, dimR = r1 - rc)
      P_call1 = P_c + projection(R1, ortho = TRUE)
    }
    if (r2 == rc){
      R2 = matrix(rep(0, n), nrow = n)
      P_call2 = P_c
    }
    else{
      R2 = updateR(X2, P_rall = P_rall2, P_c, dimR = r2 - rc)
      P_call2 = P_c + projection(R2, ortho = TRUE)
    }
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    print(obj_new)
    
    print("Update S1, S2")
    # update S and full row space
    if (r1 == rr){
      S1 = matrix(rep(0, p), nrow = p)
      P_rall1 = P_r
    }
    else{
      S1 = updateS(X1, P_call = P_call1, P_r = P_r, dimS = r1 - rr)
      P_rall1 = P_r + projection(S1, ortho = TRUE)
    }
    if (r2 == rr){
      S2 = matrix(rep(0, p), nrow = p)
      P_rall2 = P_r
    }
    else{
      S2 = updateS(X2, P_call = P_call2, P_r = P_r, dimS = r2 - rr)
      P_rall2 = P_r + projection(S2, ortho = TRUE)
    }
    # calculate objective function difference
    A1 = P_call1 %*% X1 %*% P_rall1
    A2 = P_call2 %*% X2 %*% P_rall2
    obj_new = sum((X1-A1)^2) + sum((X2-A2)^2)
    print(obj_new)
    error = abs(obj_new - obj_old)
    obj_old = obj_new
    obj_vec = c(obj_vec,obj_new)
  }
  
  return(list(A1 = A1, A2 = A2, M = M, N = N, obj_vec = obj_vec))
}