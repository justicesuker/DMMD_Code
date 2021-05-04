# The main function that does double-matched matrix decomposition.
# X1 and X2 are two double-matched matrices
# r1, r2 are the total rank of r1 and r2. Default is NULL. If the total rank is known by some ways, please specify them.
# joint_rank_c, joint_rank_r are the specified joint column rank and joint row rank. Default is NULL. It they are known, please specify.
# angle_threshold: the argument that is used in calculating the joint rank, Principal angles that is greater than the threshold is not considered as joint signal. Default is 90 degree.
# variance1. Either "equal" or "unequal". Default is "unequal". This argument is used in the profile likelihood method for determining the total rank
# variance2. Either "equal" or "unequal". Default is "unequal". This argument is used in the profile likelihood method for determining the joint rank
# throw. Either F or T. Default is F. This argument is used in the profile likelihood method for determining both total rank and joint rank
# method. Either "PL" (profile likelihood) or "ED" (edge distribution). Default is "PL" for determining the rank.
# tol. Tolerence for determining convergence.
# maxiter. Default is 1000, which is used for the maximum iteration allowed in the iterative procedure.
DMMD_v2 <- function(X1, X2, r1 = NULL, r2 = NULL, joint_rank_c = NULL, joint_rank_r = NULL, angle_threshold = 90 * pi/180, variance1 = "equal", variance2 = "equal", throw = FALSE, method = "PL", tol = .Machine$double.eps^0.5, maxiter = 1e3){
  # Check the input of method
  if (method != "PL" & method != "ED"){
    stop("Method must be either 'ED' or 'PL'.")
  }
  # Check if the column names are equal
  if (!identical(colnames(X1), colnames(X2))){
    warning("This is an algorithm for double matched matrices. The column names of given matrices do not match")
  }
  # Check if the row names are equal
  if (!identical(rownames(X1), rownames(X2))){
    warning("This is an algorithm for double matched matrices. The row names of given matrices do not match")
  }
  n = dim(X1)[1]
  p = dim(X1)[2]
  # Check if the specified ranks are legal
  if (!is.null(r1) | !is.null(r2)){
    if (max(r1,r2) > min(n, p)){
      stop("The specified rank is not legal, please check.")
    }
  }
  # Save the svd result of the original matrices
  svd_x1 = svd(X1)
  svd_x2 = svd(X2)
  # Get the estimated total rank of X1 and X2. Store it as r1 and r2.
  if (is.null(r1)){
    if (method == "PL"){
      r1 = ProfileLikCluster(svd_x1$d, variance = variance1)$index
    }
    if (method == "ED"){
      r1 = Select_ED_Rank(svd_x1$d, maxiter = maxiter)
    }
  }
  if (is.null(r2)){
    if (method == "PL"){
      r2 = ProfileLikCluster(svd_x2$d, variance = variance1)$index
    }
    if (method == "ED"){
      r2 = Select_ED_Rank(svd_x2$d, maxiter = maxiter)
    }
  }
  # Check if the specified joint rank is legal
  if (!is.null(joint_rank_c) | !is.null(joint_rank_r)){
    if (max(joint_rank_c,joint_rank_r) > min(r1, r2)){
      stop("The specified joint rank is not legal, please check.")
    }
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
  # If the specified joint column rank is NULL. Calculate it using the PL or ED method specified.
  if (is.null(joint_rank_c)){
    joint_rank_c = joint_angle_cluster(
      principal_angle_c, angle_threshold = angle_threshold, variance = variance2, throw = throw, maxiter = maxiter)$joint_rank
  }
  # Get the estimated column space by averaging the smallest joint_rank_c number of principal vectors  
  if (joint_rank_c > 0){
    joint_space_c = (pv1_c[,1:joint_rank_c] + pv2_c[,1:joint_rank_c])/2
    P_c = projection(joint_space_c)
  } 

  # Calculate joint row space
  angle_result_r = angle_cal(X1_est_r, X2_est_r, tol = tol)
  # Get the principal angles
  principal_angle_r = angle_result_r$angle
  # Get the principal vectors
  pv1_r = angle_result_r$principal_vector1
  pv2_r = angle_result_r$principal_vector2
  # If the specified joint row rank is NULL. Calculate it using the PL or ED method specified.
  if (is.null(joint_rank_r)){
    joint_rank_r = joint_angle_cluster(
      principal_angle_r, angle_threshold = angle_threshold, variance = variance2, throw = throw, maxiter = maxiter)$joint_rank
  }
  # Get the estimated row space by averaging the smallest joint_rank_r number of principal vectors 
  if (joint_rank_r > 0){
    joint_space_r = (pv1_r[,1:joint_rank_r] + pv2_r[,1:joint_rank_r])/2
    P_r = projection(joint_space_r)
  }

  # Consider the extreme cases when joint rank is 0
  if (joint_rank_c == 0 || joint_rank_r == 0){
    if (joint_rank_c == 0){
      # Both joint column and row rank are 0: joint structure is 0.
      if (joint_rank_r == 0){
        signal_mat1 = svd_recover(X_1, svd_result = svd_x1, r1)
        signal_mat2 = svd_recover(X_2, svd_result = svd_x2, r2)
        J1_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        J2_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        J1_r = matrix(rep(0,n*p), nrow = n, ncol = p)
        J2_r = matrix(rep(0,n*p), nrow = n, ncol = p)
        I1_c = signal_mat1 
        I2_c = signal_mat2
        I1_r = signal_mat1
        I2_r = signal_mat2
        E1 = X1 - signal_mat1
        E2 = X2 - signal_mat2
      }
      # Joint column rank is 0. Joint row rank > 0 
      else{
        # # Added 2020/9/4, this condition takes care of the situation when one of the individual strucuture is 0
        # # Then we simply use the estimated signal to be the joint structure.
        # if (joint_rank_r == min(r1,r2)){
        #   if (joint_rank_r == r1){
        #     # The 1st individual row structure is 0, then the joint row structure is simply svd estimate of X1
        #     joint_space_r = X1_est_r
        #     signal_mat1 = svd_recover(X_1, svd_result = svd_x1, r1)
        #     J1_r = signal_mat1
        #     J1_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        #     I1_r = matrix(rep(0,n*p), nrow = n, ncol = p)
        #     I1_c = signal_mat1
        #     signal_mat2 = t(FindOpt_SM(t(X2), joint_space_r, r2))
        #     J2_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        #     I2_c = signal_mat2
        #     J2_r = signal_mat2 %*% projection(joint_space_r)
        #     I2_r = signal_mat2 - J2_r
        #     E1 = X1 - signal_mat1
        #     E2 = X2 - signal_mat2
        #   }
        #   else{
        #     joint_space_r = X2_est_r
        #     signal_mat2 = svd_recover(X_2, svd_result = svd_x2, r2)
        #     J2_r = signal_mat2
        #     J2_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        #     I2_r = matrix(rep(0,n*p), nrow = n, ncol = p)
        #     I2_c = signal_mat2
        #     
        #     signal_mat1 = t(FindOpt_SM(t(X1), joint_space_r, r1))
        #     J1_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        #     I1_c = signal_mat1
        #     J1_r = signal_mat1 %*% projection(joint_space_r)
        #     I1_r = signal_mat1 - J1_r
        #     E1 = X1 - signal_mat1
        #     E2 = X2 - signal_mat2
        #   }
        # }
        # else{
        
        # }
        
        # Use the function of simplified version with the joint row space to find out the solution of signal matrices. 
        signal_mat1 = t(FindOpt_SM(t(X1), joint_space_r, r1))
        signal_mat2 = t(FindOpt_SM(t(X2), joint_space_r, r2))
        J1_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        J2_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        J1_r = signal_mat1 %*% P_r
        J2_r = signal_mat2 %*% P_r
        I1_c = signal_mat1
        I2_c = signal_mat2
        I1_r = signal_mat1 - J1_r
        I2_r = signal_mat2 - J2_r 
        E1 = X1 - signal_mat1
        E2 = X2 - signal_mat2
      }
    }
    # Joint column rank is NOT 0, joint row rank is 0
    else{
      # # Added 2020/9/4, this condition takes care of the situation when one of the individual column strucuture is 0
      # # Then we simply use the estimated signal to be the joint structure.
      # if (joint_rank_c == min(r1,r2)){
      #   if (joint_rank_c == r1){
      #     # The 1st individual column structure is 0, then the joint column structure is simply svd estimate of X1
      #     joint_space_c = X1_est_c
      #     signal_mat1 = svd_recover(X_1, svd_result = svd_x1, r1)
      #     J1_c = signal_mat1
      #     J1_r = matrix(rep(0,n*p), nrow = n, ncol = p)
      #     I1_c = matrix(rep(0,n*p), nrow = n, ncol = p)
      #     I1_r = signal_mat1
      #     
      #     signal_mat2 = FindOpt_SM(X2, joint_space_c, r2)
      #     J2_r = matrix(rep(0,n*p), nrow = n, ncol = p)
      #     I2_r = signal_mat2
      #     J2_c = projection(joint_space_c) %*% signal_mat2 
      #     I2_c = signal_mat2 - J2_c
      #     E1 = X1 - signal_mat1
      #     E2 = X2 - signal_mat2
      #   }
      #   else{
      #     joint_space_c = X2_est_c
      #     signal_mat2 = svd_recover(X_2, svd_result = svd_x2, r2)
      #     J2_c = signal_mat2
      #     J2_r = matrix(rep(0,n*p), nrow = n, ncol = p)
      #     I2_c = matrix(rep(0,n*p), nrow = n, ncol = p)
      #     I2_r = signal_mat2
      #     
      #     signal_mat1 = FindOpt_SM(X1, joint_space_c, r1)
      #     J1_r = matrix(rep(0,n*p), nrow = n, ncol = p)
      #     I1_r = signal_mat1
      #     J1_c = projection(joint_space_c) %*% signal_mat1 
      #     I1_c = signal_mat1 - J1_c
      #     E1 = X1 - signal_mat1
      #     E2 = X2 - signal_mat2
      #   }
      # }
      # else{
      
      # }
      
      # Use the function of simplified version with the joint column space to find out the solution of signal matrices. 
      signal_mat1 = FindOpt_SM(X1, joint_space_c, r1)
      signal_mat2 = FindOpt_SM(X2, joint_space_c, r2)
        
      J1_c = P_c %*% signal_mat1
      J2_c = P_c %*% signal_mat2
      J1_r = matrix(rep(0,n*p), nrow = n, ncol = p)
      J2_r = matrix(rep(0,n*p), nrow = n, ncol = p)
        
      I1_c = signal_mat1 - J1_c
      I2_c = signal_mat2 - J2_c
      I1_r = signal_mat1 
      I2_r = signal_mat2 
      E1 = X1 - signal_mat1
      E2 = X2 - signal_mat2
    }
  }
  else{
    # For general cases use iterative algorithm to solve the optimization problem
    result1 = FindOpt_DM_Iterative(X1, joint_space_c, joint_space_r, r1, maxiter = maxiter, tol = tol)
    result2 = FindOpt_DM_Iterative(X2, joint_space_c, joint_space_r, r2, maxiter = maxiter, tol = tol)
    signal_mat1 = result1$result
    signal_mat2 = result2$result
    
    # Get the decomposition
    J1_c = P_c %*% signal_mat1
    J2_c = P_c %*% signal_mat2
    J1_r = signal_mat1 %*% P_r
    J2_r = signal_mat2 %*% P_r
    I1_c = signal_mat1 - J1_c
    I2_c = signal_mat2 - J2_c
    I1_r = signal_mat1 - J1_r
    I2_r = signal_mat2 - J2_r 
    E1 = X1 - signal_mat1
    E2 = X2 - signal_mat2
  }
  
  # Make sure that the rownames as well as the colnames stay the same
  names_r = rownames(X1)
  names_c = colnames(X1)
  rownames(J1_c) = rownames(J2_c) = rownames(J1_r) = rownames(J2_r) = names_r
  rownames(I1_c) = rownames(I2_c) = rownames(I1_r) = rownames(I2_r) = names_r
  rownames(E1) = rownames(E2) = names_r 
  
  colnames(J1_c) = colnames(J2_c) = colnames(J1_r) = colnames(J2_r) = names_c
  colnames(I1_c) = colnames(I2_c) = colnames(I1_r) = colnames(I2_r) = names_c
  colnames(E1) = colnames(E2) = names_c
  # Prepare for output
  rank_information = list("Rank 1" = r1, "Rank 2" = r2, "Joint Column Rank" = joint_rank_c, "Joint Row Rank" = joint_rank_r)
  column_decomposition = list("Joint Column 1" = J1_c, "Individual Column 1" = I1_c, 
                              "Joint Column 2" = J2_c, "Individual Column 2" = I2_c)
  row_decomposition = list("Joint Row 1" = J1_r, "Individual Row 1" = I1_r, 
                           "Joint Row 2" = J2_r, "Individual Row 2" = I2_r)
  error = list("Error1" = E1, "Error2" = E2)
  
  return(list("Rank" = rank_information, "Column Decomposition" = column_decomposition, "Row Decomposition" = row_decomposition, "Error" = error))
}

# This function finds the decomposition based on the input of two data matrices and possible given individual structure
# If Ind is not given, the function is the same as DMMD_v2
# X1 and X2 are two double-matched matrices
# r1, r2 are the total rank of r1 and r2. Default is NULL. If the total rank is known by some ways, please specify them.
# joint_rank_c, joint_rank_r are the specified joint column rank and joint row rank. Default is NULL. It they are known, please specify.
# angle_threshold: the argument that is used in calculating the joint rank, Principal angles that is greater than the threshold is not considered as joint signal. Default is 90 degree.
# variance1. Either "equal" or "unequal". Default is "equal". This argument is used in the profile likelihood method for determining the total rank
# variance2. Either "equal" or "unequal". Default is "equal". This argument is used in the profile likelihood method for determining the joint rank
# throw. Either F or T. Default is F. This argument is used in the profile likelihood method for determining both total rank and joint rank
# method. Either "PL" (profile likelihood) or "ED" (edge distribution). Default is "PL" for determining the rank.
# tol. Tolerence for determining convergence.
# maxiter. Default is 1000, which is used for the maximum iteration allowed in the iterative procedure.
# Ind is the list of input for given individual structute. 
# I1c = Ind[[1]]: the individual column sturcture for the 1st matrix 
# I2c = Ind[[2]]: the individual column sturcture for the 2nd matrix 
# I1r = Ind[[3]]: the individual row sturcture for the 1st matrix 
# I2r = Ind[[4]]: the individual row sturcture for the 2nd matrix 
# Notice that if individual structure is given, you need also input the joint column rank and joint row rank.
Find_Decom <- function(X1, X2, Ind = NULL, r1 = NULL, r2 = NULL, joint_rank_c = NULL, joint_rank_r = NULL, angle_threshold = 90 * pi/180, 
                       variance1 = "equal", variance2 = "equal", throw = FALSE, method = "PL", tol = .Machine$double.eps^0.5, maxiter = 1e3){
  # Get the dimension of the signal matrices
  # Check input
  if (!is.null(Ind) & (is.null(joint_rank_c) | is.null(joint_rank_r))){
    stop("If individual structure is known, you need also input the joint column rank and joint row rank.")
  }
  n = dim(X1)[1]
  p = dim(X1)[2]
  zero_matrix = matrix(rep(0,n*p),nrow = n)
  # If there is no input for Ind, initialize them to be all zero.
  # Save whether Ind is null, will be used later on 
  isnull_Ind = is.null(Ind)
  if (isnull_Ind){
    Ind = list(zero_matrix,zero_matrix,zero_matrix,zero_matrix)
  }
  I1c = Ind[[1]]
  I2c = Ind[[2]]
  I1r = Ind[[3]]
  I2r = Ind[[4]]
  
  # Get the proxy of the joint structure by subtracting individual structure from estimated signals
  # Notice that if there is no given individual structure, the proxy joint strucutre is same as the X1 and X2
  J1c_proxy = X1 - I1c
  J2c_proxy = X2 - I2c
  J1r_proxy = X1 - I1r
  J2r_proxy = X2 - I2r
  
  # Save the svd result of the original matrices
  svd_x1 = svd(X1)
  svd_x2 = svd(X2)
  # Get the rank if r1 or r2 is not specified.
  if (is.null(r1)){
    if (method == "PL"){
      r1 = ProfileLikCluster(svd_x1$d, variance = variance1)$index
    }
    if (method == "ED"){
      r1 = Select_ED_Rank(svd_x1$d, maxiter = maxiter)
    }
  }
  if (is.null(r2)){
    if (method == "PL"){
      r2 = ProfileLikCluster(svd_x2$d, variance = variance1)$index
    }
    if (method == "ED"){
      r2 = Select_ED_Rank(svd_x2$d, maxiter = maxiter)
    }
  }
  
  # Save the column space and joint space of the proxy joint structure
  # They will be used for estimating new joint structure.
  svd_j1c = svd(J1c_proxy)
  svd_j2c = svd(J2c_proxy)
  svd_j1r = svd(J1r_proxy)
  svd_j2r = svd(J2r_proxy)
  
  # Get the estimated column space or row space.
  # If Ind is null, that means we are not sure about the joint structure
  if (isnull_Ind){
    J1c_est = as.matrix(svd_j1c$u[,1:r1])
    J2c_est = as.matrix(svd_j2c$u[,1:r2])
    J1r_est = as.matrix(svd_j1r$v[,1:r1])
    J2r_est = as.matrix(svd_j2r$v[,1:r2])
    # For column space
    # Get the principal angles
    angle_result_c = angle_cal(J1c_est, J2c_est, tol = tol)
    principal_angle_c = angle_result_c$angle
    # Get the principal vectors
    pv1_c = angle_result_c$principal_vector1
    pv2_c = angle_result_c$principal_vector2
    # If the specified joint column rank is NULL. Calculate it using the PL or ED method specified.
    if (is.null(joint_rank_c)){
      joint_rank_c = joint_angle_cluster(
        principal_angle_c, angle_threshold = angle_threshold, variance = variance2, throw = throw, maxiter = maxiter)$joint_rank
    }
    
    # Get the estimated column space by averaging the smallest joint_rank_c number of principal vectors  
    if (joint_rank_c > 0){
      joint_space_c = (pv1_c[,1:joint_rank_c] + pv2_c[,1:joint_rank_c])/2
      P_c = projection(joint_space_c)
    } 
    
    # For row space
    # Get the principal angles
    angle_result_r = angle_cal(J1r_est, J2r_est, tol = tol)
    principal_angle_r = angle_result_r$angle
    # Get the principal vectors
    pv1_r = angle_result_r$principal_vector1
    pv2_r = angle_result_r$principal_vector2
    # If the specified joint row rank is NULL. Calculate it using the PL or ED method specified.
    if (is.null(joint_rank_r)){
      joint_rank_r = joint_angle_cluster(
        principal_angle_r, angle_threshold = angle_threshold, variance = variance2, throw = throw, maxiter = maxiter)$joint_rank
    }
    # Get the estimated row space by averaging the smallest joint_rank_r number of principal vectors 
    if (joint_rank_r > 0){
      joint_space_r = (pv1_r[,1:joint_rank_r] + pv2_r[,1:joint_rank_r])/2
      P_r = projection(joint_space_r)
    }
  }
  # However, if Ind is not null, then we are pretty sure about the joint structure
  # We only use the first r_r or r_c number of basis.
  # The basis vectors of joint space are acquired by the SVD of the concatenated matrices
  else{
    c_str = cbind(J1c_proxy,J2c_proxy)
    r_str = rbind(J1r_proxy,J2r_proxy)
    svd_c = svd(c_str)
    svd_r = svd(r_str)
    joint_space_c = as.matrix(svd_c$u[,1:joint_rank_c])
    joint_space_r = as.matrix(svd_r$v[,1:joint_rank_r])
    P_c = projection(joint_space_c)
    P_r = projection(joint_space_r)
  }
  
  # Consider the extreme cases when joint rank is 0
  if (joint_rank_c == 0 || joint_rank_r == 0){
    if (joint_rank_c == 0){
      # Both joint column and row rank are 0: joint structure is 0.
      if (joint_rank_r == 0){
        signal_mat1 = svd_recover(X_1, svd_result = svd_x1, r1)
        signal_mat2 = svd_recover(X_2, svd_result = svd_x2, r2)
        J1_c = zero_matrix
        J2_c = zero_matrix
        J1_r = zero_matrix
        J2_r = zero_matrix
        I1_c = signal_mat1 
        I2_c = signal_mat2
        I1_r = signal_mat1
        I2_r = signal_mat2
        E1 = X1 - signal_mat1
        E2 = X2 - signal_mat2
      }
      # Case when joint column rank is 0. Joint row rank > 0 
      else{
        # Use the function of simplified version with the joint row space to find out the solution of signal matrices. 
        signal_mat1 = t(FindOpt_SM(t(X1), joint_space_r, r1))
        signal_mat2 = t(FindOpt_SM(t(X2), joint_space_r, r2))
        J1_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        J2_c = matrix(rep(0,n*p), nrow = n, ncol = p)
        J1_r = signal_mat1 %*% P_r
        J2_r = signal_mat2 %*% P_r
        I1_c = signal_mat1
        I2_c = signal_mat2
        I1_r = signal_mat1 - J1_r
        I2_r = signal_mat2 - J2_r 
        E1 = X1 - signal_mat1
        E2 = X2 - signal_mat2
      }
    }
    # Joint column rank is NOT 0, joint row rank is 0
    else{
      # Use the function of simplified version with the joint column space to find out the solution of signal matrices. 
      signal_mat1 = FindOpt_SM(X1, joint_space_c, r1)
      signal_mat2 = FindOpt_SM(X2, joint_space_c, r2)
      
      J1_c = P_c %*% signal_mat1
      J2_c = P_c %*% signal_mat2
      J1_r = matrix(rep(0,n*p), nrow = n, ncol = p)
      J2_r = matrix(rep(0,n*p), nrow = n, ncol = p)
      
      I1_c = signal_mat1 - J1_c
      I2_c = signal_mat2 - J2_c
      I1_r = signal_mat1 
      I2_r = signal_mat2 
      E1 = X1 - signal_mat1
      E2 = X2 - signal_mat2
    }
  }
  else{
    # For general cases use iterative algorithm to solve the optimization problem
    result1 = FindOpt_DM_Iterative(X1, joint_space_c, joint_space_r, r1, maxiter = maxiter, tol = tol)
    result2 = FindOpt_DM_Iterative(X2, joint_space_c, joint_space_r, r2, maxiter = maxiter, tol = tol)
    signal_mat1 = result1$result
    signal_mat2 = result2$result
    
    # Get the decomposition
    J1_c = P_c %*% signal_mat1
    J2_c = P_c %*% signal_mat2
    J1_r = signal_mat1 %*% P_r
    J2_r = signal_mat2 %*% P_r
    I1_c = signal_mat1 - J1_c
    I2_c = signal_mat2 - J2_c
    I1_r = signal_mat1 - J1_r
    I2_r = signal_mat2 - J2_r 
    E1 = X1 - signal_mat1
    E2 = X2 - signal_mat2
  }
  
  # Make sure that the rownames as well as the colnames stay the same
  names_r = rownames(X1)
  names_c = colnames(X1)
  rownames(J1_c) = rownames(J2_c) = rownames(J1_r) = rownames(J2_r) = names_r
  rownames(I1_c) = rownames(I2_c) = rownames(I1_r) = rownames(I2_r) = names_r
  rownames(E1) = rownames(E2) = names_r 
  
  colnames(J1_c) = colnames(J2_c) = colnames(J1_r) = colnames(J2_r) = names_c
  colnames(I1_c) = colnames(I2_c) = colnames(I1_r) = colnames(I2_r) = names_c
  colnames(E1) = colnames(E2) = names_c
  # Prepare for output
  rank_information = list("Rank 1" = r1, "Rank 2" = r2, "Joint Column Rank" = joint_rank_c, "Joint Row Rank" = joint_rank_r)
  column_decomposition = list("Joint Column 1" = J1_c, "Individual Column 1" = I1_c, 
                              "Joint Column 2" = J2_c, "Individual Column 2" = I2_c)
  row_decomposition = list("Joint Row 1" = J1_r, "Individual Row 1" = I1_r, 
                           "Joint Row 2" = J2_r, "Individual Row 2" = I2_r)
  error = list("Error1" = E1, "Error2" = E2)
  
  return(list("Rank" = rank_information, "Column Decomposition" = column_decomposition, "Row Decomposition" = row_decomposition, "Error" = error))
}
