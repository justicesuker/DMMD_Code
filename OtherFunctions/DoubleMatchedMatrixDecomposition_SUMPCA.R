# DMMD function that uses SUM-PCA to calculate joint structures. The goal of this function is to compare SUM-PCA with averaging principal vectors. We found it inferior.
# So please use DMMD_v2 or DMMD_i
DMMD_sumPCA <- function(X1, X2, r1 = NULL, r2 = NULL, joint_rank_c = NULL, joint_rank_r = NULL, variance1 = "equal", variance2 = "equal", method = "PL", tol = .Machine$double.eps^0.5, maxiter = 1e3){
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
  X_c = cbind(X1, X2)
  X_r = cbind(t(X1), t(X2))
  c_result = SUM_PCA_joint(X_c, k = 2, joint_rank = joint_rank_c)
  r_result = SUM_PCA_joint(X_r, k = 2, joint_rank = joint_rank_r)
  joint_space_c = c_result$result
  joint_rank_c = c_result$joint_rank
  joint_space_r = r_result$result
  joint_rank_r = r_result$joint_rank
  P_c = projection(joint_space_c)
  P_r = projection(joint_space_r)
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