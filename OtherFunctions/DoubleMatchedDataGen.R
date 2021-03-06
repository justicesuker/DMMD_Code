# New function that generates double-matched data with different singular values.
DoubleDataGen3 <- function(n = 1000, p = 400, rank = c(130,150), joint_rank_col = 30, joint_rank_row = 50, nrep = 20, noise = TRUE, std1 = 0.1, std2 = 0.1, lb = 0.5, ub = 1.5){
  # Check if the inputs meet some requirements.
  if (n %% 4 != 0){stop("4 must be a factor of n.")}
  if (p %% 4 != 0){stop("4 must be a factor of p.")}
  if (max(rank) > min(n, p)){stop("Total rank exceeds the size of matrices.")}
  if (max(joint_rank_col, joint_rank_row) > max(rank)){stop("Joint rank exceeds the total rank.")}
  # Initialize the output of (noisy) double-matched matrices
  X1_list = list()
  X2_list = list()
  # Initialize the output of correponding signal matrices 
  Signal1_list = list()
  Signal2_list = list()
  for (i in 1:nrep){
    # Get the individual ranks
    individual_rank_col = rank - joint_rank_col
    individual_rank_row = rank - joint_rank_row
    # Get the index where 1 occurs in the standard basis for the joint column space
    joint_index_col = sample(1:(n/2),joint_rank_col)
    # Get the index where 1 occurs in the standard basis for the joint row space
    joint_index_row = sample(1:(p/2),joint_rank_row)
    # Get the index where 1 occurs in the standard basis for the individual column space 1.
    individual_index1_col = sample(1:(n/4),individual_rank_col[1])
    # Get the index where 1 occurs in the standard basis for the individual column space 2.
    individual_index2_col = sample(1:(n/4),individual_rank_col[2])
    # Get the index where 1 occurs in the standard basis for the individual row space 1.
    individual_index1_row = sample(1:(p/4),individual_rank_row[1])
    # Get the index where 1 occurs in the standard basis for the individual row space 2.
    individual_index2_row = sample(1:(p/4),individual_rank_row[2])
    
    # Initialize column space
    # The following steps only concerns the first basis for corresponding column space, which will be used in the possible later cbind function.
    # Joint column space
    tempv_col = rep(0,n)
    tempv_col[joint_index_col[1]] = 1
    joint_matrix_col = as.matrix(tempv_col)
    
    # Individual column space 1
    tempv21_col = rep(0,n)
    # the index should add n/2 so that the individual space 1 is not shared with joint space
    tempv21_col[individual_index1_col[1]+(n/2)] = 1
    individual_matrix1_col = as.matrix(tempv21_col)
    
    # Individual column space 2
    tempv22_col = rep(0,n)
    # the index should add 3*n/4 so that the individual space 2 is not shared with joint space
    tempv22_col[individual_index2_col[1]+ (3/4*n)] = 1
    individual_matrix2_col = as.matrix(tempv22_col)
    
    # Full joint column space
    # If the joint column rank is not 1, we need to do combine all the basis to be the matrix. 
    if (joint_rank_col >= 2){
      for (i in 2:joint_rank_col){
        tempv_col = rep(0,n)
        tempv_col[joint_index_col[i]] = 1
        joint_matrix_col = cbind(joint_matrix_col,tempv_col)
      }
    }
    # Full individual column space 1
    # If the individual column rank 1 is not 1, we need to do combine all the basis to be the matrix. 
    if (individual_rank_col[1] >= 2){
      for (i in 2:individual_rank_col[1]){
        tempv21_col = rep(0,n)
        tempv21_col[individual_index1_col[i]+(n/2)] = 1
        individual_matrix1_col = cbind(individual_matrix1_col,tempv21_col)
      }
    }
    # Full individual column space 2
    # If the individual column rank 2 is not 1, we need to do combine all the basis to be the matrix. 
    if (individual_rank_col[2] >= 2){
      for (i in 2:individual_rank_col[2]){
        tempv22_col = rep(0,n)
        tempv22_col[individual_index2_col[i]+(3/4*n)] = 1
        individual_matrix2_col = cbind(individual_matrix2_col,tempv22_col)
      }
    }
    
    # Some extreme cases.
    if (joint_rank_col == 0){
      col_space1 = individual_matrix1_col
      col_space2 = individual_matrix2_col
    }
    if (individual_rank_col[1] == 0){
      col_space1 = joint_matrix_col
    }
    if (individual_rank_col[2] == 0){
      col_space2 = joint_matrix_col
    }
    if (individual_rank_col[1] != 0 & joint_rank_col != 0){
      col_space1 = cbind(joint_matrix_col,individual_matrix1_col)
    }
    if (individual_rank_col[2] != 0 & joint_rank_col != 0){
      col_space2 = cbind(joint_matrix_col,individual_matrix2_col)  
    }
    
    # Initialize row space
    # The following steps only concerns the first basis for corresponding row space, which will be used in the possible later cbind function.
    # Joint row space
    tempv_row = rep(0,p)
    tempv_row[joint_index_row[1]] = 1
    joint_matrix_row = as.matrix(tempv_row)
    
    # Individual row space 1
    tempv21_row = rep(0,p)
    # the index should add p/2 so that the individual space 1 is not shared with joint space
    tempv21_row[individual_index1_row[1]+(p/2)] = 1
    individual_matrix1_row = as.matrix(tempv21_row)
    
    # Individual row space 2
    tempv22_row = rep(0,p)
    # the index should add 3*p/4 so that the individual space 2 is not shared with joint space
    tempv22_row[individual_index2_row[1]+ (3/4*p)] = 1
    individual_matrix2_row = as.matrix(tempv22_row)
    
    # Full joint row space
    # If the joint row rank is not 1, we need to do combine all the basis to be the matrix. 
    if (joint_rank_row >= 2){
      for (i in 2:joint_rank_row){
        tempv_row = rep(0,p)
        tempv_row[joint_index_row[i]] = 1
        joint_matrix_row = cbind(joint_matrix_row,tempv_row)
      }
    }
    
    # Full individual row space 1
    # If the individual row rank 1 is not 1, we need to do combine all the basis to be the matrix. 
    if (individual_rank_row[1] >= 2){
      for (i in 2:individual_rank_row[1]){
        tempv21_row = rep(0,p)
        tempv21_row[individual_index1_row[i]+(p/2)] = 1
        individual_matrix1_row = cbind(individual_matrix1_row,tempv21_row)
      }
    }
    
    # Full individual row space 2
    # If the individual row rank 2 is not 1, we need to do combine all the basis to be the matrix. 
    if (individual_rank_row[2] >= 2){
      for (i in 2:individual_rank_row[2]){
        tempv22_row = rep(0,p)
        tempv22_row[individual_index2_row[i]+(3/4*p)] = 1
        individual_matrix2_row = cbind(individual_matrix2_row,tempv22_row)
      }
    }
    
    # Some extreme cases.
    if (joint_rank_row == 0){
      row_space1 = individual_matrix1_row
      row_space2 = individual_matrix2_row
    }
    if (individual_rank_row[1] == 0){
      row_space1 = joint_matrix_row
    }
    if (individual_rank_row[2] == 0){
      row_space2 = joint_matrix_row
    }
    if (individual_rank_row[1] != 0 & joint_rank_row != 0){
      row_space1 = cbind(joint_matrix_row,individual_matrix1_row)
    }
    if (individual_rank_row[2] != 0 & joint_rank_row != 0){
      row_space2 = cbind(joint_matrix_row,individual_matrix2_row)  
    }
    
    # Prepare the orthonormal transformation matrix
    R1_col = GenOrthoMatrix(rank[1])
    R1_row = GenOrthoMatrix(rank[1])
    R2_col = GenOrthoMatrix(rank[2])
    R2_row = GenOrthoMatrix(rank[2])
    
    # Vecotr of Singular values
    # They are from a uniform distribution to make sure they are different
    Dvec1 = runif(rank[1],min = lb, max = ub)
    # Scale the singular values so that their sum of square equals to the rank.
    Dvec1 = Dvec1/sqrt(sum(Dvec1^2))*sqrt(rank[1])
    if (rank[1] >= 2){
      Dvec1 = diag(sort(Dvec1, decreasing = TRUE))
    }
    if (rank[1] < 2){
      Dvec1 = as.matrix(Dvec1)
    }
    Dvec2 = runif(rank[2],min = lb, max = ub)
    Dvec2 = Dvec2/sqrt(sum(Dvec2^2))*sqrt(rank[2])
    if (rank[2] >= 2){
      Dvec2 = diag(sort(Dvec2, decreasing = TRUE))
    }
    if (rank[2] < 2){
      Dvec2 = as.matrix(Dvec2)
    }
    
    # Transform the original standard basis containing only 0,1 to non-standard basis and get the signal matrices
    Signal1 = col_space1 %*% R1_col %*% Dvec1 %*% R1_row %*% t(row_space1)
    Signal2 = col_space2 %*% R2_col %*% Dvec2 %*% R2_row %*% t(row_space2)
    
    # Append the output list
    Signal1_list = append(Signal1_list, list(Signal1))
    Signal2_list = append(Signal2_list, list(Signal2))
    
    # Add noise if noise = TRUE
    if (noise){
      E1 = matrix(rnorm(n*p, mean = 0, sd = std1), nrow = n)
      E2 = matrix(rnorm(n*p, mean = 0, sd = std2), nrow = n)
    }
    
    # Get the noisy matrices
    X1 = Signal1 + E1
    X2 = Signal2 + E2
    
    # Append the output list
    X1_list = append(X1_list, list(X1))
    X2_list = append(X2_list, list(X2))
  }
  return(list("X1_list" = X1_list, "X2_list" = X2_list, 
              "Signal1_list" = Signal1_list, "Signal2_list" = Signal2_list,
              "joint_row_space" = joint_matrix_row, "joint_col_space" = joint_matrix_col,
              "ind_row_space1" = individual_matrix1_row, "ind_row_space2" = individual_matrix2_row,
              "ind_col_space1" = individual_matrix1_col, "ind_col_space2" = individual_matrix2_col))
}