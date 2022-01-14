#' Function that calculates the optimum matrix in single matched case
#'  
#' @param X The original matrix
#' @param M The column space that the output matrix should have 
#' @param r The rank of output matrix
#' 
#' @details The resulting matrix is the solution the following problem:
#' \deqn{min_{Y}{||X - Y||^2_F}}
#' subject to \deqn{M is a subset of column space of Y}
#' 
#' @return A matrix that is the solution to the optimization problem.
#'
#' @examples 
#' X = matrix(c(2,1,1,3,2,2),nrow = 3)
#' M = matrix(c(1,1,0),nrow = 3)
#' r = 2
#' FindOpt_SM(X,M,r)
FindOpt_SM <- function(X, M, r){
  M = cbind(M)
  s = dim(M)[2]
  n = dim(M)[1]
  if (n != dim(X)[1]){stop("The dimension of M and X does not match.")}
  if (s > r){stop("The number of columns in M should not exceed r.")}
  temp = X - projection(M) %*% X
  svd_result = svd(temp)$u
  if (r == s){
    Mtilde = M
  }
  else{
    Mtilde = cbind(M,svd_result[,1:(r-s)])
  }
  result = projection(Mtilde) %*% X 
  return(result)
}

#' This function only works when either the rank of M or N equal to the specified rank r. It's a prerequisite of FindOpt_DM_Iterative
#' @param X The original matrix
#' @param M The column space that the output matrix should have 
#' @param N The row space that the solution should have
#' @param r The rank of output matrix
#' 
#' @details The resulting matrix is the solution the following problem:
#' \deqn{min_{Y}{||X - Y||^2_F}}
#' subject to \deqn{M is a subset of column space of Y}
#' 
#' @return A matrix that is the solution to the optimization problem.
FindOpt_SimplifiedCase<- function(X, M, N, r){
  r_m = dim(M)[2]
  r_n = dim(N)[2]
  m = dim(M)[1]
  n = dim(N)[1]
  projM = projection(M)
  projN = projection(N)
  if (m != dim(X)[1]){stop("The dimension of M and X does not match.")}
  if (n != dim(X)[2]){stop("The dimension of N and X does not match.")}
  if (r != max(r_m,r_n)){stop("This is not a simplied case, please use FindOpt_DM_Iterative function")}
  # Simple case
  if (r_m == r_n){
    return(projM %*% X %*% projN)
  }
  # If r_m > r_n, we need to find the rest of the row space and form them into a new matrix that contains the full row space and then calculate projM %*% X %*% proj(Nnew)
  else if (r_m > r_n){
    svd_result = svd(projM %*% X %*% (diag(n) - projN))
    S = svd_result$v[,1:(r_m - r_n)]
    sol = projM %*% X %*% (projN + projection(S,ortho = TRUE))
    return(sol)
  }
  # If r_m < r_n, we need to find the rest of the column space and form them into a new matrix that contains the full column space and then calculate proj(Mnew) %*% X %*% projN
  else{
    svd_result = svd((diag(m) - projM) %*% X %*% projN)
    R = svd_result$u[,1:(r_n - r_m)]
    sol = (projM + projection(R,ortho = TRUE)) %*% X %*% projN
    return(sol)
  }
}

# An iterative function that solves the optimization problem of (1) in the manuscript.
# M is the column space that the solution should have
# N is the row space that the solution should have
# r is the rank
FindOpt_DM_Iterative <- function(X, M, N, r, maxiter = 1e4, tol = .Machine$double.eps^0.5){
  M = as.matrix(M)
  N = as.matrix(N)
  X = as.matrix(X)
  r_m = dim(M)[2]
  r_n = dim(N)[2]
  m = dim(M)[1]
  n = dim(N)[1]
  # Output error message if something does not match.
  if (m != dim(X)[1]){stop("The dimension of M and X does not match.")}
  if (n != dim(X)[2]){stop("The dimension of N and X does not match.")}
  if (r_m > r){stop("The number of columns in M should not exceed r.")}
  if (r_n > r){stop("The number of columns in N should not exceed r.")}
  # If this is a simplified case, simply call the FindOpt_SimplifiedCase function
  if (r == max(r_m,r_n)){
    return(list(result = FindOpt_SimplifiedCase(X, M, N, r)))
  }
  # Main body of the function
  else{
    projN = projection(N)
    projM = projection(M)
    # Warm start. Initialize
    temp_init = (diag(m) - projM) %*% X
    # Initialize the rest of column space.
    R = svd(temp_init)$u[,1:(r - r_m)]
    # Initialize the Frobenius norm difference of the previous matrix and the matrix after one iteration as +inf to make sure it does not converge at the first step. 
    norm_diff = +Inf
    # Record the number of iterations
    num_iter = 1
    # Initialize the values of the objective function
    opt_function_val = c()
    while(norm_diff > tol){
      # Stop the loop if the current number of iteration exceeeds the specified maximum iteration number 
      if (num_iter > maxiter){break}
      # Calculate the best possible rest of row space for the specified column space 
      temp1 = (projM + projection(R, ortho = TRUE)) %*% X %*% (diag(n) - projN)
      Snew = svd(temp1)$v[,1:(r - r_n)]
      # Calculate the best possible rest of column space for the calculated row space 
      temp2 = (diag(m) - projM) %*% X %*% (projN + projection(Snew, ortho = TRUE))
      Rnew = svd(temp2)$u[,1:(r - r_m)]
      # If the number of iteration is 1, there is no S. Force the norm_diff = +inf for the second iteration.
      if (num_iter == 1){
        norm_diff = +Inf
      }
      # After one iteration, record Frobenius norm difference of the new column/row space and the old column/row space
      # The maximum value of the two will be used to check convergence.
      else{
        norm_diff1 = Fnorm(Snew - S)
        norm_diff2 = Fnorm(Rnew - R)
        norm_diff = max(norm_diff1,norm_diff2)
      }
      # After one iteration
      S = Snew
      R = Rnew
      result = (projM + projection(R, ortho = TRUE)) %*% X %*% (projN + projection(S, ortho = TRUE))
      # Calculate the objective value
      opt_val = Fnorm(result - X)
      opt_function_val = append(opt_function_val,opt_val)
      # iteration number + 1
      num_iter = num_iter + 1
    }
    return(list(result = result, R = R, S = S, num_iter = num_iter,opt_function_val = opt_function_val))
  }
}
