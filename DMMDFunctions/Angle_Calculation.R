#' Function that calculates principal angles between column spaces of two matrices
#'
#' @param X The first matrix.
#' @param Y The second matrix.
#' @param tol Tolerence, default is the square root of machine precision.
#'
#' @return A list that contains the following:
#' \item{angle}{A vector of principal angles with increasing order.}
#' \item{cos_angle}{A vector of cosine principal angles}
#' \item{principal_vector1}{Principal vectors of matrix \code{X}}
#' \item{principal_vector2}{Principal vectors of matrix \code{Y}}
#'
#' @examples
#' X = matrix(c(1,1,1,1,1,0),nrow = 3, ncol = 2)
#' Y = matrix(c(1,1,1,2,1,0),nrow = 3, ncol = 2)
#' angle_cal(X,Y)
#' 
angle_cal <- function(X, Y, tol = .Machine$double.eps^0.5){
  X = as.matrix(X)
  Y = as.matrix(Y)
  Xnorm = svd(X)$u
  Ynorm = svd(Y)$u
  M = crossprod(Xnorm,Ynorm)
  # Extreme case when both X and Y only contain one number
  if (dim(M)[1] == 1 && dim(M)[2] == 1){
    cos_angle = abs(M)
    principal_angle = NA
    if (cos_angle >= 1){principal_angle = 0}
    if (cos_angle <= 0){principal_angle = pi/2}
    if (cos_angle > 0 && cos_angle < 1){
      principal_angle = acos(cos_angle)
    }
    principal_mat1 = Xnorm 
    principal_mat2 = Ynorm 
    return(list("angle" = principal_angle, "cos_angle" = cos_angle,
              "principal_vector1" = principal_mat1, 
              "principal_vector2" = principal_mat2))
  }
  # Normal case when X and Y are matrices (data frames)
  else{
    svd_result = svd(M)
    cos_angle = svd_result$d
    l = length(cos_angle)
    principal_angle = rep(NA, l)
    for (i in 1:l){
      if (cos_angle[i] >= 1){principal_angle[i] = 0}
      if (cos_angle[i] <= 0){principal_angle[i] = pi/2}
      if (cos_angle[i] > 0 && cos_angle[i] < 1){
        principal_angle[i] = acos(cos_angle[i])
      }
    }
    principal_mat1 = Xnorm %*% svd_result$u
    principal_mat2 = Ynorm %*% svd_result$v
    return(list("angle" = principal_angle, "cos_angle" = cos_angle,
                "principal_vector1" = principal_mat1, 
                "principal_vector2" = principal_mat2))
  }
}

#' Function that estimates joint rank by method of profile likelihood.
#' @param angle_vec A vector of principal angles
#' @param angle_threshold Optional Threshold (radians) that is used to truncate principal angles. Default is 45 degrees.
#' @param variance Either "equal" or "unequal", i.e whether the assumption is equal variance or unequal variance. Default is unequal.
#' @param cos_method TRUE or FALSE, i.e., whether use angles or cosine angles as the data to cluster. Default is FALSE.
#'
#' @return A list with the following elements:
#' \item{joint_rank}{Estimated joint rank.}
#' \item{profileloglikvec}{Profile log likelihood calculated at each index.
#' The function will return NA with a warning message, if less or equal to 2 principal angles are smaller than the threshold.}
#'
#' @examples
#' 
#' data = datagen(n = 100, p = 20, joint_rank = 5, individual_rank = c(5,7), nrep = 1)
#' X1 = data$X1_list[[1]]
#' X2 = data$X2_list[[1]]
#' angle_vec = angle_cal(X1,X2)$angle
#' joint_angle_cluster(angle_vec)
#' 
joint_angle_cluster <- function(angle_vec, angle_threshold = 90 * pi/180, variance = "unequal", throw = FALSE, maxiter = 100){
  if (angle_threshold > 90 * pi/180){
    stop("Angle threshold cannot exceed 90 degrees.")
  }
  if (angle_threshold < 0){
    stop("Angle threshold cannot be less than 0 degree.")
  }
  L = length(angle_vec)
  # Get the angles that are smaller than the specified threshold. 
  small_angle_vec = angle_vec[angle_vec < angle_threshold]
  l = length(small_angle_vec)
  # Extreme cases 
  if (l <= 2){
    # No small angle. Joint rank = 0
    if (l == 0){
      warning("There is no angle that is smaller than threshold, the joint rank is 0.")
      return(list("joint_rank" = 0, "angle_vec_considered" = NA, "profileloglikvec" = NA))
    }
    # Only one small angle. Joint rank = 1.
    if (l == 1){
      if (small_angle_vec < pi/4){
        warning("There is only one angle that is smaller than 45 degrees, the joint rank is recommended as 1 without using profile likelihood method.")
        return(list("joint_rank" = 1, "angle_vec_considered" = NA, "profileloglikvec" = NA))
      }
      else{
        warning("There is only one angle that is smaller than threshold, but the angle is greater than 45 degrees. The joint rank is recommended as 0 without using profile likelihood method.")
        return(list("joint_rank" = 0, "angle_vec_considered" = NA, "profileloglikvec" = NA))
      }

    }
    # When there are only two angles that are small
    if (l == 2){
      # Two really small angles, joint rank = 2
      if (small_angle_vec[1] <= pi/4 & small_angle_vec[2] <= pi/4){
        warning("There are only two angles that are smaller than threshold, and they both are smaller than 45 degrees. The joint rank is recommended to be 2 without using profile likelihood method.")
        return(list("joint_rank" = 2, "angle_vec_considered" = NA, "profileloglikvec" = NA))
      }
      # One angle > pi/4, one < pi/4. Automatically output two clusters with one angle in each, joint rank = 1.
      if (small_angle_vec[1] <= pi/4 & small_angle_vec[2] > pi/4){
        warning("There are only two angles that are smaller than threshold. One of them is smaller than 45 degrees, while the other is greater than 45 degrees. 
                The joint rank is recommended to be 1 without using profile likelihood method.")
        return(list("joint_rank" = 1, "angle_vec_considered" = NA, "profileloglikvec" = NA))
      }
      # Two angles > pi/4, joint rank = 0
      if (small_angle_vec[1] > pi/4 & small_angle_vec[2] > pi/4){
        warning("There are only two angles that are smaller than threshold, but both of them are greater than 45 degrees. 
                The joint rank is recommended to be 0 without using profile likelihood method.")
        return(list("joint_rank" = 0, "angle_vec_considered" = NA, "profileloglikvec" = NA))
      }
    }
  }
  # When there are 3 or more angles that could possibly considered joint signal,
  else{
    if (var(small_angle_vec) < .Machine$double.eps^0.5){
      warning("The variance of angles that are smaller than threshold is too small. It is treated that all of them are joint signals. Please check manually.")
      return(list("joint_rank" = l, "angle_vec_considered" = NA, "profileloglikvec" = NA))
    }
    else{
      # If throw is T, we need to only consider the small angles, and apply profile likelihood (PL) or edge distribution (ED) method on them.
      if (throw){
        # Add 0 and angle_threshold artificially to include 0 joint rank case and full joint rank case
        small_angle_vec = c(0,small_angle_vec,angle_threshold)
        result = ProfileLikCluster(small_angle_vec, variance = variance)
        index = result$index - 1
        profileloglikvec = result$profileloglikvec
        return(list("joint_rank" = index, "angle_vec_considered" = small_angle_vec, "profileloglikvec" = profileloglikvec))
      }
      # If throw is F, we need to consider all the principal angles, and apply profile likelihood (PL) or edge distribution (ED) method on them.
      # But the joint rank could not exceed the number of small angles.
      else{
        # Similarly, add 0 artificially to include 0 joint rank case and full joint rank case
        angle_vec = c(0,angle_vec)
        # If all the principal angles are small, we need to add the threshold as well to cover the case when joint rank = total rank.
        if (L == l){angle_vec = c(angle_vec, angle_threshold)}
        # Apply PL on all principal angles
        # Notice maxindex is l+1 since we add 0 angle.
        result = ProfileLikCluster(angle_vec, variance = variance, maxindex = l + 1)
        index = result$index - 1
        profileloglikvec = result$profileloglikvec
        return(list("joint_rank" = index, "angle_vec_considered" = small_angle_vec, "profileloglikvec" = profileloglikvec))
      }
    }
  }
}

#' Function that estimates joint row and column space by SUM-PCA.
#' X is the concatenated matrices (matched by rows)
SUM_PCA_joint <- function(X, k = 2, joint_rank = NULL){
  n = dim(X)[1]
  p = dim(X)[2]/k
  if (!is.null(joint_rank)){
    if (joint_rank > min(n,p)){
      stop("The joint rank input is greater than n or p, which is illegal.")
    }
  }
  svd_x = svd(X)
  if (is.null(joint_rank)){
    joint_rank = ProfileLikCluster(svd_x$d, variance = 'equal', maxindex = min(n,p))$index
  }
  result = svd_x$u[,1:joint_rank, drop = FALSE]
  return(list('result' = result, 'joint_rank' = joint_rank))
}