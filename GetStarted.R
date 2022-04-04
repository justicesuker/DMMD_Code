rm(list = ls())
source("OtherFunctions/DoubleMatchedDataGen.R")
source("DMMDFunctions/Preliminary_Functions.R")
set.seed(37)
# Generate data
data = DoubleDataGen3(n = 20, p = 16, rank = c(4, 3), joint_rank_col = 2, joint_rank_row = 1, nrep = 1, std1 = 0.01, std2 = 0.01)
X1 = data$X1_list[[1]]
X2 = data$X2_list[[1]]

# Source all necessary functions to apply DMMD
source("DMMDFunctions/Angle_Calculation.R")
source("DMMDFunctions/Profile_Likelihood_Rank_Selection.R")
source("DMMDFunctions/DoubleMatchedMatrixDecomposition.R")
source("DMMDFunctions/DMMD_iterative.R")
source("DMMDFunctions/FindOptMatrix.R")
# Apply DMMD
DMMD_result = DMMD_v2(X1, X2)
# Extract estimated ranks
DMMD_result$"Rank"
# Extract common and individual structures for 1st signal
J1c = DMMD_result$"Column Decomposition"$"Joint Column 1"
I1c = DMMD_result$"Column Decomposition"$"Individual Column 1"
J1r = DMMD_result$"Row Decomposition"$"Joint Row 1"
I1r = DMMD_result$"Row Decomposition"$"Individual Row 1"
# Verify that both decompositions give the same estimated signal
sum((J1c + I1c - J1r - I1r)^2)
# Get estimated 1st signal matrix of DMMD
A1_est = X1 - DMMD_result$Error$Error1
# Get the true signal
A1 = data$Signal1_list[[1]]
# Check the difference between estimated signal of DMMD with true signal 
sum((A1_est - A1)^2)

# Apply DMMD-i
DMMD_i_result = DMMD_i(X1, X2)
# Get estimated signal matrices of DMMD-i
A1_est_i = DMMD_i_result$A1
A2_est_i = DMMD_i_result$A2
# Get estimated basis vectors for joint subspaces
M_est = DMMD_i_result$M # joint column space (matched samples)
N_est = DMMD_i_result$N # joint row space (matched features)
# Get the joint column structure for view 1. 
J1c = M_est %*% t(M_est) %*% X1
# Get individual column structure
I1c = A1_est_i - J1c
# Get the joint row structure for view 1
J1r =  X1 %*% N_est %*% t(N_est) 
# Get individual crow structure
I1r = A1_est_i - J1r
# Check the difference between estimated signal of DMMD-i with true signal
sum((A1_est_i - A1)^2)