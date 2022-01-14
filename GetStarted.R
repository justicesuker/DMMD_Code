function_path1 = "DMMDFunctions/"
function_path2 = "OtherFunctions/"
source(paste(function_path1,"Angle_Calculation.R",sep=''))
source(paste(function_path1,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path1,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path1,"DMMD_iterative.R",sep=''))
source(paste(function_path1,"FindOptMatrix.R",sep=''))
source(paste(function_path1,"Preliminary_Functions.R",sep=''))
source(paste(function_path2,"DoubleMatchedDataGen.R",sep=''))
set.seed(37)
# Generate data
data = DoubleDataGen3(n = 20, p = 16, rank = c(5, 4), joint_rank_col = 2, joint_rank_row = 1, nrep = 1, std1 = 0.01, std2 = 0.01)
X1 = data$X1_list[[1]]
X2 = data$X2_list[[1]]
# Apply DMMD
DMMD_result = DMMD_v2(X1, X2)
# Apply DMMD-i
DMMD_i_result = DMMD_i(X1, X2)
# Get estimated 1st signal matrix of DMMD
A1_est = X1 - DMMD_result$Error$Error1
A1 = data$Signal1_list[[1]]
# Check the difference between estimated signal of DMMD with true signal 
head(A1_est - A1)
# Get estimated 1st signal matrix of DMMD-i
A1_est_i = DMMD_i_result$A1
# Check the difference between estimated signal of DMMD-i with true signal 
head(A1_est_i - A1)
