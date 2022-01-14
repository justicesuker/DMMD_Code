rm(list=ls())
function_path1 = "DMMDFunctions/"
function_path2 = "OtherFunctions/"
source(paste(function_path1,"Preliminary_Functions.R",sep =''))
source(paste(function_path2,"DoubleMatchedDataGen.R",sep =''))
set.seed(37)

n = 240
p = 200

# Number of Replicates
nrep = 140

# Fix the total rank to be 20, 18
# In this way, the default maximum rank assumption of Edge Distribution which is 0.1*min(n,p) is satisfied.
# Joint column rank 4 
# Joint row rank 3
total_rank1 = rep(20,nrep)
total_rank2 = rep(18,nrep)
min_total = rep(18,nrep)

# Denote r_c, r_r to be joint column rank and joint row rank respectively
# We separate 140 replications as four special cases: 

# 1) 1-35: r_c, r_r = 0; both joint row rank and column rank 0
# 2) 36-70: r_c = 0, r_r = min_total; joint column rank 0, one of individual row rank 0
# 3) 71_105: r_c = min_total, r_r = 0; joint row rank 0, one of individual column rank 0
# 4) 106-140: r_c, r_r = min_total; one of individual column and row rank 0

# Get the joint rank
joint_rank_col = c(rep(0, nrep/2), min_total[(nrep/2+1):nrep])
joint_rank_row = c(rep(0, nrep/4), min_total[(nrep/4+1):(nrep/2)], rep(0, nrep/4), min_total[(3*nrep/4+1):nrep])

# Individual rank vectors
ind_rank_col1 = total_rank1 - joint_rank_col
ind_rank_col2 = total_rank2 - joint_rank_col
ind_rank_row1 = total_rank1 - joint_rank_row
ind_rank_row2 = total_rank2 - joint_rank_row

# The variance of noise matrix that makes noise signal ratio 1.
std_est1 = sqrt(total_rank1/(n*p))
std_est2 = sqrt(total_rank2/(n*p))

# Initialize
signal1_list = list()
signal2_list = list()
X1_list = list()
X2_list = list()

J1c_list = list()
J2c_list = list()
J1r_list = list()
J2r_list = list()

zero_matrix = matrix(rep(0, n*p), nrow = n, ncol = p)

# Main for loop
for(i in 1:nrep){
  # Generate data
  data = DoubleDataGen3(n = n, p = p, rank = c(total_rank1[i], total_rank2[i]), joint_rank_col = joint_rank_col[i], joint_rank_row = joint_rank_row[i], nrep = 1, std1 = std_est1[i], std2 = std_est2[i])
  
  # Get noisy double matched data
  X1 = data$X1_list[[1]]
  X2 = data$X2_list[[1]]
  
  # Get true signals
  signal1 =  data$Signal1_list[[1]]
  signal2 =  data$Signal2_list[[1]]
  
  # Based on the replication, decide whether joint structure is 0 or projection of signal onto joint space.
  # 1-35: r_c, r_r = 0; both joint row rank and column structure 0
  # 36-70: r_c = 0, r_r = min_total; joint column structure 0, joint row struture is a simple projection.
  # 71_105: r_c = min_total, r_r = 0; joint row structure 0, joint column struture is a simple projection.
  # 106-140: r_c, r_r = min_total; both joint row and column struture are simple projections.
  if (i <= nrep/4){
    J1r = zero_matrix
    J2r = zero_matrix
    J1c = zero_matrix
    J2c = zero_matrix
  }
  else if (i <= nrep/2){
    J1r = signal1 %*% projection(data$joint_row_space)
    J2r = signal2 %*% projection(data$joint_row_space)
    J1c = zero_matrix
    J2c = zero_matrix
  }
  
  else if (i <= 3*nrep/4){
    J1r = zero_matrix
    J2r = zero_matrix
    J1c = projection(data$joint_col_space) %*% signal1
    J2c = projection(data$joint_col_space) %*% signal2
  }
  else{
    J1r = signal1 %*% projection(data$joint_row_space)
    J2r = signal2 %*% projection(data$joint_row_space)
    J1c = projection(data$joint_col_space) %*% signal1
    J2c = projection(data$joint_col_space) %*% signal2
  }
  
  # Append the list
  signal1_list = append(signal1_list, list(signal1))
  signal2_list = append(signal2_list, list(signal2))
  X1_list = append(X1_list, list(X1))
  X2_list = append(X2_list, list(X2))
  
  J1c_list = append(J1c_list, list(J1c))
  J2c_list = append(J2c_list, list(J2c))
  J1r_list = append(J1r_list, list(J1r))
  J2r_list = append(J2r_list, list(J2r))
  
}

# Save the data
# Separate into two. If save them together, the size will exceed 100mb
save(signal1_list, signal2_list, X1_list, file = "Simulations/SimulationData_Setting6/Data1.RData")

save(X2_list, J1r_list, J2r_list, J1c_list, J2c_list, total_rank1, total_rank2, joint_rank_col, joint_rank_row, file = "Simulations/SimulationData_Setting6/Data2.RData")

