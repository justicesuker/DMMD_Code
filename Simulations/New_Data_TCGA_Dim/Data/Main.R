rm(list=ls())
function_path1 = "DMMDFunctions/"
function_path2 = "OtherFunctions/"
source(paste(function_path1,"Preliminary_Functions.R",sep=''))
source(paste(function_path2,"DoubleMatchedDataGen.R",sep=''))

set.seed(37)

n = 88
p = 736
temp_min = min(n,p)

# Decide to make the total_upper to be 20
# In this way, the default maximum rank assumption of Edge Distribution which is 0.1*min(n,p) is satisfied.

# Number of Replicates
nrep = 140

# Fix the total rank to be 20, 18
# Joint column rank 4 
# Joint row rank 3
total_rank1 = rep(8,nrep)
total_rank2 = rep(6,nrep)

joint_rank_col = rep(0,nrep)
joint_rank_row = rep(2,nrep)

ind_rank_col1 = total_rank1 - joint_rank_col
ind_rank_col2 = total_rank2 - joint_rank_col
ind_rank_row1 = total_rank1 - joint_rank_row
ind_rank_row2 = total_rank2 - joint_rank_row

# Signal to noise ratio 1 to 1
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
  
  # Get true joint structure. we dont record individual signals since it's simply "signal - joint"
  J1r = signal1 %*% projection(data$joint_row_space)
  J2r = signal2 %*% projection(data$joint_row_space)
  
  J1c = zero_matrix
  J2c = zero_matrix
  
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
save(X1_list, file = "Simulations/New_Data_TCGA_Dim/Data/Data1.RData")
save(X2_list, signal1_list, signal2_list, J1r_list, J2r_list, J1c_list, J2c_list, file = "Simulations/New_Data_TCGA_Dim/Data/Data2.RData")
