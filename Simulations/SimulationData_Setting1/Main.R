rm(list=ls())
function_path = "Simulations/MyFunction/"
source(paste(function_path,"DoubleMatchedDataGen.R",sep=''))
source(paste(function_path,"Preliminary_Functions.R",sep=''))
set.seed(37)

n = 240
p = 200
temp_min = min(n,p)

# Decide to make the total_upper to be 20
# In this way, the default maximum rank assumption of Edge Distribution which is 0.1*min(n,p) is satisfied.

# Upperbound for total rank. 
total_upper = 0.1*min(n,p)

# Number of Replicates
nrep = 140

# Total rank vectors
total_rank1 = rep(NA,nrep)
total_rank2 = rep(NA,nrep)

# Joint rank vectors
joint_rank_col = rep(NA,nrep)
joint_rank_row = rep(NA,nrep)

# Get the total rank and joint rank that satisfies the constraint when generating data.
# Refer to the manuscript for those constraints. Section 3.1
for (i in 1:nrep){
  temp_lower = -Inf
  temp_upper = +Inf
  while(temp_lower <= (temp_upper - temp_min/4 + 1)){
    temp1 = sample(2:total_upper, 1)
    temp2 = sample(2:total_upper, 1)
    temp_lower = min(temp1, temp2)
    temp_upper = max(temp1, temp2)
  }
  total_rank1[i] = temp1
  total_rank2[i] = temp2
  
  # Let's make the upper bound for joint rank be 5. Since we dont need to hard threshold it according to ED.
  temp_joint_upper = min((temp_lower - 1),temp_min/2, 5)
  temp_joint_col_lower = max(1, (temp_upper - n/4))
  temp_joint_row_lower = max(1, (temp_upper - p/4))

  joint_rank_col[i] = sample(temp_joint_col_lower:temp_joint_upper,1)
  joint_rank_row[i] = sample(temp_joint_row_lower:temp_joint_upper,1)
}

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
  J1c = projection(data$joint_col_space) %*% signal1
  J2c = projection(data$joint_col_space) %*% signal2
  
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
save(signal1_list, signal2_list, X1_list, X2_list, J1r_list, J2r_list, 
     J1c_list, J2c_list, total_rank1, total_rank2, joint_rank_col, joint_rank_row, file = "Simulations/SimulationData_Setting1/Data.RData")

