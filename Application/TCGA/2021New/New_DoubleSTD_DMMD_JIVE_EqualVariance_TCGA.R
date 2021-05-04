rm(list = ls())
load("Application/TCGA/TCGA_BRCA_New_Data.RData")
function_path = "Simulations/MyFunction/"

source(paste(function_path,"Angle_Calculation.R",sep=''))
source(paste(function_path,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path,"DoubleMatchedDataGen.R",sep=''))
source(paste(function_path,"FindOptMatrix.R",sep=''))
source(paste(function_path,"Preliminary_Functions.R",sep=''))
source(paste(function_path,"Select_ED_Rank.R",sep=''))
library(r.jive)

# Before doing center and scale, check if some of the rows are of variance 0.
zero_variance_index = c()
for (i in 1:dim(cancer_data)[1]){
  if (sd(normal_data[i,]) == 0 | sd(cancer_data[i,]) == 0){
    zero_variance_index = append(zero_variance_index,i)
  }
}
length(zero_variance_index)
# Get rid of those rows
cancer_data = cancer_data[-zero_variance_index,]
normal_data = normal_data[-zero_variance_index,]

# Check variance of columns
zero_variance_index_col = c()
for (i in 1:dim(cancer_data)[2]){
  if (sd(normal_data[,i]) == 0 | sd(cancer_data[,i]) == 0){
    zero_variance_index_col = append(zero_variance_index_col,i)
  }
}
zero_variance_index_col
# Double Standardized the cancer data
dcs_cancer = DoubleStandardize(cancer_data)$Result

# Check the singular values
plot(svd(dcs_cancer)$d)

# Double Standardized the normal data
dcs_normal = DoubleStandardize(normal_data)$Result

# Check the singular values
plot(svd(dcs_normal)$d)
# Make the column names match
colnames(dcs_normal) = colnames(dcs_cancer)

# DMMD result on the dataset -- all equal variance assumption:
# $`Rank 1`
# [1] 8
# 
# $`Rank 2`
# [1] 6
# 
# $`Joint Column Rank`
# [1] 2
# 
# $`Joint Row Rank`
# [1] 0

# Use the default assumption of equal variance.
result = DMMD_v2(dcs_cancer,dcs_normal)

data_log_row = list(dcs_cancer,dcs_normal)
data_log_col = list(t(dcs_cancer),t(dcs_normal))

set.seed(37)
# JIVE requires the match of the column dimension, here the subjects.
# This calculates the row decomposition of the p by n matrices we use here. 
# Corresponding to the column decomposition of n by p matrices.
# Quick, joint rank 0, individual rank 11 9
result_row = jive(data_log_row, scale = FALSE, center = FALSE)

# time: ~ 1 min. Final joint rank: 2 , final individual ranks: 12 9 
result_col = jive(data_log_col, scale = FALSE, center = FALSE)

save(result,result_row,result_col,dcs_cancer,dcs_normal,file = "Application/TCGA/2021New/DoubleSTD_DMMD_JIVE_EqualVaraince_Result.RData")


