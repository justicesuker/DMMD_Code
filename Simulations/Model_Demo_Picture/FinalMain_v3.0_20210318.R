rm(list=ls())
function_path1 = "DMMDFunctions/"
function_path2 = "OtherFunctions/"
source(paste(function_path1,"Angle_Calculation.R",sep=''))
source(paste(function_path1,"Profile_Likelihood_Rank_Selection.R",sep=''))
source(paste(function_path1,"DoubleMatchedMatrixDecomposition.R",sep=''))
source(paste(function_path1,"FindOptMatrix.R",sep=''))
source(paste(function_path1,"Preliminary_Functions.R",sep=''))
source(paste(function_path2,"DoubleMatchedDataGen.R",sep=''))
library(r.jive)
set.seed(37)

n = 80
p = 40
temp_min = min(n,p)

# Decide to make the total_upper to be 20
# In this way, the default maximum rank assumption of Edge Distribution which is 0.1*min(n,p) is satisfied.

# Upperbound for total rank. 
total_upper = 0.1*min(n,p)

# Number of Replicates
nrep = 1

# Total rank vectors
total_rank1 = 15
total_rank2 = 12

# Joint rank vectors
joint_rank_col = 7
joint_rank_row = 5

# Individual rank vectors
ind_rank_col1 = total_rank1 - joint_rank_col
ind_rank_col2 = total_rank2 - joint_rank_col
ind_rank_row1 = total_rank1 - joint_rank_row
ind_rank_row2 = total_rank2 - joint_rank_row

# The variance of noise matrix that makes noise signal ratio 1.
std_est1 = sqrt(total_rank1/(n*p))
std_est2 = sqrt(total_rank2/(n*p))

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
  
  I1r = signal1 - J1r
  I2r = signal2 - J2r
  I1c = signal1 - J1c
  I2c = signal2 - J2c
  
  E1 = X1 - signal1
  E2 = X2 - signal2
}

my_result = DMMD_v2(X1, X2, r1 = total_rank1, r2 = total_rank2, joint_rank_c = joint_rank_col, joint_rank_r = joint_rank_row)

E1_est = my_result$Error$Error1
E2_est = my_result$Error$Error2

J1r_est = my_result$`Row Decomposition`$`Joint Row 1`
J2r_est = my_result$`Row Decomposition`$`Joint Row 2`

J1c_est = my_result$`Column Decomposition`$`Joint Column 1`
J2c_est = my_result$`Column Decomposition`$`Joint Column 2`

I1r_est = my_result$`Row Decomposition`$`Individual Row 1`
I2r_est = my_result$`Row Decomposition`$`Individual Row 2`

I1c_est = my_result$`Column Decomposition`$`Individual Column 1`
I2c_est = my_result$`Column Decomposition`$`Individual Column 2`

library(cowplot)
library(dplyr)
library(MASS)
library(forcats)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(rafalib)
library(gplots)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

range_max = max(max(X1),max(X2),max(J1c),max(J2c),
                max(I1c),max(I2c),max(J1r),max(J2r),max(I1r),max(I2r),max(E1),max(E2))

range_min = min(min(X1),min(X2),min(J1c),min(J2c),
                min(I1c),min(I2c),min(J1r),min(J2r),min(I1r),min(I2r),min(E1),min(E2))

final_range = max(abs(range_max),abs(range_min))

hmcol <- colorRamp2(c(-final_range, 0, final_range), c("blue ", "white", "red"))

fig.path = "Simulations/Model_Demo_Picture/FinalFigures2021/"

X_tot = rbind(X1,X2)
J_c = rbind(J1c,J2c)
I_c = rbind(I1c,I2c)
J_r = rbind(J1r,J2r)
I_r = rbind(I1r,I2r)
E = rbind(E1,E2)

g1 = Heatmap(X_tot, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, row_title = c(expression(X[1]), expression(X[2])), 
             show_column_names = FALSE, row_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)),row_title_rot = 0)

g2 = Heatmap(J_c, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(J[c]),
             show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g3 = Heatmap(I_c, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(I[c]),
             show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g4 = Heatmap(E, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = "Noise",
             show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g5 = Heatmap(J_r, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(J[r]),
             show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g6 = Heatmap(I_r, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(I[r]),
             show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

pdf(file = paste(fig.path,"Column_Decomposition_Demo.pdf",sep=""), width = 8, height = 6)
print(g1 + g2 + g3 + g4)
dev.off()

pdf(file = paste(fig.path,"Row_Decomposition_Demo.pdf",sep=""), width = 8, height = 6)
print(g1 + g5 + g6 + g4)
dev.off()

# Plot my estimated decomposition
J_c_est = rbind(J1c_est,J2c_est)
I_c_est = rbind(I1c_est,I2c_est)
J_r_est = rbind(J1r_est,J2r_est)
I_r_est = rbind(I1r_est,I2r_est)
E_est = rbind(E1_est,E2_est)

g1_est = Heatmap(X_tot, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, row_title = c(expression(X[1]), expression(X[2])), 
                 show_column_names = FALSE, row_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)),row_title_rot = 0)

g2_est = Heatmap(J_c_est, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(J[c]),
                 show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g3_est = Heatmap(I_c_est, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(I[c]),
                 show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g4_est = Heatmap(E_est, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = "Noise",
                 show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g5_est = Heatmap(J_r_est, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(J[r]),
                 show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

g6_est = Heatmap(I_r_est, col = hmcol, border = 'black', cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, column_title = expression(I[r]),
                 show_column_names = FALSE, column_title_gp = gpar(fontsize = 30), row_split = c(rep(1,n),rep(2,n)), row_title = NULL)

pdf(file = paste(fig.path,"Est_Column_Decomposition_Demo.pdf",sep=""), width = 8, height = 6)
print(g1_est + g2_est + g3_est + g4_est)
dev.off()

pdf(file = paste(fig.path,"Est_Row_Decomposition_Demo.pdf",sep=""), width = 8, height = 6)
print(g1_est + g5_est + g6_est + g4_est)
dev.off()
