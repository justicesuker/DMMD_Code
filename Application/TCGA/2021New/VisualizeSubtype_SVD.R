rm(list = ls())
library(cowplot)
library(dplyr)
library(MASS)
library(forcats)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(rafalib)
library(gplots)
library(d3heatmap)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Load the result of DMMD and JIVE
load('Application/TCGA/2021New/DoubleSTD_DMMD_JIVE_EqualVaraince_Result.RData')

set.seed(37)

cdata = dcs_cancer
ndata = dcs_normal

# There are 734 genes and 87 samples.

# !!!Remember: cancer data first, normal data second!!!
# Remember: row space -- miRNA gene signal
#           column space -- subject signal

# The rank estimation:
# DMMD result on the dataset -- all equal variance assumption:
# r1 = 8, r2 = 6, r_c = 2, r_r = 0

# Specify path to save figures
results.path = 'Application/TCGA/2021New'
figure.path = paste(results.path, "Cluster_SVD/",sep = "/")

hmcol <- colorRamp2(c(-3, 0, 3), c("blue ", "white", "red"))

# Load the clinical information.
load('Application/TCGA/ClinicalInformation/IndexChange.RData')
# First let's reorder the individual strucutre according to the true subtype.
# there are 7 basal, 4 HER2, 18 Luminal A, 12 Luminal B, 46 NAs no Normals
true_subtype = small_clinic$PAM50.subtype[index_change_vec]
true_subtype = as.character(true_subtype)
for (i in 1:length(true_subtype)){
  if (is.na(true_subtype[i])){
    true_subtype[i] = "Unknown"
  }
}
true_subtype 

basal_index = c(1:87)[true_subtype == "Basal-like"]
HER2_index = c(1:87)[true_subtype == "HER2-enriched"]
LuminalA_index = c(1:87)[true_subtype == "Luminal A"]
LuminalB_index = c(1:87)[true_subtype == "Luminal B"]
Unknown_index = c(1:87)[true_subtype == "Unknown"]
reorder_index = c(basal_index,HER2_index,LuminalA_index,LuminalB_index,Unknown_index)
# there are 7 basal, 4 HER2, 18 Luminal A, 12 Luminal B, 46 unknowns, no Normals.
DMMD_reordered_ind_row1 = result$`Row Decomposition`$`Individual Row 1`[,reorder_index]
DMMD_reordered_ind_row2 = result$`Row Decomposition`$`Individual Row 2`[,reorder_index]
DMMD_reordered_joint_col1 = result$`Column Decomposition`$`Joint Column 1`[,reorder_index]
DMMD_reordered_joint_col2 = result$`Column Decomposition`$`Joint Column 2`[,reorder_index]
DMMD_reordered_ind_col1 = result$`Column Decomposition`$`Individual Column 1`[,reorder_index]
DMMD_reordered_ind_col2 = result$`Column Decomposition`$`Individual Column 2`[,reorder_index]

# Let us create the annotation object first.
ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:6), 
                                        labels = c("Basal","H","LumA","LumB","Unknown"),
                                        labels_gp = gpar(col = "black", fontsize = 10,fontface = "bold")))
split = c(rep(1,7),rep(2,4),rep(3,18),rep(4,12),rep(5,46))

# Plot the heatmaps for cancer data
# Let's first look at the individual row structure for cancer data. 
# Remember
# The rank estimation:
# DMMD result on the dataset -- all equal variance assumption:
# r1 = 8, r2 = 6, r_c = 2, r_r = 0
# JIVE column
# joint rank: 2 , individual ranks: 12 9 
# JIVE row
# joint rank 0, individual rank 11 9

# Let's first look at inidividual row for cancer data
Heatmap(DMMD_reordered_ind_row1, column_split = split, top_annotation = ha, col = hmcol, cluster_rows = T, cluster_columns = FALSE, name = "Individual Signal", row_title = "miRNA", column_title = "Primary Tumor Tissue", show_row_names = FALSE, show_column_names = FALSE) 

# Create a function that adds % to each number in a vector
# Input: vec, the numeric singular value vector.
# Output: a string vector that has the percentage information.
AddPercentage <- function(vec,digit = 1){
  new_vec = round(vec^2/sum(vec^2)*100,digit)
  new_vec = format(new_vec,nsmall = digit)
  new_vec = as.character(new_vec)
  for (i in 1:length(vec)){
     new_vec[i] = paste(new_vec[i],"%",sep = "")
  }
  return(new_vec)
}

# Let's check the SVD result of this matrix. Reordering the column should not affect the rank.
svd_result = svd(DMMD_reordered_ind_row1)
plot(svd_result$d)
variation_vec = AddPercentage(svd_result$d)[1:8]
print(variation_vec)
# We are looking at the row space, so the first 8 columns of V matrix are of our interest.
row_space = svd_result$v[,1:8]
colnames(row_space) = as.character(variation_vec)
# max(row_space) -- 0.32
# min(row_space) -- -0.31
# Let's make a heatmap out of it with the true subtypes along with it.
hmcol2 <- colorRamp2(c(-0.4, 0, 0.4), c("blue ", "white", "red"))
Heatmap(t(row_space), column_split = split, top_annotation = ha, col = hmcol2, cluster_rows = FALSE, cluster_columns = FALSE, name = "Individual row signal", row_title = "Variation Percentage", column_title = "Primary Tumor Tissue", show_row_names = TRUE, show_column_names = FALSE) 

pdf(paste(figure.path, "DMMD_Ind_Row_Cancer.pdf", sep = ""), onefile = T)
Heatmap(t(row_space), column_split = split, top_annotation = ha, col = hmcol2, cluster_rows = FALSE, cluster_columns = FALSE, name = "Individual row signal", row_title = "Variation Percentage", column_title = "Primary Tumor Tissue", show_row_names = TRUE, show_column_names = FALSE) 
dev.off()
# The plot seems good, let's use the same format and plot the rest.
# Individual row for normal data:
svd_result2 = svd(DMMD_reordered_ind_row2)
plot(svd_result2$d)
variation_vec2 = AddPercentage(svd_result2$d)[1:6]
print(variation_vec2)
# We are looking at the row space, so the first 8 columns of V matrix are of our interest.
row_space2 = svd_result2$v[,1:6]
colnames(row_space2) = as.character(variation_vec2)
# max(row_space2) -- 0.27
# min(row_space2) -- -0.26
# Let's make a heatmap out of it with the true subtypes along with it.
Heatmap(t(row_space2), column_split = split, top_annotation = ha, col = hmcol2, cluster_rows = FALSE, cluster_columns = FALSE, name = "Individual row signal", row_title = "Variation Percentage", column_title = "Normal Tissue", show_row_names = TRUE, show_column_names = FALSE) 

pdf(paste(figure.path, "DMMD_Ind_Row_Normal.pdf", sep = ""), onefile = T)
Heatmap(t(row_space2), column_split = split, top_annotation = ha, col = hmcol2, cluster_rows = FALSE, cluster_columns = FALSE, name = "Individual row signal", row_title = "Variation Percentage", column_title = "Normal Tissue", show_row_names = TRUE, show_column_names = FALSE) 
dev.off()

# Let's plot a heatmap with the most important vector.
# We multiply the first vector by -1, which will not affect the variation.
# Because we want to make it consistent with what we find in the individual structure. 
individual_vec1_cancer = -t(row_space[,1])
rownames(individual_vec1_cancer) = colnames(row_space)[1]
Heatmap(individual_vec1_cancer, column_split = split, top_annotation = ha, col = hmcol2, cluster_rows = FALSE, cluster_columns = FALSE, row_title = "Variation Percentage", column_title = "Primary Tumor Tissue", show_row_names = TRUE, show_column_names = FALSE,
        height = 1,width = 7, name = "Individual signal") 

pdf(paste(figure.path, "DMMD_Ind_Row_Cancer_One.pdf", sep = ""), onefile = T, height = 1)
Heatmap(individual_vec1_cancer, column_split = split, top_annotation = ha, col = hmcol2, cluster_rows = FALSE, row_names_side = "left", cluster_columns = FALSE, column_title = "Primary Tumor Tissue", show_row_names = TRUE, show_column_names = FALSE, name = "Individual signal",row_title_gp = gpar(fontsize = 5)) 
dev.off()

individual_vec1_normal = -t(row_space2[,1])
rownames(individual_vec1_normal) = colnames(row_space2)[1]

pdf(paste(figure.path, "DMMD_Ind_Row_Normal_One.pdf", sep = ""), onefile = T, height = 1)
Heatmap(individual_vec1_normal, column_split = split, top_annotation = ha, col = hmcol2, cluster_rows = FALSE, row_names_side = "left", cluster_columns = FALSE, column_title = "Normal Tissue", show_row_names = TRUE, show_column_names = FALSE, name = "Individual signal",row_title_gp = gpar(fontsize = 5)) 
dev.off()
