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
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Load the result of DMMD and JIVE
load('Application/TCGA/2021New/DoubleSTD_DMMD_JIVE_EqualVaraince_Result.RData')

set.seed(37)

cdata = dcs_cancer
ndata = dcs_normal

# There are 734 miRNAs and 87 samples.

# !!!Remember: cancer data first, normal data second!!!
# Remember: row space -- miRNA miRNA signal
#           column space -- subject signal

# The rank estimation:
# DMMD result on the dataset -- all equal variance assumption:
# r1 = 8, r2 = 6, r_c = 2, r_r = 0
# JIVE column
# joint rank: 2 , individual ranks: 12 9 
# JIVE row
# joint rank 0, individual rank 11 9

# Specify path to save figures
figure.path = "Application/TCGA/2021New/JointmiRNA/"

hmcol <- colorRamp2(c(-2.5, 0, 2.5), c("blue ", "white", "red"))

# DMMD
# The order of rows and columns is based on the joint column structure 1 of DMMD
# The ordering for cancer data and normal data is the same so we could plot in a single file.
dist_miRNA_DMMD <- dist(result$`Column Decomposition`$`Joint Column 1`, method = "euclidean", diag = FALSE, upper = FALSE)
clus_miRNA_DMMD <- hclust(dist_miRNA_DMMD, method = "ward.D2", members = NULL)
dist_sample_DMMD <- dist(t(result$`Column Decomposition`$`Joint Column 1`), method = "euclidean", diag = FALSE, upper = FALSE)
clus_sample_DMMD <- hclust(dist_sample_DMMD, method = "ward.D2", members = NULL)

pdf(paste(figure.path, "DMMD_Joint_Col_Reordered_Tumor.pdf", sep = ""), onefile = T)
Heatmap(result$`Column Decomposition`$`Joint Column 1`, col = hmcol, cluster_rows = clus_miRNA_DMMD, cluster_columns = clus_sample_DMMD, name = "Joint Signal", row_title = "miRNA", column_title = "Primary Tumor Tissue", show_row_names = FALSE, show_column_names = FALSE)
dev.off()

pdf(paste(figure.path, "DMMD_Joint_Col_Reordered_Normal.pdf", sep = ""), onefile = T)
Heatmap(result$`Column Decomposition`$`Joint Column 2`, col = hmcol, cluster_rows = clus_miRNA_DMMD, cluster_columns = clus_sample_DMMD, name = "Joint Signal", row_title = "miRNA", column_title = "Normal Tissue", show_row_names = FALSE, show_column_names = FALSE) 
dev.off()
# If we want to have n by p matrices as described in draft, we need to plot:
pdf(paste(figure.path, "DMMD_Joint_Col_Reordered_Tumor_n_by_p.pdf", sep = ""), onefile = T)
Heatmap(t(result$`Column Decomposition`$`Joint Column 1`), col = hmcol, cluster_rows = clus_sample_DMMD, 
        cluster_columns = clus_miRNA_DMMD, row_title = "Primary Tumor Tissue", column_title = "miRNA", 
        show_row_names = FALSE, show_column_names = FALSE, name = "Joint Signal")
dev.off()

pdf(paste(figure.path, "DMMD_Joint_Col_Reordered_Normal_n_by_p.pdf", sep = ""), onefile = T)
Heatmap(t(result$`Column Decomposition`$`Joint Column 2`), col = hmcol, cluster_rows = clus_sample_DMMD, 
        cluster_columns = clus_miRNA_DMMD, row_title = "Normal Tissue", column_title = "miRNA", 
        show_row_names = FALSE, show_column_names = FALSE, name = "Joint Signal") 
dev.off()