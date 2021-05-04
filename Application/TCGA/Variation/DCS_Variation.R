rm(list = ls())
# Load the result of DMMD and JIVE
load('Application/TCGA/2021New/DoubleSTD_DMMD_JIVE_EqualVaraince_Result.RData')
function_path = "Simulations/MyFunction/"
source(paste(function_path,"Preliminary_Functions.R",sep=''))

library(ggplot2)
library(scales)
set.seed(37)

cdata = dcs_cancer
ndata = dcs_normal

# There are 734 genes and 87 samples.

# !!!Remember: cancer data first, normal data second!!!
# Remember: row space -- miRNA gene signal
#           column space -- subject signal

# Check if JIVE gets the total rank the same using different ways:
row_rank_approx_jive = result_row$rankJ + result_row$rankA
print(row_rank_approx_jive)

col_rank_approx_jive = result_col$rankJ + result_col$rankA
print(col_rank_approx_jive)

# Transpose the jive result so that the dimension matches

for (i in 1:2){
  result_col$joint[[i]] = t(result_col$joint[[i]])
  result_col$individual[[i]] = t(result_col$individual[[i]])
}

# Make the colnames and rownames consistent
names_r = rownames(cdata)
names_c = substr(colnames(cdata),1,12)

for (i in 1:2){
  rownames(result_col$joint[[i]]) = names_r
  rownames(result_col$individual[[i]]) = names_r
  rownames(result_row$joint[[i]]) = names_r
  rownames(result_row$individual[[i]]) = names_r
  colnames(result_col$joint[[i]]) = names_c
  colnames(result_col$individual[[i]]) = names_c
  colnames(result_row$joint[[i]]) = names_c
  colnames(result_row$individual[[i]]) = names_c
}

DMMD_J_c1 = result$`Column Decomposition`$`Joint Column 1`
DMMD_J_r1 = result$`Row Decomposition`$`Joint Row 1`
DMMD_I_c1 = result$`Column Decomposition`$`Individual Column 1`
DMMD_I_r1 = result$`Row Decomposition`$`Individual Row 1`
DMMD_Signal1 = cdata - result$Error$Error1
DMMD_Error1 = result$Error$Error1

DMMD_J_c2 = result$`Column Decomposition`$`Joint Column 2`
DMMD_J_r2 = result$`Row Decomposition`$`Joint Row 2`
DMMD_I_c2 = result$`Column Decomposition`$`Individual Column 2`
DMMD_I_r2 = result$`Row Decomposition`$`Individual Row 2`
DMMD_Signal2 = ndata - result$Error$Error2
DMMD_Error2 = result$Error$Error2


JIVE_J_c1 = result_col$joint[[1]]
JIVE_J_r1 = result_row$joint[[1]]
JIVE_I_c1 = result_col$individual[[1]]
JIVE_I_r1 = result_row$individual[[1]]
JIVE_Signal_c1 = JIVE_J_c1 + JIVE_I_c1
JIVE_Signal_r1 = JIVE_J_r1 + JIVE_I_r1
JIVE_Error_c1 = cdata - JIVE_Signal_c1
JIVE_Error_r1 = cdata - JIVE_Signal_r1

JIVE_J_c2 = result_col$joint[[2]]
JIVE_J_r2 = result_row$joint[[2]]
JIVE_I_c2 = result_col$individual[[2]]
JIVE_I_r2 = result_row$individual[[2]]
JIVE_Signal_c2 = JIVE_J_c2 + JIVE_I_c2
JIVE_Signal_r2 = JIVE_J_r2 + JIVE_I_r2
JIVE_Error_c2 = ndata - JIVE_Signal_c2
JIVE_Error_r2 = ndata - JIVE_Signal_r2

# Check if joint and individual are perpendicular by checking the Fnorm^2
Fnorm(DMMD_Signal1)^2
Fnorm(DMMD_J_c1)^2 + Fnorm(DMMD_I_c1)^2
Fnorm(DMMD_J_r1)^2 + Fnorm(DMMD_I_r1)^2

Fnorm(JIVE_Signal_c1)^2
Fnorm(JIVE_J_c1)^2 + Fnorm(JIVE_I_c1)^2

variation_frame1 = as.data.frame(matrix(rep(NA, 36) ,nrow = 12, ncol = 3))
colnames(variation_frame1) = c("SSE", "Section", "Method")
variation_frame1[,"Section"] = rep(c("Joint","Individual","Error"),4)

variation_frame1[,"Method"] = c(rep("DMMD (subjects)",3),rep("DMMD (miRNAs)",3),
                               rep("JIVE (subjects)",3),rep("JIVE (miRNAs)",3))
variation_frame1[,1] = c(Fnorm(DMMD_J_r1)^2,Fnorm(DMMD_I_r1)^2,Fnorm(DMMD_Error1)^2,
                         Fnorm(DMMD_J_c1)^2,Fnorm(DMMD_I_c1)^2,Fnorm(DMMD_Error1)^2,
                         Fnorm(JIVE_J_r1)^2,Fnorm(JIVE_I_r1)^2,Fnorm(JIVE_Error_r1)^2,
                         Fnorm(JIVE_J_c1)^2,Fnorm(JIVE_I_c1)^2,Fnorm(JIVE_Error_c1)^2)
  

variation_frame2 = as.data.frame(matrix(rep(NA, 36) ,nrow = 12, ncol = 3))
colnames(variation_frame2) = c("SSE", "Section", "Method")
variation_frame2[,"Section"] = rep(c("Joint","Individual","Error"),4)

variation_frame2[,"Method"] = c(rep("DMMD (subjects)",3),rep("DMMD (miRNAs)",3),
                                rep("JIVE (subjects)",3),rep("JIVE (miRNAs)",3))
variation_frame2[,1] = c(Fnorm(DMMD_J_r2)^2,Fnorm(DMMD_I_r2)^2,Fnorm(DMMD_Error2)^2,
                         Fnorm(DMMD_J_c2)^2,Fnorm(DMMD_I_c2)^2,Fnorm(DMMD_Error2)^2,
                         Fnorm(JIVE_J_r2)^2,Fnorm(JIVE_I_r2)^2,Fnorm(JIVE_Error_r2)^2,
                         Fnorm(JIVE_J_c2)^2,Fnorm(JIVE_I_c2)^2,Fnorm(JIVE_Error_c2)^2)

figure.path = 'Application/TCGA/Variation/'

pdf(paste(figure.path, "VariationExplained_Cancer.pdf", sep = ""), onefile = T)
s1 <- ggplot(variation_frame1, aes(Method, SSE, fill = Section))
s1 <- s1 + geom_bar(position = "fill", stat = "identity") + scale_y_continuous(labels = percent) +
  ggtitle("Primary Tumor Tissue") + ylab("Percentage Variance Explained") + theme(axis.text = element_text(size = 10),legend.text = element_text(size = 10), axis.title = element_text(size = 12))
print(s1)
dev.off()

pdf(paste(figure.path, "VariationExplained_Normal.pdf", sep = ""), onefile = T)
s2 <- ggplot(variation_frame2, aes(Method, SSE, fill = Section))
s2 <- s2 + geom_bar(position = "fill", stat = "identity") + scale_y_continuous(labels = percent) + 
  ggtitle("Normal Tissue") + ylab("Percentage Variance Explained") + theme(axis.text = element_text(size = 10),legend.text = element_text(size = 10), axis.title = element_text(size = 12))
print(s2)
dev.off()



