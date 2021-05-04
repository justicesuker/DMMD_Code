library(stringr)
rm(list=ls())
load("Data/TCGA-BRCA/TCGA_Assembly/miRNA/New_TCGA.RData")
dim(Hiseq_hg18_Cancer)
dim(Hiseq_hg18_Normal)
head(Hiseq_hg18_Cancer)
# Check if the miRNA name matches.
# We use the assay from Illumina Hiseq Assay on Hg18 reference genome. 
# Data is raw count of miRNA. See https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/miRNA_Pipeline/
# For more information, check manual of the TCGA-Assemble package.

identical(rownames(Hiseq_hg18_Cancer),rownames(Hiseq_hg18_Normal))
p1 = dim(Hiseq_hg18_Cancer)[2]
p2 = dim(Hiseq_hg18_Normal)[2]

# We need to make sure the samples match, too.
# Since the number of samples in normal data (87) is greatly less than the number of samples in cancer data (772),
# We only need to find the mathced sample in cancer data, so we could form a duoble-matched matrices with size 1046*87.
# Check the pattern of the sample name
colnames(Hiseq_hg18_Cancer)[1:5]
colnames(Hiseq_hg18_Normal)[1:5]
# The first 12 strings are the patient ID. e.g. TCGA.3C.AAAU
# TCGA_SampleTypeID = 01 stands for "Primary Solid Type" which is corresponding to tissueType = TP 
# TCGA_SampleTypeID = 11 stands for "Solid Tissue Normal" which is corresponding to tissueType = NT
# The (e.g.) 11R after 01A is the portion, which represents the order of portion in a sequence of 100-120 mg sample portions
# The R in 11R represents the molecular type of analyte for analysis. 
# For more information, go to https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

# We simply focus on the first 12 letters and make sure they match.
substr("TCGA.3C.AAAU.01A.11R.A41G.13",1,12)
str_sub("TCGA.3C.AAAU.01A.11R.A41G.13",-7,-1)
cancer_index = c()
normal_index = c()
for (i in 1:p2){
  for (j in 1:p1){
    if (substr(colnames(Hiseq_hg18_Cancer)[j],1,12) == substr(colnames(Hiseq_hg18_Normal)[i],1,12) &&
      str_sub(colnames(Hiseq_hg18_Cancer)[j],-7,-1) == str_sub(colnames(Hiseq_hg18_Normal)[i],-7,-1))
      break
  }
  cancer_index = append(cancer_index,j)
  normal_index = append(normal_index,i)
}

temp_cancer_data = Hiseq_hg18_Cancer[,cancer_index]
temp_normal_data = Hiseq_hg18_Normal[,normal_index]

for (i in 1:length(normal_index)){
  temp_cancer_data[,i] = as.numeric(as.character(temp_cancer_data[,i]))
  temp_normal_data[,i] = as.numeric(as.character(temp_normal_data[,i]))
}

# Take the log-transformation.
cancer_data = log(temp_cancer_data + 1)
normal_data = log(temp_normal_data + 1)

save(cancer_data, normal_data, file = "TCGA_BRCA_New_Data.RData")

