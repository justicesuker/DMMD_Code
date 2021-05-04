# This is the file that I use to download TCGA_BRCA data mentioned in my draft paper. 

# packages <- c("HGNChelper", "httr", "RCurl", "rjson", "stringr")
# install.packages(packages, dependencies = T)
rm(list = ls())
source("Data/TCGA-BRCA/TCGA_Assembly/miRNA/Module_A.R")
sCancer <- "BRCA"

path1 <- "Data/TCGA-BRCA/TCGA_Assembly/miRNA/UseInDraft"

# Download the miRNA data
miRNAseq_Raw_Cancer <- DownloadmiRNASeqData(cancerType = sCancer, tissueType = 'TP', 
                                            saveFolderName = path1, outputFileName = "miRNA_Cancer")


miRNAseq_Raw_Normal <- DownloadmiRNASeqData(cancerType = sCancer, tissueType = 'NT', 
                                            saveFolderName = path1, outputFileName = "miRNA_Normal")
