setwd("Data/TCGA-BRCA/TCGA_Assembly/miRNA/UseInDraft")
list.files()
### 8 files:
# "miRNA_Cancer__BRCA__mir_GA.hg18__TP__20210106135139.txt"             
# "miRNA_Cancer__BRCA__mir_GA.hg19.mirbase20__TP__20210106135225.txt"   
# "miRNA_Cancer__BRCA__mir_HiSeq.hg18__TP__20210106135404.txt"          
# "miRNA_Cancer__BRCA__mir_HiSeq.hg19.mirbase20__TP__20210106135546.txt"
# "miRNA_Normal__BRCA__mir_GA.hg18__NT__20210106140355.txt"             
# "miRNA_Normal__BRCA__mir_GA.hg19.mirbase20__NT__20210106140402.txt"   
# "miRNA_Normal__BRCA__mir_HiSeq.hg18__NT__20210106140415.txt"          
# "miRNA_Normal__BRCA__mir_HiSeq.hg19.mirbase20__NT__20210106140430.txt"

# In the end, we find that the dataset of hg18 and hg19 are very similar (we doubt if they are the same).
# We also think HiSeq.hg18 is the same as HiSeq.hg19

#### For cancer data: ####
# miRNA from HiSeq.hg19 - #1
TEMP_Hiseq.hg19_Cancer = read.table("miRNA_Cancer__BRCA__mir_HiSeq.hg19.mirbase20__TP__20210106135546.txt",header = T)
TEMP_Hiseq.hg19_Cancer[1:5,1:5]
colnames(TEMP_Hiseq.hg19_Cancer)
# We found that for each sample, there are two data points for each miRNA related to it.
# One is the raw count and the other is reads_per_million_miRNA_mapped. 
# We decide to use the raw count.

# Get the miRNA ID.
rowname_vec_19_hiseq = TEMP_Hiseq.hg19_Cancer[,TEMP_Hiseq.hg19_Cancer[1,] == "miRNA_ID"]
# Get rid of the first redundent name.
rowname_vec_19_hiseq = as.character(rowname_vec_19_hiseq[-1])
rowname_vec_19_hiseq
# Get the read count data set. Call it Hiseq_hg19_Cancer.
Hiseq_hg19_Cancer = TEMP_Hiseq.hg19_Cancer[,TEMP_Hiseq.hg19_Cancer[1,] == "read_count"]
Hiseq_hg19_Cancer = Hiseq_hg19_Cancer[-1,]
rownames(Hiseq_hg19_Cancer) = rowname_vec_19_hiseq
Hiseq_hg19_Cancer[1:5,1:5]
# Similar data processing procedure is applied to the other data.

# miRNA from HiSeq.hg18 - #2
TEMP_Hiseq.hg18_Cancer = read.table("miRNA_Cancer__BRCA__mir_HiSeq.hg18__TP__20210106135404.txt",header = T)
TEMP_Hiseq.hg18_Cancer[1:5,1:5]
rowname_vec_18_hiseq = TEMP_Hiseq.hg18_Cancer[,TEMP_Hiseq.hg18_Cancer[1,] == "miRNA_ID"]
rowname_vec_18_hiseq = as.character(rowname_vec_18_hiseq[-1])
rowname_vec_18_hiseq[1:5]
# Get the read count data set. Call it Hiseq_hg19_Cancer.
Hiseq_hg18_Cancer = TEMP_Hiseq.hg18_Cancer[,TEMP_Hiseq.hg18_Cancer[1,] == "read_count"]
Hiseq_hg18_Cancer = Hiseq_hg18_Cancer[-1,]
rownames(Hiseq_hg18_Cancer) = rowname_vec_18_hiseq
Hiseq_hg18_Cancer[1:5,1:5]

# miRNA from GA_hg18 - #3
TEMP_GA.hg18_Cancer = read.table("miRNA_Cancer__BRCA__mir_GA.hg18__TP__20210106135139.txt",header = T)
TEMP_GA.hg18_Cancer[1:5,1:5]
rowname_vec18_cancer = TEMP_GA.hg18_Cancer[,TEMP_GA.hg18_Cancer[1,] == "miRNA_ID"]
rowname_vec18_cancer = as.character(rowname_vec18_cancer[-1])
GAhg18_cancer = TEMP_GA.hg18_Cancer[,TEMP_GA.hg18_Cancer[1,] == "read_count"]
GAhg18_cancer = GAhg18_cancer[-1,]
rownames(GAhg18_cancer) = rowname_vec18_cancer
GAhg18_cancer[1:5,1:5]

# miRNA from GA_hg19 - #4
TEMP_GA.hg19_Cancer = read.table("miRNA_Cancer__BRCA__mir_GA.hg19.mirbase20__TP__20210106135225.txt",header = T)
TEMP_GA.hg19_Cancer[1:5,1:5]
rowname_vec19_cancer = TEMP_GA.hg19_Cancer[,TEMP_GA.hg19_Cancer[1,] == "miRNA_ID"]
rowname_vec19_cancer = as.character(rowname_vec19_cancer[-1])
GAhg19_cancer = TEMP_GA.hg19_Cancer[,TEMP_GA.hg19_Cancer[1,] == "read_count"]
GAhg19_cancer = GAhg19_cancer[-1,]
rownames(GAhg19_cancer) = rowname_vec19_cancer
GAhg19_cancer[1:5,1:5]


#### For normal data: ####
# miRNA from HiSeq.hg18 - #5
TEMP_Hiseq.hg18_Normal = read.table("miRNA_Normal__BRCA__mir_HiSeq.hg18__NT__20210106140415.txt",header = T)
TEMP_Hiseq.hg18_Normal[1:5,1:5]
rowname_vec_18_hiseq_normal = TEMP_Hiseq.hg18_Normal[,TEMP_Hiseq.hg18_Normal[1,] == "miRNA_ID"]
rowname_vec_18_hiseq_normal = as.character(rowname_vec_18_hiseq_normal[-1])
Hiseq_hg18_Normal = TEMP_Hiseq.hg18_Normal[,TEMP_Hiseq.hg18_Normal[1,] == "read_count"]
Hiseq_hg18_Normal = Hiseq_hg18_Normal[-1,]
rownames(Hiseq_hg18_Normal) = rowname_vec_18_hiseq_normal
Hiseq_hg18_Normal[1:5,1:5]

# miRNA from HiSeq.hg19 - #6
TEMP_Hiseq.hg19_Normal = read.table("miRNA_Normal__BRCA__mir_HiSeq.hg19.mirbase20__NT__20210106140430.txt",header = T)
TEMP_Hiseq.hg19_Normal[1:5,1:5]
rowname_vec_19_hiseq_normal = TEMP_Hiseq.hg19_Normal[,TEMP_Hiseq.hg19_Normal[1,] == "miRNA_ID"]
rowname_vec_19_hiseq_normal = as.character(rowname_vec_19_hiseq_normal[-1])
Hiseq_hg19_Normal = TEMP_Hiseq.hg19_Normal[,TEMP_Hiseq.hg19_Normal[1,] == "read_count"]
Hiseq_hg19_Normal = Hiseq_hg19_Normal[-1,]
rownames(Hiseq_hg19_Normal) = rowname_vec_19_hiseq_normal
Hiseq_hg19_Normal[1:5,1:5]

# miRNA from GA.hg18 - #7
TEMP_GA.hg18_Normal = read.table("miRNA_Normal__BRCA__mir_GA.hg18__NT__20210106140355.txt",header = T)
TEMP_GA.hg18_Normal[1:5,1:5]
rowname_vec18_normal = TEMP_GA.hg18_Normal[,TEMP_GA.hg18_Normal[1,] == "miRNA_ID"]
rowname_vec18_normal = as.character(rowname_vec18_normal[-1])
GAhg18_normal = TEMP_GA.hg18_Normal[,TEMP_GA.hg18_Normal[1,] == "read_count"]
GAhg18_normal = GAhg18_normal[-1,]
rownames(GAhg18_normal) = rowname_vec18_normal
GAhg18_normal[1:5,1:5]

# miRNA from GA.hg19 - #8
TEMP_GA.hg19_Normal = read.table("miRNA_Normal__BRCA__mir_GA.hg19.mirbase20__NT__20210106140402.txt",header = T)
TEMP_GA.hg19_Normal[1:5,1:5]
rowname_vec19_normal = TEMP_GA.hg19_Normal[,TEMP_GA.hg19_Normal[1,] == "miRNA_ID"]
rowname_vec19_normal = as.character(rowname_vec19_normal[-1])
GAhg19_normal = TEMP_GA.hg19_Normal[,TEMP_GA.hg19_Normal[1,] == "read_count"]
GAhg19_normal = GAhg19_normal[-1,]
rownames(GAhg19_normal) = rowname_vec19_normal
GAhg19_normal[1:5,1:5]

save(list = c("GAhg18_cancer","GAhg19_cancer","Hiseq_hg18_Cancer","Hiseq_hg19_Cancer",
              "GAhg18_normal","GAhg19_normal","Hiseq_hg18_Normal","Hiseq_hg19_Normal"),
     file = "../New_TCGA.RData")
