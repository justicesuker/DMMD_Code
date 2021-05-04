list.files()
clinical = read.table(file = 'Application/TCGA/ClinicalInformation/brca_tcga_pub_clinical_data.tsv', sep = '\t', header = TRUE)
# Explore basic information
head(clinical)
colnames(clinical)
# The following information is the most improtant to us:
# Patient.ID and PAM50.subtype.
# Except that, we might find some relationship between sex and developing cancer?
interest_col = c("Patient.ID","PAM50.subtype","Sex")
small_clinic = clinical[,interest_col]
head(small_clinic)
save(small_clinic,file = "Application/TCGA/ClinicalInformation/Subtype_Sex_Information.RData")
