getwd()
setwd("/Users/minwoocho/Desktop/qbio_data_analysis_minwooc/analysis_data/")
getwd()
BiocManager::install("SummarizedExperiment")
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
GDCdownload(query) # only need to download the data once! Comment this out once you have completed it once
sum_exp <- GDCprepare(query)
str(sum_exp)
assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]
rowData(sum_exp)
head(rowData(sum_exp))
colData(sum_exp)
colData(sum_exp)[1:5, 25:29]
metadata(sum_exp)

sum_exp_copy <- sum_exp

# 2.1
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")

str(colData(sum_exp))
# rows correspond to columns of assays

colnames(colData(sum_exp))
colData(sum_exp)$age_at_diagnosis[1:10]
is.na(sum_exp_copy$age_at_diagnosis) <- 0
sum_exp_copy$age_at_diagnosis <- sum_exp_copy$age_at_diagnosis / 365
sum_exp_copy$age_category <- ifelse(sum_exp_copy$age_at_diagnosis>=50, "Old", "Young")

# 2.2
head(rowData(sum_exp_copy))
dim(rowData(sum_exp_copy))
"KRAS" %in% rowData(sum_exp)$external_gene_name
"FGR" %in% rowData(sum_exp)$external_gene_name

# 2.3
assays(sum_exp)$"HTSeq - Counts"[20:25, 30:35]
geneA_id_mask <- rowData(sum_exp_copy)$external_gene_name == "KRAS"
sum(geneA_id_mask)
ensembl_geneA <- rowData(sum_exp_copy)$ensembl_gene_id[geneA_id_mask]

geneB_id_mask <- rowData(sum_exp_copy)$external_gene_name == "FGR"
sum(geneB_id_mask)
ensembl_geneB <- rowData(sum_exp_copy)$ensembl_gene_id[geneB_id_mask]
# GeneA = ENSG00000000938
# GeneB = ENSG00000133703

min(assays(sum_exp_copy)$"HTSeq - Counts"[ensembl_geneA, ])
max(assays(sum_exp_copy)$"HTSeq - Counts"[ensembl_geneA, ])
summary(assays(sum_exp_copy)$"HTSeq - Counts"[ensembl_geneB, ])

data_for_geneA <- assays(sum_exp_copy)$"HTSeq - Counts"[ensembl_geneA, ]
data_for_geneB <- assays(sum_exp_copy)$"HTSeq - Counts"[ensembl_geneB, ]
plot(data_for_geneA,
     data_for_geneB, 
     xlab = "KRAS",
     ylab = "FGR"
  
)

bool_age_na <- is.na(colData(sum_exp_copy)$age_category)
num_na <- sum(bool_age_na)                     
age_cat_no_NAs = !bool_age_na
length(assays(sum_exp_copy)$"HTSeq - Counts"[age_cat_no_NAs])
sum(age_cat_no_NAs)
sum_exp_copy = age_cat_no_NAs + bool_age_na

identical((rownames(colData(sum_exp_copy))), colnames(assays(sum_exp_copy)$"HTSeq - Counts"))
gene_counts = dataframe[bool_age_na, ensembl_geneA]
boxplot(gene_counts ~ age_cat_no_NAs, xlab = "", ylab = "")