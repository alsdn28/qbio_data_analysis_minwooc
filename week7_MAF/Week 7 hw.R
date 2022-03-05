getwd()
setwd("/Users/minwoocho/Desktop/qbio_data_analysis_minwooc/analysis_data/")
BiocManager::install("maftools")
library(TCGAbiolinks)
library(maftools)

clinic <- data.table::fread("/Users/minwoocho/Desktop/qbio_data_analysis_minwooc/coad_clinical_data.csv", data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)
getwd()
setwd("/Users/minwoocho/Desktop/qbio_data_analysis_minwooc/analysis_data/GDCdata/")

maf_dataframe <- data.table::fread("./GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", data.table = F)
clinic <- data.table::fread("../coad_clinical_data.csv",
                            data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)
#2 explore MAF
maf_object
str(maf_object)
maf_object@data

# MAF analysis
oncoplot(maf = maf_object,
         top = 10)
# APC is the most mutated which is a tumor supressor protein
clinic <- maf_object@clinical.data
clinic$age_category <- ifelse(clinic$age_at_initial_pathologic_diagnosis <= 50, "Young", "Old")

name <- ifelse(clinic$age_category=="Young", TRUE, FALSE)
young_patient_ids <- clinic[name, Tumor_Sample_Barcode]
young_maf = subsetMaf(maf = young_maf,
                      tsb = young_patient_ids)

name2 <- ifelse(clinic$age_category=="Old", TRUE, FALSE)
old_patient_ids <- clinic[name2, Tumor_Sample_Barcode]
old_maf = subsetMaf(maf = old_maf,
                      tsb = old_patient_ids)

coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Young Patients", 
           m2Name = "Old Patients")
# less % of genes mutated in younger patients

#3.3
lollipopPlot(maf_object, gene = "APC")

lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name = "Young Patients",
              m2_name = "Old Patients",
              gene = "APC")
  # nonsense mutation most common, genes less mutated in younger patients

#3.9
# b = 7, c = 2, d = 35, e = 37, f = 42
# A and B likely independent

geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")
geneB_maf <- subsetMaf(maf = maf_object,
                       genes = "KRAS")
#3.11 
# subsetMaf() subsetted the maf object dataframe into data containing certain genes
# there is more than one mutation in gene A because data dimension is much greater in maf_object 

mut_bc_geneA <- c(geneA_maf@clinical.data$Tumor_Sample_Barcode)
num_mut_geneA <- length(mut_bc_geneA)

mut_bc_geneB <- c(geneB_maf@clinical.data$Tumor_Sample_Barcode)
num_mut_geneB <- length(mut_bc_geneB)

mut_bc_geneAB <- intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB <- length(mut_bc_geneAB)
#78 patients with mutations in both genes

num_mut_geneA_only <- num_mut_geneA - num_mut_geneAB
num_mut_geneB_only <- num_mut_geneB - num_mut_geneAB
num_neither_mutation <- length(maf_object@clinical.data$Tumor_Sample_Barcode) - num_mut_geneA - num_mut_geneB + num_mut_geneAB

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_neither_mutation), 
                       nrow=2)
contig_table

fe_results <- fisher.test(contig_table)
fe_results


















