#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
setwd("/Users/minwoocho/Desktop/qbio_data_analysis_minwooc/analysis_data/")
clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
GDCdownload(clin_query)
#downloaded GDC 
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"
str(clinic)
colnames(clinic)
clinic$bcr_patient_barcode
plot(clinic$age_at_initial_pathologic_diagnosis, clinic$weight, xlab="age", ylab="weight")
boxplot(clinic$age_at_initial_pathologic_diagnosis~clinic$race_list)
gsub("^$", "No data", clinic$race_list)
#2.2
min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)
young_old <- ifelse(clinic$age_at_initial_pathologic_diagnosis>=50, TRUE, FALSE)
sum(young_old)
#2.6
young_patient_ids <- which(young_old==FALSE)
young_patient_ids <- clinic[young_patient_ids, 14]
old_patient_ids <- which(young_old==TRUE)
old_patient_ids <- clinic[old_patient_ids, 14]
#2.7
clinic$age_category <- ifelse(young_old, TRUE, FALSE)
#2.8
#leaving a blank means it accounts all of the column or row 
clinic[1,1]
#2.9
young_clinic <- clinic[which(clinic$age_category==FALSE), ]
old_clinic <- clinic[which(clinic$age_category==TRUE), ]
young_clinic_one_line <- clinic[which(clinic$age_at_initial_pathologic_diagnosis<50), ]

#3
install.packages("survival")
install.packages("survminer")
library(survival)
library(survminer)
library(dyplyr)

         