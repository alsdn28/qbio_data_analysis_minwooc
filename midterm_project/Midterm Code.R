setwd("/Users/minwoocho/Desktop/qbio_data_analysis_minwooc/analysis_data/")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)

maf_dataframe <- data.table::fread("./GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)
# maf gender analysis
clinic <- maf_object@clinical.data

name <- ifelse(clinic$gender=="MALE", TRUE, FALSE)
male_patient_ids <- clinic[name, Tumor_Sample_Barcode]
male_maf <- subsetMaf(maf = maf_object,
                      tsb = male_patient_ids)

name2 <- ifelse(clinic$gender=="FEMALE", TRUE, FALSE)
female_patient_ids <- clinic[name2, Tumor_Sample_Barcode]
female_maf = subsetMaf(maf = maf_object,
                    tsb = female_patient_ids)

coOncoplot(m1 = male_maf, 
           m2 = female_maf, 
           m1Name = "Male Patients", 
           m2Name = "Female Patients")

# over/under expression of gene y in gender
library(DESeq2)

counts <- assays(sum_exp)$"HTSeq - Counts"
patient_data <- colData(sum_exp)

gender_na <- (is.na(patient_data$gender))
patient_data <- patient_data[!gender_na, ]
counts <- counts[ , !gender_na]

patient_data$gender = factor(patient_data$gender, levels = c("male", "female"))

if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
  
counts_row_sums <- rowSums(counts)
low_counts_mask <- ifelse(counts_row_sums>=10, TRUE, FALSE)
sum(low_counts_mask)
counts <- counts[low_counts_mask, ]

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patient_data,
                             design = ~gender)

dds_obj = DESeq(dds)
resultsNames(dds_obj) 

results = results(dds_obj, format = "DataFrame", contrast = c("gender", "male", "female"))

row_order <- order(results$padj)
results <- results[row_order, ]

log2FoldChange_threshold <- 1
padj_threshold <- 0.05

log2_mask <- ifelse(results$log2FoldChange > log2FoldChange_threshold | 
  results$log2FoldChange < -log2FoldChange_threshold, TRUE, FALSE)

padj_mask <- results$padj < padj_threshold

results <- results[log2_mask&&padj_mask, ]

fc_threshold = 2  
p_threshold = 0.05  


plot(x = -log10(results$padj),
     y = results$log2FoldChange,
     xlab = "male/female", # be sure the specify that it's young over old!
     ylab = "log2foldchange",
     pch = 20) # smaller solid circles
abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05


library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in males",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in males", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (male/female)",
       y = "-log10 Adjusted p-value")


volcano_plot

# KM curve
library(survival)
library(survminer)
clinic_copy <- clinic
clinic_copy$days_to_death <- ifelse(is.na(clinic_copy$days_to_death), clinic_copy$days_to_last_follow_up, clinic_copy$days_to_death)

clinic_copy$survival_time <- as.numeric(clinic_copy$days_to_death)

clinic_copy$death_event = as.integer(clinic$vital_status == "Dead")


surv_object <- Surv(time = clinic_copy$survival_time, 
                    event = clinic_copy$death_event)
gender_fit <- surv_fit( surv_object ~ clinic_copy$gender, data = clinic_copy )

survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
survplot
