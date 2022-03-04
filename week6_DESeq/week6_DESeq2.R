getwd()
setwd("/Users/minwoocho/Desktop/qbio_studentresearch_minwooc/analysis_data/")
BiocManager::install("SummarizedExperiment")
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)
BiocManager::install("DESeq2")
library(DESeq2)
counts <- assays(sum_exp)$"HTSeq - Counts"
patient_data <- colData(sum_exp)

age_na <- (is.na(patient_data$age_at_index))
patient_data <- patient_data[!age_na, ]
counts <- counts[ , !age_na]





patient_data$age_category <- ifelse(patient_data$age_at_index<50, "young", "old")
patient_data$age_category = factor(patient_data$age_category, levels = c("young", "old"))

if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
counts_row_sums <- rowSums(counts)
low_counts_mask <- ifelse(counts_row_sums>=10, TRUE, FALSE)
sum(low_counts_mask)
counts <- counts[low_counts_mask, ]

#2
load("DESeq2")
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patient_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

# Analysis

my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))

order_indices = order(my_df$y)

row_order <- order(results$padj)

log2FoldChange_threshold <- 1
padj_threshold <- 0.05

log_mask <- ifelse(results$log2FoldChange > log2FoldChange_threshold | 
                     results$log2FoldChange < -log2FoldChange_threshold, TRUE, FALSE)
sum(log_mask)
padj_no_na <- is.na(results$padj)
results <- results[!padj_no_na, ]

padj_mask <- ifelse(results$padj < padj_threshold, TRUE, FALSE)
sum(padj_mask)

# Volcano Plots
fc_threshold = 2 
p_threshold = 0.05

plot(x = log2FoldChange_threshold,
     y = -log10(padj),
     xlab = "young/old", # be sure the specify that it's young over old!
     ylab = "p-value",
     pch = 20) # smaller solid circles

library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in young",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")


volcano_plot




