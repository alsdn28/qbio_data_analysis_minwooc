getwd()
setwd("/Users/minwoocho/Desktop/qbio_data_analysis_minwooc/")
library(TCGAbiolinks)
clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"
clinic_copy <- clinic

# 1 
# continuous variable: quantitative values such as integers, decimals, etc
# categorical variable: stored as factor, qualitative data, gender types etc
# discrete variable: quantitative values but is finite, number of people present

# 2
# chose days_to_birth, 2 na values
na_age <- is.na(clinic_copy$days_to_birth)
sum(na_age)

# 3
# age variable is continuous

# 4 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6788938/
# colorectal cancer in adults over 50 increase in US due to epidemiology 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6885269/
# age disparaties in colorectal cancer, especially in young adults

# 5 
# chose whether kras mutations found. alternations in the kras gene classified as mutations into a categorical variable
colnames(clinic_copy)
na_mutations <- is.na(clinic_copy$kras_mutation_codon)
sum(na_mutations)

# 6 
# 1. the greater age, the more chance that mutations may have happened, but unlikely 
# 2. the greater the age, survival rate with colorectal cancer is likely to be low
# 3. if there is a mutation, surivival rate is also low as colorectal cancer will develop

ifelse(clinic_copy$age_at_initial_pathologic_diagnosis, "Young", "Old")
is.na(clinic_copy$kras_mutation_found)
na_loci <- is.na(clinic_copy$number_of_abnormal_loci)
sum(na_loci)
mean_age <- mean(clinic_copy$age_at_initial_pathologic_diagnosis)
clinic_copy$age_category <- ifelse(clinic_copy$age_at_initial_pathologic_diagnosis>=67.83, TRUE, FALSE)
sum(clinic_copy$age_category)
plot(clinic_copy$age_at_initial_pathologic_diagnosis, clinic_copy$number_of_abnormal_loci, 
     xlab="age", ylab="abnormal loci",
     col = ifelse(clinic_copy$age_at_initial_pathologic_diagnosis<mean_age, 'red', 'green'))

# coding 
library(survival)
library(survminer)
na_daystodeath <- is.na(clinic_copy$days_to_death)
clinic_copy$days_to_death[is.na(clinic_copy$days_to_death)] <- clinic_copy$days_to_last_followup[is.na(clinic_copy$days_to_death)]
clinic_copy$survival_time <- clinic_copy$days_to_death - clinic_copy$days_to_birth
clinic_copy$death_event <- ifelse(clinic$vital_status=="Alive", 1, 0)

surv_object <- Surv(time = clinic_copy$days_to_death, 
                    event = clinic_copy$death_event)
race_fit <- surv_fit( surv_object ~ clinic_copy$race_list, data = clinic_copy )

survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
