#### Overview ####

## correlation test between the selected CpGs and to age as well

#### Preparation ####
library(tibble)
library(dplyr)
library(ggplot2)

#### loading data

## methyl values
load("/workspace/cfngle/results-data/05_shared_methyl_values/AC_meth_values_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/AS_meth_values_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/EH_meth_values_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/JM_meth_values_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/ZF_meth_values_JM.Rdata")

## methyl sites
load("/workspace/cfngle/results-data/05_shared_methyl_values/AC_methyl_sites_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/AC_selected_CpGs")
AC_methyl_sites_JM <- AC_methyl_sites_JM[AC_selected_CpGs,]
load("/workspace/cfngle/results-data/05_shared_methyl_values/AS_methyl_sites_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/EH_methyl_sites_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/JM_methyl_sites_JM.Rdata")
load("/workspace/cfngle/results-data/05_shared_methyl_values/ZF_methyl_sites_JM.Rdata")

### methyl values for all samples including age
## AC
xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/meth-corrected-batchcorrected-cod.Rdata")
assign("AC_meth_data", get(xx))
AC_meth_data <- as.data.frame(AC_meth_data)
AC_age <- AC_meth_data$age

##AS
xx <- load("/workspace/cfngle/raw-data/AS/zzz_methyl_data/Meth-complete-snapper.RData")
assign("AS_meth_data", get(xx))
AS_age <- AS_meth_data$age

##EH
xx <- load("/workspace/cfngle/raw-data/EH/zzz-methyl_data/Meth-complete-hake.RData")
assign("EH_meth_data", get(xx))
EH_age <- EH_meth_data$age
EH_metadata_samples <- read.csv("/workspace/cfngle/raw-data/EH/zzz-methyl_data/hake-samples.txt", sep = "\t")
EH_sex <- EH_metadata_samples$sex

##JM
JM_meth_data <- load("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_243285_CpGs.Rdata")
JM_meth_data <- JM_24_methyl_data
JM_age <- JM_meth_data$age

##ZF
ZF_meth_data <- load("/workspace/cfngle/raw-data/ZF/zzz_methyldata/ZF_methyldata_88.RData")
ZF_meth_data <- ZF_methyl_data
ZF_age <- ZF_meth_data$age


#### functions ####

cor.test.age <- function(methyl_values, age, SMR = "not_defined", species = "undefined", method = "pearson") {
  correlation_results <- list()
  print(paste0("Running correlation test against age with ", method, " method. Results are stored in tibble."))
  # Loop through each methylation site
  for (i in 1:ncol(methyl_values)) {
    site_name <- colnames(methyl_values)[i]
    # Perform correlation test with age
    test_result <- cor.test(methyl_values[,i], age, method = method) # Use "spearman" or "kendall" if more appropriate
    
    # Store the results
    correlation_results[[site_name]] <- list(
      correlation_coefficient = test_result$estimate,
      p_value = test_result$p.value
    )
  }
  
  # Optionally, convert the results list to a more convenient format like a dataframe
  correlation_summary <- tibble(
    Site = names(correlation_results),
    Correlation = sapply(correlation_results, function(x) x$correlation_coefficient),
    P_value = sapply(correlation_results, function(x) x$p_value),
    SMR = SMR,
    species = species
  )
  return(correlation_summary)
}

cor.test.age.filter <- function(input, p_value = 0.05) {
  significant_vector <- as.vector(ifelse(input$P_value <= p_value, TRUE, FALSE))
  input$significant <- significant_vector
  return(input)
}


#### Correlation tests ####

AC_cor_age_pearson <- cor.test.age(AC_meth_values_JM, AC_age, AC_methyl_sites_JM$SMR, species = "AC")
AC_cor_age_filtered_pearson <- cor.test.age.filter(AC_cor_age_pearson, 0.05)

# AC_cor_age_spearman <- cor.test.age(AC_meth_values_JM, AC_meth_data$age, method = "spearman")
# AC_cor_age_filtered_spearman <- cor.test.age.filter(AC_cor_age_spearman, 0.05)
# 
# AC_cor_age_kendall <- cor.test.age(AC_meth_values_JM, AC_meth_data$age, method = "kendall")
# AC_cor_age_filtered_kendall <- cor.test.age.filter(AC_cor_age_kendall, 0.05)

AS_cor_age_pearson <- cor.test.age(AS_meth_values_JM, AS_age, AS_methyl_sites_JM$SMR, species = "AS")
AS_cor_age_filtered_pearson <- cor.test.age.filter(AS_cor_age_pearson, 0.05)

EH_cor_age_pearson <- cor.test.age(EH_meth_values_JM, EH_age, EH_methyl_sites_JM$SMR, species = "EH")
EH_cor_age_filtered_pearson <- cor.test.age.filter(EH_cor_age_pearson, 0.05)

JM_cor_age_pearson <- cor.test.age(JM_meth_values_JM, JM_age, JM_methyl_sites_JM$SMR, species = "JM")
JM_cor_age_filtered_pearson <- cor.test.age.filter(JM_cor_age_pearson, 0.05)

ZF_cor_age_pearson <- cor.test.age(ZF_meth_values_JM, ZF_age, ZF_methyl_sites_JM$SMR, species = "ZF")
ZF_cor_age_filtered_pearson <- cor.test.age.filter(ZF_cor_age_pearson, 0.05)

cor_all <- rbind(AC_cor_age_filtered_pearson,
                 AS_cor_age_filtered_pearson,
                 EH_cor_age_filtered_pearson, 
                 JM_cor_age_filtered_pearson,
                 ZF_cor_age_filtered_pearson)

## selecting CpGs

ncol(AC_meth_values_JM[AC_cor_age_filtered_pearson$significant])
ncol(AS_meth_values_JM[AS_cor_age_filtered_pearson$significant])
ncol(EH_meth_values_JM[EH_cor_age_filtered_pearson$significant])
ncol(JM_meth_values_JM[JM_cor_age_filtered_pearson$significant])
ncol(ZF_meth_values_JM[ZF_cor_age_filtered_pearson$significant])


#### plotting ####
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector()
colpal <- hcl.colors(7, "SunsetDark") 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## all
ggplot(cor_all, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species, alpha = significant)) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  # facet_row(~SMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot(cor_all, aes()) +
  geom_point(aes(x = Site, y = Correlation, color = species, alpha = significant)) +
  # geom_line(aes(x = c(-1,1), y = log2(0.05), color = "#CC79A7")) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  facet_wrap(~SMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only significant
ggplot(subset(cor_all, significant == TRUE), aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species)) +
  facet_wrap(~SMR) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
