#### Overview ####
# retrieving age data on all samples for all species 

#### Preparation ####
library(ggplot2)
library(ggforce)
library(dplyr)

#### Retrieving data ####

##AC
xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/meth-corrected-batchcorrected-cod.Rdata")
# xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/Meth-complete-nobatchcorrection-cod.RData")
assign("AC_meth_data", get(xx))
AC_meth_data <- as.data.frame(AC_meth_data)
AC_meta_data <- read.csv("/workspace/cfngle/raw-data/AC/zzz_methyl_data/cod-samples.txt", sep = "\t")

# tail(colnames(AC_meth_data))
AC_age <- AC_meth_data$age
AC_sample_names <- rownames(AC_meth_data) 
AC_names <- sprintf("AC_%03d", 1:length(AC_sample_names))
AC_max_age <- 25 

AC_sample_sex <- AC_meta_data$sex %>% 
  gsub("M", "Male", .) %>%
  gsub("F", "Female", .) 

AC_sample_sex[is.na(AC_sample_sex)] <- "unknown"

AC_df <- data.frame(sample_name = AC_sample_names,
                    name = AC_names,
                    age = AC_age,
                    sex = AC_sample_sex,
                    max_age = AC_max_age,
                    rel_age = AC_age / AC_max_age, 
                    species = "AC",
                    scientific_name = "Gadus morhua")

##AS
xx <- load("/workspace/cfngle/raw-data/AS/zzz_methyl_data/Meth-complete-snapper.RData")
assign("AS_meth_data", get(xx))
# tail(colnames(AS_meth_data))
AS_age <- AS_meth_data$age
AS_sample_names <- rownames(AS_meth_data)
AS_names <- sprintf("AS_%03d", 1:length(AS_sample_names))
AS_max_age <- 54 
AS_sample_sex <- "unknown"

AS_df <- data.frame(sample_name = AS_sample_names,
                    name = AS_names,
                    age = AS_age,
                    sex = AS_sample_sex,
                    max_age = AS_max_age,
                    rel_age = AS_age / AS_max_age, 
                    species = "AS",
                    scientific_name = "Chrysophrys auratus")

##EH
xx <- load("/workspace/cfngle/raw-data/EH/zzz-methyl_data/Meth-complete-hake.RData")
assign("EH_meth_data", get(xx))
tail(colnames(EH_meth_data))

EH_meta_data <- read.csv("/workspace/cfngle/raw-data/EH/zzz-methyl_data/hake-samples.txt", sep = "\t")

EH_age <- EH_meta_data$age
EH_sample_names <- rownames(EH_meth_data)
EH_names <- sprintf("EH_%03d", 1:length(EH_sample_names))
EH_max_age <- 20 
EH_sample_sex <- EH_meta_data$sex %>% 
  gsub("M", "Male", .) %>%
  gsub("F", "Female", .)
  
EH_df <- data.frame(sample_name = EH_sample_names,
                    name = EH_names,
                    age = EH_age,
                    sex = EH_sample_sex,
                    max_age = EH_max_age,
                    rel_age = EH_age / EH_max_age, 
                    species = "EH",
                    scientific_name = "Merluccius merluccius")

##JM
# JM_meth_data <- read.csv("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_66079_CpGs.csv")
JM_meth_data <- load("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_179818_CpGs.Rdata")
JM_meth_data <- JM_24_methyl_data

JM_age <- JM_meth_data$age / 365
JM_sample_names <- rownames(JM_meth_data)
JM_names <- sprintf("JM_%03d", 1:length(JM_sample_names))
JM_max_age <- 5 
JM_sample_sex <- "Male"
JM_df <- data.frame(sample_name = JM_sample_names,
                    name = JM_sample_names,
                    age = JM_age,
                    sex = JM_sample_sex,
                    max_age = JM_max_age,
                    rel_age = JM_age / JM_max_age, 
                    species = "JM",
                    scientific_name = "Oryzias latipes")
  
##ZF
ZF_meta_data <- read.csv("/workspace/cfngle/raw-data/ZF/metadata/ZF_metadata.csv")
ZF_age <- ZF_meta_data$Age_Weeks / 52
ZF_sample_names <- ZF_meta_data$name_new
ZF_max_age <- 5.5
ZF_sample_sex <- ZF_meta_data$Sex

ZF_df <- data.frame(sample_name = ZF_sample_names,
                    name = ZF_sample_names,
                    age = ZF_age,
                    sex = ZF_sample_sex,
                    max_age = ZF_max_age,
                    rel_age = ZF_age / ZF_max_age, 
                    species = "ZF",
                    scientific_name = "Danio rerio")

#### Data manipulation ####

df_all <- rbind(AC_df, AS_df, EH_df, JM_df, ZF_df)

# save files
save(AC_df, AS_df, EH_df, JM_df, ZF_df, file = "000_data/000_metadata/metadata_age.RData")

write.csv(AC_df, file = "000_data/000_metadata/AC_age.csv")
write.csv(AS_df, file = "000_data/000_metadata/AS_age.csv")
write.csv(EH_df, file = "000_data/000_metadata/EH_age.csv")
write.csv(JM_df, file = "000_data/000_metadata/JM_age.csv")
write.csv(ZF_df, file = "000_data/000_metadata/ZF_age.csv")

#### Plotting ####

ggplot(df_all, aes(y = rel_age, color = species)) +
  # geom_point(aes(cex = 4, y = rel_age, x = 0.5)) +
  geom_sina(aes(x = species, cex = 1, shape = sex, alpha = 0.5)) +
  geom_boxplot(aes(x = species), fill = NA) +
  theme_classic() +
  labs(title = "Relative age distribution of Cod (AC), Snapper (AS), Hake (EH), Medaka (JM)", x = "Species", y = "Relative age")

ggplot(df_all) +
  # geom_point(aes(cex = 4, y = rel_age, x = 0.5)) +
  geom_sina(aes(x = species, y = rel_age, shape = sex, alpha = 0.5, color = species), cex = 1) +
  # geom_histogram(aes(x = rel_age), position = "dodge") +
  # geom_density(aes(x = rel_age, alpha = 0.3)) +
  theme_classic() +
  labs(title = "Relative age distribution of Cod (AC), Snapper (AS), Hake (EH), Medaka (JM)", x = "Species", y = "Relative age")
