#### Overview ####
# extraction of metyltion values for all species in all samples

#### Preparation ####
library(methylKit)
library(tidyverse)

setwd("/workspace/cfngle/raw-data/")

#### Data preparation ####
##JM
JM_file_names <- list.files("JM/008.methylkit/", pattern = "*.txt")
JM_sample_names <- as.list(substr(JM_file_names, start = 1, stop = 6))
JM_file_names <- as.list(paste0("JM/008.methylkit/", JM_file_names))
# creating arbitrary vector of zeros and ones for the treatment parameter in methRead
JM_treatment <- c(rep(c(0, 1), length.out = length(JM_file_names)))

##ZF
ZF_file_names <- list.files("ZF/008.methylkit_BM_local/", pattern = "*.txt")
ZF_sample_names <- as.list(substr(ZF_file_names, start = 1, stop = 6))
ZF_file_names <- as.list(paste0("ZF/008.methylkit_BM_local/", ZF_file_names))
ZF_treatment <- c(rep(c(0, 1), length.out = length(ZF_file_names)))

### preparing methylKitList objects (does only work in Jupyter notebook and R version < 3.5.3)

# does not work for some reason. >>> This script ran in jupyter notebook and in version R 3.5.2  
# methyl_data_JM <- methRead(JM_file_names,
#                             sample.id = JM_sample_names,
#                             assembly="ASM223467v1",
#                             header=TRUE, mincov = 1,
#                             treatment = JM_treatment)

# 
# methyl_data_ZF <- methRead(ZF_file_names,
#                             sample.id = ZF_sample_names,
#                             assembly="GRCz11",
#                             header=TRUE, mincov = 1,
#                             treatment = ZF_treatment)

## loading methylKitList objects which were saved in jupyter notebook
## JM
load("/workspace/cfngle/results-data/03_extracted_methyl/methyl_data_JM_bt2_local.RData")

## ZF
load("/workspace/cfngle/results-data/03_extracted_methyl/methyl_data_ZF_bt2_local.RData")

# Check the objects
methyl_data_JM
methyl_data_ZF


total_CpGs_ZF <- sapply(methyl_data_ZF, function(x) nrow(x))

boxplot(total_CpGs_ZF)

## JM
# for JM, the data did not contain enough CpGs and had lower read coverage in the samples from number 24 (second batch - see paper). Thus only the first batch was selected which still contains three different age groups
methyl_data_JM_24 <- reorganize(methyl_data_JM, sample.ids = unlist(JM_sample_names)[1:24], treatment = JM_treatment[1:24])


#### Data manipulation ####
# all the steps below are included in one single function for ease of use (jump to section "function"). Here, the included steps + some graphical representation of the data are explained

## prior distribution
# before manipulating the data, the distribution can be looked at to identify PCR bias

getMethylationStats(methyl_data_JM_24[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(methyl_data_JM[[1]],plot=TRUE,both.strands=FALSE)

## filtering 
# samples with certain amounts of reads are filtered out 
# reads with very low number of reads (<10) and exceeding the 99.9% percentile are filtered out in this case

filtered_JM_24 = filterByCoverage(methyl_data_JM_24, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

# in the JM paper, they used lo.count = 5 and hi.count = 100. It seems that the option above gives better results (graphic comparison, and more CpGs)
filtered_JM_24_100 = filterByCoverage(methyl_data_JM_24, lo.count=10, lo.perc=NULL, hi.count=100, hi.perc=NULL)

getMethylationStats(filtered_JM_24[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(filtered_JM_24[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(filtered_JM_24_100[[1]],plot=TRUE,both.strands=FALSE)

## normalizing
# normalize coverage across samples again to eliminate PCR bias
norm_JM_24 <- normalizeCoverage(filtered_JM_24)

getMethylationStats(norm_JM_24[[10]],plot=TRUE,both.strands=FALSE)
getCoverageStats(norm_JM_24[[1]],plot=TRUE,both.strands=FALSE)

## getting shared methylation sites
# this allows to extract CpGs which are present in all samples. The parameter min.per.group, lets us define how many samples per treatment group have to have the CpG site in order to keep it. The lower the number, the higher the CpGs. 

# in this case CpGs have to be present in at least 20 out of 24 samples (83%)
meth_JM_24 = methylKit::unite(norm_JM_24, destrand=FALSE, min.per.group = 10L)

# check how many CpGs were kept
dim(meth_JM_24)

# Check object
save(meth_JM_24, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_179818_CpGs.RData")

### getting methylation amount in percent
# Obtain percent methylation values

##JM
perc.meth_JM_24 <-percMethylation(meth_JM_24)

# Transpose dataframe to have samples as rows and CpGs as columns
perc.meth.df_JM_24 <- as.data.frame(t(perc.meth_JM_24))

# Obtain the unique names of CpGs in the form of chromosome.start
meth_df_JM_24 <- as.data.frame(meth_JM_24)
cpg_df_JM <- bind_cols(chr=meth_JM_24$chr, start=meth_JM_24$start)
cpgs_JM <- tidyr::unite(cpg_df_JM, cpgs_JM, chr:start, sep=":")

# Add the unique CpG names as column names in the dataframe
colnames(perc.meth.df_JM_24) <- t(cpgs_JM)
# Add the variables of interest. In this case we add "age" and "batch". Read the file containing this information in the following format with the samples ordered as in our dataframe (perc.meth.df):

samples.age_JM <- read.csv("JM/raw-reads/00_metadata/00_JM_metadata.csv", sep=",") 
head(samples.age)

meth.age.df <- perc.meth.df %>% mutate(age=samples.age$age, sex = samples.age$sex)

write.csv(meth.age.df, file = "JM/zzz-methyldata/00_JM_methyldata_66079_CpGs.csv")

# Obtain a dataframe listing all the CpGs and their position on the genome
meth_pos <- as.data.frame(bind_cols(chr = meth$chr, chr_pos = meth$start, strand = meth$strand, CpG_name = cpgs))
write.csv(meth_pos, file = "JM/zzz-methyldata/01_JM_methylpos_179818_CpGs.csv")

#### function ####

meth.extraction <- function(methyl_data, min_per_group = NULL, samples_age, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) {
  print("filtering...")
  temp = filterByCoverage(methyl_data, lo.count=lo.count, lo.perc=lo.perc, hi.count=hi.count, hi.perc=hi.perc) 
  print("filtered, normalizing...")
  temp = normalizeCoverage(temp)
  print("normalized, uniting all samples...")
  meth = methylKit::unite(temp, destrand=FALSE, min.per.group=min_per_group)
  print("united, further downstream processes...")
  # Obtain percent methylation values
  perc.meth=percMethylation(meth)
  # Transpose dataframe to have samples as rows and CpGs as columns
  perc.meth.df <- as.data.frame(t(perc.meth))
  
  # Obtain the unique names of CpGs in the form of chromosome.start
  meth.df <- as.data.frame(meth)
  cpg.df <- bind_cols(chr=meth$chr, start=meth$start)
  cpgs <- tidyr::unite(cpg.df, cpgs, chr:start, sep=":")
  # Add the unique CpG names as column names in the dataframe
  colnames(perc.meth.df) <- t(cpgs)
  
  meth.age.df <- perc.meth.df %>% mutate(age=samples_age$age, sex = samples_age$sex)
  
  # Obtain a dataframe listing all the CpGs and their position on the genome
  meth_pos <- as.data.frame(bind_cols(chr = meth$chr, chr_pos = meth$start, strand = meth$strand, CpG_name = cpgs))
  print("Done!")
  print("List will be returned with first element being a dataframe containing the methylation positions and the second one a dataframe containing the metadata and names")
  return(list(meth_pos,meth.age.df))
}

## JM
# getting metadata
samples.age_JM <- read.csv("JM/raw-reads/00_metadata/00_JM_metadata.csv", sep=",") 
samples_age_JM_24 <- samples.age_JM[1:24,]

# extracting methyl values and positions
JM_24_methyl_data_pos <- meth.extraction(methyl_data_JM_24, samples_age = samples_age_JM_24, min_per_group = 8L)

JM_24_methyl_data <- JM_24_methyl_data_pos[[2]]
JM_24_methyl_pos <- JM_24_methyl_data_pos[[1]]

save(JM_24_methyl_data, file = "JM/zzz-methyldata/00_JM_methyldata_243285_CpGs.Rdata")
save(JM_24_methyl_pos, file = "JM/zzz-methyldata/00_JM_methylpos_243285_CpGs.Rdata")

write.csv(JM_24_methyl_data_pos[[2]], file = "JM/zzz-methyldata/00_JM_methyldata_179818_CpGs.csv")

write.csv(JM_24_methyl_data_pos[[1]], file = "JM/zzz-methyldata/01_JM_methylpos_179818_CpGs.csv")

## ZF
# File size too large, was done with HPC 
# metadata
samples.age_ZF <- read.csv("ZF/metadata/ZF_metadata.csv", sep=",") 
colnames(samples.age_ZF) <- gsub("Age_Weeks", "age", colnames(samples.age_ZF))

ZF_methyl_data_pos <- meth.extraction(methyl_data_ZF, samples_age = samples.age_ZF, min_per_group = 48L)

ZF_methyl_data <- ZF_methyl_data_pos[[2]] 
ZF_methyl_pos <- ZF_methyl_data_pos[[1]]

save(ZF_methyl_data, file = "ZF/zzz_methyldata/ZF_methyldata.RData")
save(ZF_methyl_pos, file = "ZF/zzz_methyldata/ZF_methylpos.RData")

load("ZF/zzz_methyldata/ZF_methyldata.RData")
load("ZF/zzz_methyldata/ZF_methylpos.RData")

#### corr testing and plotting ####
meth_ <- meth
meth_@sample.ids <- as.character(JM_metadata$age)

getCorrelation(methylKit::select(meth,c(1:5)))

clusterSamples(meth_, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth_JM_23, screeplot=TRUE)
PCASamples(meth_JM_23)

JM_metadata <- read.csv("JM/raw-reads/00_metadata/00_JM_metadata.csv")
JM_age <- data.frame(age = JM_metadata$age)


PC_age_JM <- assocComp(mBase=meth,JM_age)
sort(PC_age_JM$association)
##38 , 33

df_temp_PCA <- data.frame(ID = JM_metadata$JM_names, PC1 = as.vector(PC_age_JM$pcs[,1]),PC10 = as.vector(PC_age_JM$pcs[,10]), PC2 = as.vector(PC_age_JM$pcs[,2]),PC3 = as.vector(PC_age_JM$pcs[,3]), PC33 = as.vector(PC_age_JM$pcs[,33]), PC38 = as.vector(PC_age_JM$pcs[,38]), age = JM_metadata$age)

meth@treatment <- rep(c(1,2,3), 10)

tiles = tileMethylCounts(meth,win.size=1000,step.size=1000,cov.bases = 10)

myDiff=calculateDiffMeth(meth)


library(ggplot2)


temp_df <- data.frame(chr = as.factor(meth$chr), pos = meth$start, coverage = meth$coverage1)

ggplot(temp_df, aes(x = pos, y = coverage, alpha = 0.5)) +
  geom_point() +  
  facet_wrap(~ chr, scales = "free_x") +
  labs(x = "Position", y = "", title = "JM medaka CpGs associated with age (in weeks)") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))

ggplot(df_temp_PCA, aes(x = PC10, y = PC3, color = age, cex = 2)) +
  geom_point() +  
  labs(title = "JM medaka CpGs associated with age (in weeks)") +
  scale_color_gradient(low = "red", high = "darkblue") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))
