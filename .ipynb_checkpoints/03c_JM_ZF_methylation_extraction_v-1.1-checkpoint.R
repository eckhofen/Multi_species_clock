# Preparation of a single methylation matrix from several individual sample files. The process follows the standard procedure suggested by methylKit. After uniting the samples we find that very few CpGs are left and we re-iterate the process using only samples with initial >100.000 CpGs.
# NOTE: for my setup methylKit seems to run properly only in R version below 3.5.3.
# Input data: individual sample files in methylKit format.
# Output data: a methylBase object from the methylKit package containing methylation data from all samples.

# 1. Prepare the environment
# Load required packages
library(methylKit)
library(tidyverse)

# Set working directory: CHANGE TO YOUR WORKING DIRECTORY
setwd("/workspace/cfngle/raw-data/")

# 2. Prepare data
# INPUT: Read individual sample  files
# List files names
# JM_file_names <- list.files("JM/008.methylkit/", pattern = "*.txt") 
# JM_sample_names <- as.list(substr(JM_file_names, start = 1, stop = 6))
# JM_file_names <- as.list(paste0("JM/008.methylkit/", JM_file_names))
# 
# # creating arbitrary vector of zeros and ones for the treatment parameter in methRead
# JM_treatment <- c(rep(c(0, 1), length.out = length(JM_file_names)))

# test <- read.csv(JM_file_names[[1]][[1]], sep = "\t")
# test_methyl <- methRead(JM_file_names[[1]], sample.id = JM_sample_names[1], assembly="ASM223467v1", header=TRUE, mincov = 1, treatment = JM_treatment[1])
# 
# str(JM_file_names)

# 3. Read files. List contains samples names and treatment vector is arbitrary since we don't have 2 groups only.

# does not work for some reason. >>> This script ran in jupyter notebook and in version R 3.5.2  
# methyl_data_JM <- methRead(JM_file_names,
#                             sample.id = JM_sample_names,
#                             assembly="ASM223467v1",
#                             header=TRUE, mincov = 1,
#                             treatment = JM_treatment)

load("/workspace/cfngle/results-data/03_extracted_methyl/methyl_data_JM_JN.RData")

# Check the object
methyl_data_JM

# 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
filtered = filterByCoverage(methyl_data_JM, lo.count=5, lo.perc=NULL, hi.count=100, hi.perc=NULL) 

# Check the object
filtered

# 5. Normalize coverage across samples
norm = normalizeCoverage(filtered)
# Check the object
norm
save(norm, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm.RData")

load("/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm.RData")

# 6. Keep CpGs present in all samples. This option has to be evaluated for each experimental design, here we start with the most conservative approach which is to only keep CpGs present in all samples.
meth = methylKit::unite(norm, destrand=FALSE)

# Check object
meth
dim(meth)
save(meth, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_78621_CpGs.RData")

# creating df for overview 
overview_CpGs <- c(6254, 54772, 66079, 78621)
overview_min_grp <- c(23.5, 18, 17, 16) 

metadata_CpGs <- data.frame(CpGs = overview_CpGs, included_samples = overview_min_grp*2, total_samples = 47, ratio_samples = (overview_min_grp*2)/47 )
write.csv(metadata_CpGs, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_metadata.csv")


# 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
filtered = filterByCoverage(m.data, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) 
# Check the object
filtered

# 5. Normalize coverage across samples
norm = normalizeCoverage(filtered)
# Check the object
norm

# 6. Keep CpGs present in 48 samples per group. We relaxed this criterion to obtain more CpGs but later we will need to missing data.
meth = methylKit::unite(norm, destrand=FALSE, min.per.group=19L)
# Check object


load("/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_66079_CpGs.RData")
meth
dim(meth)
# 85735   334

# 7. Save the object
# save(meth, file="meth-10cov-100000CpGs.Rdata")


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
# Add the variables of interest. In this case we add "age" and "batch". Read the file containing this information in the following format with the samples ordered as in our dataframe (perc.meth.df):
samples.age <- read.csv("JM/raw-reads/00_metadata/00_JM_metadata.csv", sep=",") 
head(samples.age)

meth.age.df <- perc.meth.df %>% mutate(age=samples.age$age, sex = samples.age$sex)

write.csv(meth.age.df, file = "JM/zzz-methyldata/00_JM_methyldata_66079_CpGs.csv")

# Obtain a dataframe listing all the CpGs and their postion on the genome
meth_pos <- as.data.frame(bind_cols(chr = meth$chr, chr_pos = meth$start, strand = meth$strand, CpG_name = cpgs))
write.csv(meth_pos, file = "JM/zzz-methyldata/01_JM_methylpos_66079_CpGs.csv")

#### function ####

meth.extraction <- function(methyl_data, min_per_group = NULL, samples_age, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) {
  temp = filterByCoverage(methyl_data, lo.count=lo.count, lo.perc=lo.perc, hi.count=hi.count, hi.perc=hi.perc) 

  temp = normalizeCoverage(temp)
  
  meth = methylKit::unite(temp, destrand=FALSE, min.per.group=paste0(min_per_group, "L"))
  
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
  meth.percent.to

  meth.age.df <- perc.meth.df %>% mutate(age=samples_age$age, sex = samples_age$sex)
  
  # Obtain a dataframe listing all the CpGs and their position on the genome
  meth_pos <- as.data.frame(bind_cols(chr = meth$chr, chr_pos = meth$start, strand = meth$strand, CpG_name = cpgs))
  return(list(meth_pos,meth.age.df))
}


meth.extraction(methyl_data_JM, samples_age = samples.age)



write.csv(meth.age.df, file = "JM/zzz-methyldata/00_JM_methyldata_66079_CpGs.csv")

