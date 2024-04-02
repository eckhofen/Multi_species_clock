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
JM_file_names <- list.files("JM/008.methylkit/", pattern = "*.txt")
JM_sample_names <- as.list(substr(JM_file_names, start = 1, stop = 6))
JM_file_names <- as.list(paste0("JM/008.methylkit/", JM_file_names))

# creating arbitrary vector of zeros and ones for the treatment parameter in methRead
JM_treatment <- c(rep(c(0, 1), length.out = length(JM_file_names)))

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

# List files names
ZF_file_names <- list.files("ZF/008.methylkit/", pattern = "*.txt")
ZF_sample_names <- as.list(substr(ZF_file_names, start = 1, stop = 6))
ZF_file_names <- as.list(paste0("ZF/008.methylkit/", ZF_file_names))
ZF_treatment <- c(rep(c(0, 1), length.out = length(ZF_file_names)))

methyl_data_ZF <- methRead(ZF_file_names,
                            sample.id = ZF_sample_names,
                            assembly="GRCz11",
                            header=TRUE, mincov = 1,
                            treatment = ZF_treatment)

# creating arbitrary vector of zeros and ones for the treatment parameter in methRead
JM_treatment <- c(rep(c(0, 1), length.out = length(JM_file_names)))

load("/workspace/cfngle/results-data/03_extracted_methyl/methyl_data_JM_bt2_local.RData")

# Check the object
methyl_data_JM

methyl_data_JM_23 <- reorganize(methyl_data_JM, sample.ids = unlist(JM_sample_names)[1:23], treatment = JM_treatment[1:23])
methyl_data_JM_24 <- reorganize(methyl_data_JM, sample.ids = unlist(JM_sample_names)[1:24], treatment = JM_treatment[1:24])

getMethylationStats(methyl_data_JM[[24]],plot=TRUE,both.strands=FALSE)
getMethylationStats(methyl_data_JM_23[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(methyl_data_JM[[1]],plot=TRUE,both.strands=FALSE)


# 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
filtered = filterByCoverage(methyl_data_JM, lo.count=10, lo.perc=NULL, hi.count=100, hi.perc=NULL)
filtered_JM_23 = filterByCoverage(methyl_data_JM_23, lo.count=10, lo.perc=NULL, hi.count=100, hi.perc=NULL)
filtered_JM_24 = filterByCoverage(methyl_data_JM_24, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
filtered_JM_24_100 = filterByCoverage(methyl_data_JM_24, lo.count=10, lo.perc=NULL, hi.count=100, hi.perc=NULL)

getMethylationStats(filtered_JM_24[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(filtered_JM_24[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(filtered_JM_24_100[[1]],plot=TRUE,both.strands=FALSE)


# Check the object
filtered

# 5. Normalize coverage across samples
norm <-  normalizeCoverage(filtered)
norm_JM_23 <- normalizeCoverage(filtered_JM_23)
norm_JM_24 <- normalizeCoverage(filtered_JM_24)

# Check the object
norm
getMethylationStats(norm_JM_24[[10]],plot=TRUE,both.strands=FALSE)
getCoverageStats(norm_JM_24[[1]],plot=TRUE,both.strands=FALSE)

save(norm, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm_JM.RData")
save(norm_JM_23, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm_JM_23.RData")
save(norm_JM_24, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm_JM_24.RData")

load("/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm_JM.RData")
load("/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm_JM_23.RData")
load("/workspace/cfngle/results-data/03_extracted_methyl/TMP_methyl_norm_JM_24.RData")

# 6. Keep CpGs present in all samples. This option has to be evaluated for each experimental design, here we start with the most conservative approach which is to only keep CpGs present in all samples.

meth = methylKit::unite(norm, destrand=FALSE, min.per.group = 17L)
meth_JM_23 = methylKit::unite(norm_JM_23, destrand=FALSE, min.per.group = 10L)
meth_JM_24 = methylKit::unite(norm_JM_24, destrand=FALSE, min.per.group = 10L)

# Check object
meth
dim(meth)
save(meth, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_78621_CpGs.RData")
save(meth_JM_24, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_179818_CpGs.RData")

# creating df for overview 
overview_CpGs <- c(6254, 54772, 66079, 78621)
overview_min_grp <- c(23.5, 18, 17, 16) 

metadata_CpGs <- data.frame(CpGs = overview_CpGs, included_samples = overview_min_grp*2, total_samples = 47, ratio_samples = (overview_min_grp*2)/47 )
write.csv(metadata_CpGs, file ="/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_metadata.csv")

# 
# # 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
# filtered = filterByCoverage(m.data, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) 
# # Check the object
# filtered
# 
# # 5. Normalize coverage across samples
# norm = normalizeCoverage(filtered)
# # Check the object
# norm
# 
# # 6. Keep CpGs present in 48 samples per group. We relaxed this criterion to obtain more CpGs but later we will need to missing data.
# meth = methylKit::unite(norm, destrand=FALSE)
# # Check object


load("/workspace/cfngle/results-data/03_extracted_methyl/TMP_JM_methyl_66079_CpGs.RData")
meth
dim(meth)

# Obtain percent methylation values
perc.meth_JM_24 <-percMethylation(meth_JM_24)

# Transpose dataframe to have samples as rows and CpGs as columns
perc.meth.df_JM_24 <- as.data.frame(t(perc.meth_JM_24))


# Obtain the unique names of CpGs in the form of chromosome.start
meth.df_JM_24 <- as.data.frame(meth_JM_24)
cpg.df <- bind_cols(chr=meth_JM_24$chr, start=meth_JM_24$start)
cpgs <- tidyr::unite(cpg.df, cpgs, chr:start, sep=":")
# Add the unique CpG names as column names in the dataframe
colnames(perc.meth.df_JM_24) <- t(cpgs)
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
  print("List will be returned with first element being a dataframe containing the mthylation positions and the second one a dataframe containing the metadata and names")
  return(list(meth_pos,meth.age.df))
}

samples.age_JM <- read.csv("JM/raw-reads/00_metadata/00_JM_metadata.csv", sep=",") 
samples_age_JM_24 <- samples.age_JM[1:24,]

# test <- meth.extraction(methyl_data_JM_24, samples_age = samples_age_JM_24, min_per_group = 10L)

JM_24_methyl_data_pos <- meth.extraction(methyl_data_JM_24, samples_age = samples_age_JM_24, min_per_group = 10L)

JM_24_methyl_data <- JM_24_methyl_data_pos[[2]]
JM_24_methyl_pos <- JM_24_methyl_data_pos[[1]]


write.csv(JM_24_methyl_data_pos[[2]], file = "JM/zzz-methyldata/00_JM_methyldata_179818_CpGs.csv")

write.csv(JM_24_methyl_data_pos[[1]], file = "JM/zzz-methyldata/01_JM_methylpos_179818_CpGs.csv")

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
