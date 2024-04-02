library(methylKit)
library(tidyverse)

setwd("/workspace/cfngle/raw-data/")

ZF_file_names <- list.files("ZF/008.methylkit_BM_local/", pattern = "*.txt")
ZF_sample_names <- as.list(substr(ZF_file_names, start = 1, stop = 6))
ZF_file_names <- as.list(paste0("ZF/008.methylkit_BM_local/", ZF_file_names))
ZF_treatment <- c(rep(c(0, 1), length.out = length(ZF_file_names)))

load("/workspace/cfngle/results-data/03_extracted_methyl/methyl_data_ZF_bt2_local.RData")

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

# metadata
samples.age_ZF <- read.csv("ZF/metadata/ZF_metadata.csv", sep=",") 
colnames(samples.age_ZF) <- gsub("Age_Weeks", "age", colnames(samples.age_ZF))

ZF_methyl_data_pos <- meth.extraction(methyl_data_ZF, samples_age = samples.age_ZF, min_per_group = 44L)

ZF_methyl_data <- ZF_methyl_data_pos[[2]] 
ZF_methyl_pos <- ZF_methyl_data_pos[[1]]

save(ZF_methyl_data, file = "ZF/zzz_methyldata/ZF_methyldata_88.RData")
save(ZF_methyl_pos, file = "ZF/zzz_methyldata/ZF_methylpos_88.RData")

write.csv(ZF_methyl_data, file = "ZF/zzz_methyldata/ZF_methyldata_88.csv")
write.csv(ZF_methyl_pos, file = "ZF/zzz_methyldata/ZF_methylpos_88.csv")