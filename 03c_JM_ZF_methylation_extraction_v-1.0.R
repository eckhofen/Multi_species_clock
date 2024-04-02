# Preparation of a single methylation matrix from several individual sample files. The process follows the standard procedure suggested by methylKit. After uniting the samples we find that very few CpGs are left and we re-iterate the process using only samples with initial >100.000 CpGs.
# NOTE: for my setup methylKit seems to run properly only in R version below 3.5.3.
# Input data: individual sample files in methylKit format.
# Output data: a methylBase object from the methylKit package containing methylation data from all samples.

# 1. Prepare the environment
# Load required packages
library(methylKit)
# library(tidyverse)

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

# does ot work for some reason
methyl_data_JM <- methRead(JM_file_names,
                            sample.id = JM_sample_names,
                            assembly="ASM223467v1",
                            header=TRUE, mincov = 1,
                            treatment = JM_treatment)

# workaround 
JM_methylation_list <- lapply(seq_along(JM_file_names), function(i) {
  methRead(JM_file_names[[i]],
           sample.id = JM_sample_names[i],
           treatment = JM_treatment[i],
           assembly = "ASM223467v1",
           mincov = 1,
           header = TRUE)
})

methyl_data_JM <- methRead(JM_file_names[[1]][[1]], sample.id = JM_sample_names[1], assembly="ASM223467v1", treatment = JM_treatment[1])

# Check the object
methyl_data_JM

# 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
filtered = filterByCoverage(m.data, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) 
# Check the object
filtered

# 5. Normalize coverage across samples
norm = normalizeCoverage(filtered)
# Check the object
norm

# 6. Keep CpGs present in all samples. This option has to be evaluated for each experimental design, here we start with the most conservative approach which is to only keep CpGs present in all samples.
meth = methylKit::unite(norm, destrand=FALSE)
# Check object
meth
# dim(meth)
# 405 364
# RESULT: Only 405 CpGs are left after the procedure. This is likely to occur if few samples have very low number of reads and are driving down all samples CpGs after the unite function. Our strategy in this case will be to identify these samples, eliminate them and start the process again.

# START AGAIN. 2. Prepare data
# INPUT: Read individual sample  files
# List files names
file.list <- list('102F.txt',	'106F.txt',	'107F.txt',	'111F.txt',	'118F.txt',	'123F.txt',	'124F.txt',	'12F.txt',	'136F.txt',	'140F.txt',	'144F.txt',	'148F.txt',	'154F.txt',	'155F.txt',	'158F.txt',	'15F.txt',	'169F.txt',	'175F.txt',	'181F.txt',	'183F.txt',	'190F.txt',	'195F.txt',	'196F.txt',	'214F.txt',	'223F.txt',	'228F.txt',	'22F.txt',	'238F.txt',	'242F.txt',	'247F.txt',	'248F.txt',	'249F.txt',	'250F.txt',	'37F.txt',	'43F.txt',	'48F.txt',	'53F.txt',	'57F.txt',	'58F.txt',	'60F.txt',	'65F.txt',	'6F.txt',	'71F.txt',	'72F.txt',	'74F.txt',	'75F.txt',	'76F.txt',	'7F.txt',	'81F.txt',	'82F.txt',	'92F.txt',	'96F.txt',	'97F.txt',	'9F.txt',	'275.txt',	'1.txt',	'2.txt',	'4.txt',	'5.txt',	'8.txt',	'38.txt',	'39.txt',	'40.txt',	'41.txt',	'47.txt',	'54.txt',	'56.txt',	'61.txt',	'62.txt',	'70.txt',	'79.txt',	'88.txt',	'90.txt',	'94.txt',	'98.txt',	'117.txt',	'141.txt',	'26.txt',	'52.txt',	'55.txt',	'63.txt',	'80.txt',	'112.txt',	'128.txt',	'151.txt',	'153.txt',	'167.txt',	'191.txt',	'193.txt',	'212.txt',	'222.txt',	'225.txt',	'21.txt',	'23.txt',	'28.txt',	'29.txt',	'83.txt',	'84.txt',	'127.txt',	'161.txt',	'163.txt',	'165.txt',	'166.txt',	'188.txt',	'199.txt',	'255.txt',	'257.txt',	'258.txt',	'261.txt',	'264.txt') 

# 3. Read files. List contains samples names and treatment vector is arbitrary since we don't have 2 groups only.
m.data <- methRead(file.list, sample.id=list('y5s1', 'y1s1',	'y7s1',	'y3s1',	'y2s1',	'y1s2',	'y1s3',	'y3s2',	'y1s4',	'y1s5',	'y2s2',	'y3s3',	'y1s6',	'y1s7',	'y2s3',	'y2s4',	'y2s5',	'y1s8',	'y1s9',	'y1s10',	'y0s1',	'y0s2',	'y1s11',	'y4s1',	'y1s12',	'y1s13',	'y2s8',	'y4s2',	'y0s3',	'y4s3',	'y6s1',	'y5s3',	'y6s2',	'y1s14',	'y4s4',	'y7s2',	'y2s9',	'y2s10',	'y1s15',	'y3s6',	'y4s5',	'y0s4',	'y4s6',	'y3s7',	'y0s6',	'y4s7',	'y3s8',	'y5s4',	'y4s8',	'y5s5',	'y0s7',	'y4s9',	'y3s9',	'y3s10',	'y0s8',	'y1s16',	'y1s17',	'y1s18',	'y1s19',	'y1s20',	'y1s21',	'y1s22',	'y1s23',	'y1s24',	'y1s25',	'y1s26',	'y1s27',	'y1s28',	'y1s29',	'y1s30',	'y1s32',	'y1s33',	'y1s34',	'y1s35',	'y1s36',	'y1s37',	'y1s38',	'y2s11',	'y2s12',	'y2s13',	'y2s14',	'y2s15',	'y2s16',	'y2s17',	'y2s18',	'y2s19',	'y2s20',	'y2s21',	'y2s22',	'y2s23',	'y2s24',	'y2s25',	'y3s11',	'y3s12',	'y3s13',	'y3s14',	'y3s15',	'y3s16',	'y3s17',	'y3s19',	'y3s20',	'y3s21',	'y3s22',	'y3s24',	'y3s25',	'y4s10',	'y4s11',	'y4s12',	'y4s13',	'y4s14'),  
                                             assembly="cod", mincov = 1, treatment=c(0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1, 1))
# Check the object
m.data

# 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
filtered = filterByCoverage(m.data, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) 
# Check the object
filtered

# 5. Normalize coverage across samples
norm = normalizeCoverage(filtered)
# Check the object
norm

# 6. Keep CpGs present in 48 samples per group. We relaxed this criterion to obtain more CpGs but later we will need to missing data.
meth = methylKit::unite(norm, destrand=FALSE, min.per.group=48L)
# Check object
meth
# dim(meth)
# 85735   334

# 7. Save the object
save(meth, file="meth-10cov-100000CpGs.Rdata")


# Obtain percent methylation values
perc.meth=percMethylation(meth)
# Transpose dataframe to have samples as rows and CpGs as columns
perc.meth.df <- as.data.frame(t(perc.meth))
# Obtain the unique names of CpGs in the form of chromosome.start
meth.df <- as.data.frame(meth)
cpg.df <- bind_cols(chr=meth$chr, start=meth$start)
cpgs <- tidyr::unite(cpg.df, cpgs, chr:start, sep=".")
# Add the unique CpG names as column names in the dataframe
colnames(perc.meth.df) <- t(cpgs)
# Add the variables of interest. In this case we add "age" and "batch". Read the file containing this information in the following format with the samples ordered as in our dataframe (perc.meth.df):
samples.age <- read.table("samples.txt", sep="\t", stringsAsFactors = FALSE, quote="", header=TRUE) 
head(samples.age)
#filename sample_file age sample batch
#1 102F.txt        102F   5   y5s1     1
#2 106F.txt        106F   1   y1s1     1
#3 107F.txt        107F   7   y7s1     1
#4 111F.txt        111F   3   y3s1     1
#5 118F.txt        118F   2   y2s1     1
#6 123F.txt        123F   1   y1s2     1
# Rownames of our dataframe (perc.meth.df) are in the same order
meth.age.df <- perc.meth.df %>% mutate(age=samples.age$age)