# short script to convert cd-hit terminal logs into readable files which show the amount of clusters which were formed

library(dplyr)
library(stringr)

label_readout <- c("AC_AS", "AC_EH", "AS_EH", "AS_AC", "EH_AC", "EH_AS")
label_bp <- c(rep("100", 6),rep("200", 6),rep("500", 6),rep("1000", 6))

readout <- read.table("/workspace/cfngle/results-data/cd-hit/job_output_cd-hit4325390.txt", sep = "\t") 
readout_filtered <- readout %>% 
  filter(.,str_detect(.[,1], fixed("compared"))) %>% 
  filter(.,str_detect(.[,1], fixed("clusters"))) %>% 
  transmute(str_replace_all(.[,1],"..........        0  compared ", "")) %>% 
  transmute(str_replace_all(.[,1]," clusters", ""))
readout_filtered
readout_filtered <- cbind(readout_filtered, label_readout, label_bp)
names(readout_filtered) <- c("alignments", "species", "seq_length")

write.csv(readout_filtered, "/workspace/cfngle/results-data/cd-hit/filtered_reads.txt")
