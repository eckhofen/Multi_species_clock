library(dplyr)
# library()

setwd("/workspace/cfngle/raw-data/JM")
JM_name_list <- read.csv("raw-reads/00_metadata/name_list_PRJNA716946.txt", sep = "\t")
JM_downloaded_order <- c("SAMN21040240", "SAMN21040241", "SAMN21040242", "SAMN21040243", "SAMN21040244", "SAMN21040245", "SAMN21040246", 
                         "SAMN21040247", "SAMN21040248", "SAMN21040249", "SAMN21040250", "SAMN21040251", "SAMN21040252", "SAMN21040253", "SAMN21040254", 
                         "SAMN21040255", "SAMN21040256", "SAMN21040257", "SAMN21040258", "SAMN21040259", "SAMN21040260", "SAMN21040261", "SAMN21040262", "SAMN21040263", 
                         "SRR18462633", "SRR18462634", "SRR18462635", "SRR18462636", "SRR18462637", "SRR18462638", "SRR18462639", "SRR18462640", 
                         "SRR18462641", "SRR18462642", "SRR18462643", "SRR18462644", "SRR18462645", "SRR18462646", "SRR18462647", "SRR18462648", 
                         "SRR18462649", "SRR18462650", "SRR18462651", "SRR18462652", "SRR18462653", "SRR18462654", "SRR18462655")
JM_names <- sprintf("JM_%03d", 1:47)

JM_downloaded_order %in% JM_name_list$sample_accession

JM_downloaded_order_new <- c()
for(i in 1:length(JM_name_list$sample_accession)) {
  if(JM_downloaded_order[i] %in% JM_name_list$sample_accession){
    JM_downloaded_order_new[i] <- JM_name_list$sample_accession[JM_name_list$sample_accession == JM_downloaded_order[i]] 
    print(JM_downloaded_order_new[i])
    }
  else if(JM_downloaded_order[i] %in% JM_name_list$run_accession){
    JM_downloaded_order_new[i] <- JM_name_list$sample_accession[JM_name_list$run_accession == JM_downloaded_order[i]] 
    print(JM_downloaded_order_new[i])
    }
}

JM_name_list <- JM_name_list[match(JM_downloaded_order_new, JM_name_list$sample_accession),]
cbind(JM_names, JM_name_list)

JM_metadata_names <- list.files(path = "raw-reads/00_metadata/", pattern = "\\.txt$", full.names = TRUE)
JM_metadata_names <- JM_metadata_names[-1]

JM_metadata <- list()

# Loop through files, read each file, and combine
for(file in JM_metadata_names) {
  # Read the current file
  temp_data <- read.table(file, header = TRUE, sep = "\t") # Adjust parameters as needed
  
  # Combine with the main dataframe
  JM_metadata[[file]] <- temp_data
}

## extract identification variables
JM_metadata_indent <- data.frame(identifiers = c(sapply(JM_metadata, function(df) df[1,1])), index = c(1:47)) 
JM_metadata_indent[,1] <- gsub("/.$.txt", "", JM_metadata_indent[,1])

## Age
JM_metadata_age <- data.frame(identifiers = c(sapply(JM_metadata, function(df) df[5,1])), index = c(1:47)) 
JM_metadata_age[,1] <- gsub("/.$.txt", "", JM_metadata_age[,1])

## sex
JM_metadata_sex <- data.frame(identifiers = c(sapply(JM_metadata, function(df) df[6,1])), index = c(1:47)) 
JM_metadata_sex[,1] <- gsub("/.$.txt", "", JM_metadata_sex[,1])

temp_cat <- strsplit(JM_metadata_indent[,1], ";")
temp_cat <- sapply(temp_cat, function(x) x[1])

JM_metadata_df <- cbind(JM_name_list, as.data.frame(temp_cat))

names(JM_metadata_df)[11] <- "biosample"
JM_metadata_df[,11] <- gsub("Identifiers: BioSample: ", "", JM_metadata_df[,11])


JM_metadata_df <- cbind(JM_metadata_df, as.data.frame(JM_metadata_age$identifiers), as.data.frame(JM_metadata_sex$identifiers))
JM_metadata_df[,12] <- gsub("    /age=", "", JM_metadata_df[,12]) %>% 
  gsub(" days post hatch", "", .)

JM_metadata_df[,13] <- gsub("    /sex=", "", JM_metadata_df[,13])

names(JM_metadata_df)[12:13] <- c("age", "sex")

JM_metadata_df_final <- JM_metadata_df %>% 
  dplyr::select(-matches("ftp"))

JM_metadata_df_final$age <- as.numeric(JM_metadata_df_final$age)
JM_metadata_all <- cbind(JM_names, JM_metadata_df_final)

hist(JM_metadata_df_final$age)

write.csv(JM_metadata_all, file = "raw-reads/00_metadata/00_JM_metadata.csv")
