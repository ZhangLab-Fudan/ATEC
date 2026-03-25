rm(list = ls()) # remove environment variables
library(tidyverse)

# set working directory
setwd("/data/whu/tcga_project/script/1_quantification_SQuIRE")

# set cancer type
cancer = "BRCA"

# import data
info <- data.table::fread("/data/whu/home/te_database/results/0_tcga_information.txt", header = T, data.table = F, check.names = F)
Dir <- paste0("/home/whu/whu/tcga_project/1_bam/", cancer)
bam <- list.files(path = Dir)

rt <- data.frame()
for(b in bam){
  # The path should be validated in container
  out_dir <- file.path("/data/whu/tcga_project/3_SQuIRE", cancer, b)
  command_1 <- paste0("mkdir -p ", out_dir)
  command_2 <- paste0("squire Count -m ", file.path("/data/whu/tcga_project/1_bam", cancer, b), " -c ", "/data/whu/tcga_project/3_SQuIRE/squire_clean", 
                      " -o ", out_dir,
                      " -f ", "/data/whu/tcga_project/3_SQuIRE/squire_fetch", 
                      " -r ", info$average_read_length[info$file_id == b],
                      " -b hg38 -p 3 --verbosity")
  rt <- rbind(rt, command_1, command_2)
}

# export data
write.table(rt, paste0(cancer, "_SQuIRE.sh"), sep = "\t", row.names = F, col.names = F, quote = F)
