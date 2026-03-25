rm(list = ls()) # remove environment variables
library(tidyverse)

# set working directory
setwd("/home/whu/whu/tcga_project/script/2_quantification_TElocal")

# set cancer type
cancer = "OV"

# import data
Dir <- paste0("/home/whu/whu/tcga_project/1_bam/", cancer)
bam <- list.files(path = Dir, pattern = "rna_seq.genomic.gdc_realn.bam$", full.names = T, recursive = T)

rt <- data.frame()
for(b in bam){
  id <- gsub("......................................rna_seq.genomic.gdc_realn.bam", "", b)
  id <- gsub(".*\\/", "", id)
  
  out_dir <- file.path("/data/whu/tcga_project/2_TElocal", cancer, id)
  
  command_1 <- paste0("mkdir ", out_dir, ";  cd ", out_dir)
  command_2 <- paste0("TElocal --sortByPos --mode multi -b ", b,  
                      " --GTF /data/whu/reference/gtf/human/gencode.v35.annotation.gtf --TE /data/whu/home/te_database/data/TElocal_gtf/hg38_rmsk_TE.gtf.locInd --project ", id)
  rt <- rbind(rt, command_1, command_2)
}

# export data
write.table(rt, paste0(cancer, "_TElocal.sh"), sep = "\t", row.names = F, col.names = F, quote = F)
