# ==============================================================================
# Partial Correlation Analysis: TE Element vs Drug Response
# Control for: Age, Gender, Plate and Tumor Purity
# ==============================================================================

# Clean environment
rm(list = ls())

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(ppcor)
})

args <- commandArgs(T)

# ==============================================================================
# Configuration
# ==============================================================================

software <- "SQuIRE"
max_cores <- 50

# Set working directory
setwd("/data/whu/home/ATEC/resutls/4_drug_TE_cor_SQuIRE")

# ==============================================================================
# Load Reference Data
# ==============================================================================

# Drug response data
load("/home/whu/whu/tcga_project/TCGA_drug/pan_CTRP1_Log2.RData")
ctrp1 <- as.data.frame(pan_CTRP1_Log2)
load("/home/whu/whu/tcga_project/TCGA_drug/pan_CTRP2_Log2.RData")
ctrp2 <- as.data.frame(pan_CTRP2_Log2)
load("/home/whu/whu/tcga_project/TCGA_drug/pan_GDSC1_Log2.RData")
gdsc1 <- as.data.frame(pan_GDSC1_Log2)
load("/home/whu/whu/tcga_project/TCGA_drug/pan_GDSC2_Log2.RData")
gdsc2 <- as.data.frame(pan_GDSC2_Log2)

# Add database labels
rownames(ctrp1) <- paste0("CTRP1_", rownames(ctrp1))
rownames(ctrp2) <- paste0("CTRP2_", rownames(ctrp2))
rownames(gdsc1) <- paste0("GDSC1_", rownames(gdsc1))
rownames(gdsc2) <- paste0("GDSC2_", rownames(gdsc2))

# Merge drug data
drug <- rbind(ctrp1, ctrp2, gdsc1, gdsc2)
colnames(drug) <- gsub("\\.", "-", colnames(drug))

# Clean up
rm(ctrp1, ctrp2, gdsc1, gdsc2, pan_CTRP1_Log2, pan_CTRP2_Log2, 
   pan_GDSC1_Log2, pan_GDSC2_Log2)
gc()

# TE expression files
Dir <- paste0("/data/whu/home/te_database/results/1_merge_expression/filtered_", software)
files <- list.files(path = Dir, pattern = "_TE_fpkm_filtered.txt")
files <- files[files != "LAML_TE_fpkm_filtered.txt"]

# ==============================================================================
# Function: Calculate Partial Correlation
# ==============================================================================

# Calculate partial correlation for a single TE against all drugs
calc_partial_cor_single_te <- function(te_expr, drug_expr_matrix, covariates_matrix) {
  
  # Check if Gender and Plate have only one level
  gender_levels <- length(unique(covariates_matrix$gender))
  plate_levels <- length(unique(covariates_matrix$plate))
  
  # Prepare covariates
  covar_cols <- c("purity", "age_at_initial_pathologic_diagnosis")  # Always include purity and age
  
  # Add gender if varies
  if(gender_levels > 1) {
    covar_cols <- c(covar_cols, "gender")
  }
  
  # Extract basic covariates
  covariates_numeric <- as.matrix(covariates_matrix[, covar_cols])
  
  # Add plate dummies if needed
  if(plate_levels > 1) {
    # Create dummy variables (k-1 dummies for k levels)
    plate_dummies <- model.matrix(~ plate, data = covariates_matrix)
    # Remove intercept column
    plate_dummies <- plate_dummies[, -1, drop = FALSE]
    # Combine with other covariates
    covariates_numeric <- cbind(covariates_numeric, plate_dummies)
  }
  
  te_expr_vec <- as.numeric(te_expr)
  n_drugs <- nrow(drug_expr_matrix)
  
  results <- data.frame(
    Drug_id = rownames(drug_expr_matrix),
    cor = numeric(n_drugs),
    pvalue = numeric(n_drugs),
    stringsAsFactors = FALSE
  )
  
  # Calculate for each drug
  for (j in 1:n_drugs) {
    drug_expr_vec <- as.numeric(drug_expr_matrix[j, ])
    
    # Combine data: 1 TE + 1 drug + covariates
    data_subset <- cbind(te_expr_vec, drug_expr_vec, covariates_numeric)
    
    tryCatch({
      # Calculate partial correlation (spearman)
      pcor_result <- pcor(data_subset, method = "spearman")
      
      # Extract partial correlation between TE (column 1) and drug (column 2)
      results$cor[j] <- pcor_result$estimate[1, 2]
      results$pvalue[j] <- pcor_result$p.value[1, 2]
      
    }, error = function(e) {
      results$cor[j] <<- NA
      results$pvalue[j] <<- NA
    })
  }
  
  return(results)
}

# ==============================================================================
# Function: Batch Processing for Multiple TEs
# ==============================================================================

calc_partial_cor_batch <- function(te_expr_matrix, drug_expr_matrix, 
                                   covariates_matrix, te_chunk) {
  # Process a chunk of TEs
  
  results_list <- list()
  
  for (i in te_chunk) {
    te_name <- rownames(te_expr_matrix)[i]
    
    # Calculate correlation
    te_results <- calc_partial_cor_single_te(
      te_expr = te_expr_matrix[i, ],
      drug_expr_matrix = drug_expr_matrix,
      covariates_matrix = covariates_matrix
    )
    
    te_results$TE_id <- te_name
    results_list[[te_name]] <- te_results
  }
  
  return(do.call(rbind, results_list))
}

# ==============================================================================
# Main Analysis Loop: Process Each Cancer Type
# ==============================================================================

f <- args[1]
Cancer <- gsub("_TE_fpkm_filtered.txt", "", f)
cat("\n")
cat("==============================================================================\n")
cat(paste0("Processing: ", Cancer, "\n"))
cat("==============================================================================\n")

start_time <- Sys.time()

# ----------------------------------------------------------------------------
# Load Data
# ----------------------------------------------------------------------------

# Clinical data
pan_cli <- data.table::fread(
  "/data/whu/home/te_database/data/tcga/TCGA-CDR-SupplementalTableS1.txt", 
  header = T, data.table = F, check.names = F
)
pan_cli_c <- pan_cli[pan_cli$type == Cancer, c(1, 3:4)]

# TE element expression
exp_cancer <- data.table::fread(
  file.path(Dir, f), 
  header = T, data.table = F, check.names = F
)
rownames(exp_cancer) <- exp_cancer[, 1]
exp_cancer <- exp_cancer[, -1]
exp_cancer <- exp_cancer[, substr(colnames(exp_cancer), 14, 15) == "01"]  # primary tumors

# ----------------------------------------------------------------------------
# Sample Matching & Filtering
# ----------------------------------------------------------------------------

# Match samples between drug response and TE expression
shared_sample <- intersect(colnames(drug), colnames(exp_cancer))

exp_drug <- drug[, shared_sample]
exp_cancer <- exp_cancer[, shared_sample]

# ----------------------------------------------------------------------------
# Prepare Covariates
# ----------------------------------------------------------------------------

covariates <- data.frame(
  sample_id = shared_sample,
  patient_id = substr(shared_sample, 1, 12),
  plate = substr(shared_sample, 22, 25)
)

# Add age and gender
covariates <- left_join(covariates, pan_cli_c, 
                        by = c("patient_id" = "bcr_patient_barcode"))

# Add tumor purity
purity <- data.table::fread(
  "/data/whu/home/ATEC/data/batch/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
)
covariates$purity_id <- substr(covariates$sample_id, 1, 15)
covariates <- left_join(covariates, purity[, c(1, 4)], 
                        by = c("purity_id" = "array"))

# Remove missing values
covariates <- na.omit(covariates)

# Convert categorical variables to numeric
# Gender: MALE=1, FEMALE=0
covariates$gender <- ifelse(covariates$gender == "MALE", 1, 0)

# Plate: convert to factor
n_plates <- length(unique(covariates$plate))
cat(sprintf("Unique plates: %d\n", n_plates))
covariates$plate <- as.factor(covariates$plate)

# Filter data to samples with complete covariates
exp_drug <- exp_drug[, covariates$sample_id]
exp_cancer <- exp_cancer[, covariates$sample_id]

# ----------------------------------------------------------------------------
# Summary Statistics
# ----------------------------------------------------------------------------

cat(sprintf("Final samples: %d\n", ncol(exp_drug)))
cat(sprintf("Final TEs: %d\n", nrow(exp_cancer)))
cat(sprintf("Final drugs: %d\n", nrow(exp_drug)))

# ----------------------------------------------------------------------------
# Parallel Computation
# ----------------------------------------------------------------------------

total_TEs <- nrow(exp_cancer)

# Dynamic chunk sizing based on data size
optimal_chunk_size <- max(10, ceiling(total_TEs / max_cores))
te_chunks <- split(1:total_TEs, ceiling(seq_along(1:total_TEs) / optimal_chunk_size))

cat(sprintf("Processing %d chunks (size ~%d) with %d cores\n", 
            length(te_chunks), optimal_chunk_size, 
            min(length(te_chunks), max_cores)))

# Run parallel computation
result_list <- mclapply(
  te_chunks,
  function(chunk) {
    calc_partial_cor_batch(
      te_expr_matrix = exp_cancer,
      drug_expr_matrix = exp_drug,
      covariates_matrix = covariates,
      te_chunk = chunk
    )
  },
  mc.cores = min(length(te_chunks), max_cores)
)

# ----------------------------------------------------------------------------
# Combine & Process Results
# ----------------------------------------------------------------------------

cat("Merging results...\n")
rt_cor_total <- do.call(rbind, result_list)
rt_cor_total <- rt_cor_total[, c("TE_id", "Drug_id", "cor", "pvalue")]

# FDR correction
rt_cor_total$fdr <- p.adjust(rt_cor_total$pvalue, method = "fdr")

# Significance labels
rt_cor_total$significance <- ifelse(
  rt_cor_total$fdr < 0.05,
  ifelse(abs(rt_cor_total$cor) > 0.3, "High", "Low"),
  "NS"
)

# Time
end_time <- Sys.time()
cat(sprintf("Processing time: %.2f minutes\n", 
            as.numeric(difftime(end_time, start_time, units = "mins"))))

# Summary
cat(sprintf("Total correlations: %d\n", nrow(rt_cor_total)))
cat(sprintf("Significant (FDR < 0.05): %d (%.1f%%)\n", 
            sum(rt_cor_total$fdr < 0.05, na.rm = TRUE),
            100 * sum(rt_cor_total$fdr < 0.05, na.rm = TRUE) / nrow(rt_cor_total)))

# ----------------------------------------------------------------------------
# Export Results
# ----------------------------------------------------------------------------

# Significant results only
rt_cor_sig <- rt_cor_total[abs(rt_cor_total$cor) > 0.3 & rt_cor_total$fdr < 0.05, ]
output_file_sig <- paste0(Cancer, "_drug_TE_element_cor.txt")
write.table(rt_cor_sig, output_file_sig, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Significant results saved: %s\n", output_file_sig))

# Clean memory
rm(exp_cancer, exp_drug, covariates, rt_cor_total, result_list, 
   rt_cor_sig, te_chunks, shared_sample, pan_cli, pan_cli_c, purity)
gc()

cat(paste0("Completed: ", Cancer, "\n"))