# ==============================================================================
# Multiple Linear Regression: TE Element vs Drug Response
# Control for: Age, Gender, Plate and Tumor Purity
# ==============================================================================

# Clean environment
rm(list = ls())

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(parallel)
})

args <- commandArgs(T)

# ==============================================================================
# Configuration
# ==============================================================================

software <- "SQuIRE"
max_cores <- 60

# Set working directory
setwd("/data/whu/home/ATEC/results/4_drug_TE_regression_SQuIRE")

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

# ==============================================================================
# Function: Multiple Linear Regression (Drug ~ TE + Covariates)
# Model: Drug response ~ TE expression + Age + Gender + Plate + Purity
# ==============================================================================

calc_regression_single_te <- function(te_expr, drug_expr_matrix,
                                      X_cov, n_samples) {
  
  te_expr_vec <- as.numeric(te_expr)
  n_drugs <- nrow(drug_expr_matrix)
  
  results <- data.frame(
    Drug_id    = rownames(drug_expr_matrix),
    beta       = numeric(n_drugs),
    beta_se    = numeric(n_drugs),
    beta_pvalue = numeric(n_drugs),
    n_samples  = integer(n_drugs),
    stringsAsFactors = FALSE
  )
  
  # Build design matrix: intercept + TE + covariates (remove intercept from X_cov)
  X <- cbind(1, te_expr_vec, X_cov[, -1])
  
  # Pre-compute TE-side matrices (shared across all drugs)
  XtX <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX)) {
    results$beta        <- NA
    results$beta_se     <- NA
    results$beta_pvalue <- NA
    results$n_samples   <- n_samples
    return(results)
  }
  
  df_resid <- n_samples - ncol(X)
  
  for (j in 1:n_drugs) {
    y <- as.numeric(drug_expr_matrix[j, ])
    
    tryCatch({
      Xty       <- crossprod(X, y)
      beta_all  <- XtX %*% Xty
      
      # Residuals and variance
      residuals <- y - X %*% beta_all
      sigma2    <- sum(residuals^2) / df_resid
      
      # Standard errors
      se_all <- sqrt(diag(sigma2 * XtX))
      
      # T-statistics and p-values
      t_stat  <- beta_all[2] / se_all[2]
      p_value <- 2 * pt(-abs(t_stat), df = df_resid)
      
      results$beta[j]        <- beta_all[2]
      results$beta_se[j]     <- se_all[2]
      results$beta_pvalue[j] <- p_value
      results$n_samples[j]   <- n_samples
      
    }, error = function(e) {
      results$beta[j]        <<- NA
      results$beta_se[j]     <<- NA
      results$beta_pvalue[j] <<- NA
      results$n_samples[j]   <<- n_samples
    })
  }
  
  return(results)
}

# ==============================================================================
# Function: Batch Processing for Multiple TEs
# ==============================================================================

calc_regression_batch <- function(te_expr_matrix, drug_expr_matrix, 
                                  te_chunk, X_cov, n_samples) {
  
  results_list <- list()
  
  for (i in te_chunk) {
    te_name <- rownames(te_expr_matrix)[i]
    
    te_results <- calc_regression_single_te(
      te_expr        = te_expr_matrix[i, ],
      drug_expr_matrix = drug_expr_matrix,
      X_cov          = X_cov,
      n_samples      = n_samples
    )
    
    te_results$TE_id <- te_name
    results_list[[te_name]] <- te_results
  }
  
  return(do.call(rbind, results_list))
}

# ==============================================================================
# Main Analysis Loop: Process Each Cancer Type
# ==============================================================================

f      <- args[1]
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
exp_cancer <- log2(exp_cancer + 1)  # log2 transformation

# ----------------------------------------------------------------------------
# Sample Matching & Filtering
# ----------------------------------------------------------------------------

shared_sample <- intersect(colnames(drug), colnames(exp_cancer))

if (length(shared_sample) < 10) {
  cat(sprintf("WARNING: Only %d shared samples found. Skipping %s\n",
              length(shared_sample), Cancer))
  quit(save = "no")
}

exp_drug   <- drug[, shared_sample]
exp_cancer <- exp_cancer[, shared_sample]

# ----------------------------------------------------------------------------
# Prepare Covariates
# ----------------------------------------------------------------------------

covariates <- data.frame(
  sample_id  = shared_sample,
  patient_id = substr(shared_sample, 1, 12),
  plate      = substr(shared_sample, 22, 25)
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

# Plate: convert to factor
gender_levels <- length(unique(covariates$gender))
plate_levels  <- length(unique(covariates$plate))
cat(sprintf("Unique plates: %d\n", plate_levels))

if (gender_levels == 1) cat("WARNING: Gender has only one level - will be EXCLUDED from model\n")
if (plate_levels  == 1) cat("WARNING: Plate has only one level - will be EXCLUDED from model\n")

# Convert categorical variables
covariates$gender <- ifelse(covariates$gender == "MALE", 1, 0)
covariates$plate <- as.factor(covariates$plate)

# Filter data to samples with complete covariates
exp_drug   <- exp_drug[, covariates$sample_id]
exp_cancer <- exp_cancer[, covariates$sample_id]

# ----------------------------------------------------------------------------
# Pre-build Covariate Design Matrix (shared across all TE-drug pairs)
# ----------------------------------------------------------------------------

formula_parts <- c("purity", "age_at_initial_pathologic_diagnosis")
if (gender_levels > 1) formula_parts <- c(formula_parts, "gender")
if (plate_levels  > 1) formula_parts <- c(formula_parts, "plate")

formula_str <- paste("~", paste(formula_parts, collapse = " + "))
cat(sprintf("Model formula: Drug ~ TE + %s\n", paste(formula_parts, collapse = " + ")))

X_cov     <- model.matrix(as.formula(formula_str), data = covariates)
n_samples <- nrow(covariates)

# ----------------------------------------------------------------------------
# Summary Statistics
# ----------------------------------------------------------------------------

cat(sprintf("Final samples: %d\n", n_samples))
cat(sprintf("Final TEs: %d\n",    nrow(exp_cancer)))
cat(sprintf("Final drugs: %d\n",  nrow(exp_drug)))

# ----------------------------------------------------------------------------
# Parallel Computation
# ----------------------------------------------------------------------------

total_TEs <- nrow(exp_cancer)

optimal_chunk_size <- max(10, ceiling(total_TEs / max_cores))
te_chunks <- split(1:total_TEs, ceiling(seq_along(1:total_TEs) / optimal_chunk_size))

cat(sprintf("Processing %d chunks (size ~%d) with %d cores\n",
            length(te_chunks), optimal_chunk_size,
            min(length(te_chunks), max_cores)))

result_list <- mclapply(
  te_chunks,
  function(chunk) {
    calc_regression_batch(
      te_expr_matrix   = exp_cancer,
      drug_expr_matrix = exp_drug,
      te_chunk         = chunk,
      X_cov            = X_cov,
      n_samples        = n_samples
    )
  },
  mc.cores = min(length(te_chunks), max_cores)
)

# ----------------------------------------------------------------------------
# Combine & Process Results
# ----------------------------------------------------------------------------

cat("Merging results...\n")
rt_reg_total <- do.call(rbind, result_list)
rt_reg_total <- rt_reg_total[, c("TE_id", "Drug_id", "beta", "beta_se",
                                 "beta_pvalue", "n_samples")]

# FDR correction
rt_reg_total$fdr <- p.adjust(rt_reg_total$beta_pvalue, method = "fdr")


# Time
end_time <- Sys.time()
cat(sprintf("Processing time: %.2f minutes\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))

# ----------------------------------------------------------------------------
# Summary Statistics
# ----------------------------------------------------------------------------

cat("\n--- Summary Statistics ---\n")
cat(sprintf("Total associations tested: %d\n", nrow(rt_reg_total)))

n_sig <- sum(rt_reg_total$fdr < 0.05, na.rm = TRUE)
cat(sprintf("Significant (FDR < 0.05): %d (%.2f%%)\n",
            n_sig, 100 * n_sig / nrow(rt_reg_total)))

sig_results <- rt_reg_total[rt_reg_total$fdr < 0.05 & !is.na(rt_reg_total$fdr), ]

n_failed <- sum(is.na(rt_reg_total$beta))
if (n_failed > 0) {
  cat(sprintf("\nFailed regressions: %d (%.2f%%)\n",
              n_failed, 100 * n_failed / nrow(rt_reg_total)))
}

# ----------------------------------------------------------------------------
# Export Results
# ----------------------------------------------------------------------------

rt_reg_sig  <- rt_reg_total[rt_reg_total$fdr < 0.05 & !is.na(rt_reg_total$fdr), ]
output_file <- paste0(Cancer, "_drug_TE_element_regression_sig.txt")
write.table(rt_reg_sig, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nSignificant results saved: %s (%d associations)\n",
            output_file, nrow(rt_reg_sig)))

# Clean memory
rm(exp_cancer, exp_drug, covariates, rt_reg_total, result_list,
   rt_reg_sig, sig_results, te_chunks, shared_sample,
   pan_cli, pan_cli_c, purity, X_cov)
gc()

cat(paste0("\nCompleted: ", Cancer, "\n"))