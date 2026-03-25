# ==============================================================================
# Multiple Linear Regression: TE Element vs Regulator Genes
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
args=commandArgs(T)

# ==============================================================================
# Configuration
# ==============================================================================

software <- "SQuIRE"
max_cores <- 60

# Set working directory
setwd("/data/whu/home/ATEC/resutls/2_regressor_TE_regression_SQuIRE")

# Create output directory if not exists
if (!dir.exists(getwd())) {
  dir.create(getwd(), recursive = TRUE)
}

# ==============================================================================
# Load Reference Data
# ==============================================================================

# TE expression files
Dir <- paste0("/data/whu/home/te_database/results/1_merge_expression/filtered_", software)
files <- list.files(path = Dir, pattern = "_TE_fpkm_filtered.txt")

# Regulator gene lists
TF <- fread("/data/whu/home/te_database/data/TF list.txt", 
            header = TRUE, data.table = FALSE, check.names = FALSE)
RBP <- fread("/data/whu/home/te_database/data/RBP_union_list.txt", 
             header = TRUE, data.table = FALSE, check.names = FALSE)
regulator <- union(TF$gene_name, RBP$gene_name)

# Load purity data
purity <- fread("/data/whu/home/ATEC/data/batch/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", 
                header = TRUE, data.table = FALSE, check.names = FALSE)
purity_dt <- purity[, c("array", "purity")]
rm(purity)

# Load clinical data
pan_cli <- fread("/data/whu/home/te_database/data/tcga/TCGA-CDR-SupplementalTableS1.txt", 
                 header = TRUE, data.table = FALSE, check.names = FALSE)

# ==============================================================================
# Regression Function
# ==============================================================================

calc_regression_vectorized <- function(te_expr_matrix, gene_expr_matrix, covariates_df) {
  
  # Pre-allocate results matrix
  n_tes <- nrow(te_expr_matrix)
  n_genes <- nrow(gene_expr_matrix)
  n_total <- n_tes * n_genes
  n_samples <- ncol(te_expr_matrix)
  
  # Check if Gender and Plate have only one level
  gender_levels <- length(unique(covariates_df$gender))
  plate_levels <- length(unique(covariates_df$plate))
  
  # ========== Dynamic Model Construction ==========
  # Build model based on available covariates
  formula_parts <- c("purity", "age")
  
  if(gender_levels > 1) {
    formula_parts <- c(formula_parts, "gender")
  }
  
  if(plate_levels > 1) {
    formula_parts <- c(formula_parts, "plate")
  }
  
  # Build formula
  formula_str <- paste("~", paste(formula_parts, collapse = " + "))
  
  X_cov <- model.matrix(as.formula(formula_str), data = covariates_df)
  
  # Initialize results storage
  results <- data.table(
    TE_id = character(n_total),
    Gene_id = character(n_total),
    beta = numeric(n_total),
    beta_se = numeric(n_total),
    beta_pvalue = numeric(n_total),
    n_samples = integer(n_total)
  )
  
  idx <- 1
  
  # Process each TE
  for (i in 1:n_tes) {
    te_name <- rownames(te_expr_matrix)[i]
    y <- as.numeric(te_expr_matrix[i, ])
    
    # Process all genes for this TE at once (vectorized)
    for (j in 1:n_genes) {
      gene_name <- rownames(gene_expr_matrix)[j]
      x_gene <- as.numeric(gene_expr_matrix[j, ])
      
      # Build full design matrix
      X <- cbind(x_gene, X_cov[, -1])  # Remove intercept duplicate
      X <- cbind(1, X)  # Add intercept
      
      tryCatch({
        # Fast matrix computation
        XtX <- crossprod(X)
        Xty <- crossprod(X, y)
        
        # Solve for coefficients
        beta_all <- solve(XtX, Xty)
        
        # Residuals and variance
        y_hat <- X %*% beta_all
        residuals <- y - y_hat
        sigma2 <- sum(residuals^2) / (n_samples - ncol(X))
        
        # Standard errors
        vcov <- sigma2 * solve(XtX)
        se_all <- sqrt(diag(vcov))
        
        # T-statistics and p-values
        t_stats <- beta_all / se_all
        p_values <- 2 * pt(-abs(t_stats), df = n_samples - ncol(X))
        
        # Store results for gene coefficient (index 2, after intercept)
        n_samples_value <- n_samples
        results[idx, `:=`(
          TE_id = te_name,
          Gene_id = gene_name,
          beta = beta_all[2],
          beta_se = se_all[2],
          beta_pvalue = p_values[2],
          n_samples = n_samples_value
        )]
        
      }, error = function(e) {
        results[idx, `:=`(
          TE_id = te_name,
          Gene_id = gene_name,
          beta = NA_real_,
          beta_se = NA_real_,
          beta_pvalue = NA_real_,
          n_samples = n_samples_value
        )]
      })
      
      idx <- idx + 1
    }
  }
  return(as.data.frame(results))
}

# ==============================================================================
# Batch Processing
# ==============================================================================

calc_regression_batch_optimized <- function(te_expr_matrix, gene_expr_matrix, 
                                            covariates_df, te_chunk) {
  
  # Process only the specified chunk
  te_subset <- te_expr_matrix[te_chunk, , drop = FALSE]
  
  # Use vectorized function
  results <- calc_regression_vectorized(
    te_expr_matrix = te_subset,
    gene_expr_matrix = gene_expr_matrix,
    covariates_df = covariates_df
  )
  
  return(results)
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

# Clinical data for this cancer type
pan_cli_c <- pan_cli[pan_cli$type == Cancer, c(1, 3:4)]

# TE element expression
exp_cancer <- fread(
  file.path(Dir, f), 
  header = TRUE, data.table = FALSE, check.names = FALSE
)
rownames(exp_cancer) <- exp_cancer[, 1]
exp_cancer <- exp_cancer[, -1]
exp_cancer <- exp_cancer[, substr(colnames(exp_cancer), 14, 15) == "01"]
exp_cancer <- log2(exp_cancer + 1) # log2 transformation

# Gene expression
exp_PC <- fread(
  paste0("/data/whu/home/te_database/results/0_TCGA_geneExp/", Cancer, "_geneExp.txt"), 
  header = TRUE, data.table = FALSE, check.names = FALSE
)
exp_PC[, 1] <- gsub(".*\\|", "", exp_PC[, 1])
exp_PC <- exp_PC[!duplicated(exp_PC[, 1]), ]
rownames(exp_PC) <- exp_PC[, 1]
exp_PC <- exp_PC[, -1]
exp_PC <- exp_PC[rownames(exp_PC) %in% regulator, ]

# ----------------------------------------------------------------------------
# Sample Matching & Filtering
# ----------------------------------------------------------------------------

shared_sample <- intersect(colnames(exp_PC), colnames(exp_cancer))
exp_PC <- exp_PC[, shared_sample, drop = FALSE]
exp_cancer <- exp_cancer[, shared_sample, drop = FALSE]

# ----------------------------------------------------------------------------
# Prepare Covariates
# ----------------------------------------------------------------------------

covariates <- data.frame(
  sample_id = shared_sample,
  patient_id = substr(shared_sample, 1, 12),
  plate = substr(shared_sample, 22, 25),
  purity_id = substr(shared_sample, 1, 15)
)

# Join operations
covariates <- merge(covariates, pan_cli_c, 
                    by.x = "patient_id", by.y = "bcr_patient_barcode", 
                    all.x = TRUE)
colnames(covariates)[5] <- "age"
covariates <- merge(covariates, purity_dt, 
                    by.x = "purity_id", by.y = "array", 
                    all.x = TRUE)

# Remove missing values
covariates <- na.omit(covariates)

# Convert categorical variables
covariates$gender <- ifelse(covariates$gender == "MALE", 1, 0)

# Plate info
n_plates <- length(unique(covariates$plate))
cat(sprintf("Unique plates: %d\n", n_plates))

covariates$plate <- factor(covariates$plate)

# Check if Gender and Plate have only one level
gender_levels <- length(unique(covariates$gender))
plate_levels <- length(unique(covariates$plate))

if(gender_levels == 1) {
  cat("WARNING: Gender has only one level - will be EXCLUDED from model\n")
}
if(plate_levels == 1) {
  cat("WARNING: Plate has only one level - will be EXCLUDED from model\n")
}

# Model formula preview
formula_parts <- c("purity", "age")
if(gender_levels > 1) formula_parts <- c(formula_parts, "gender")
if(plate_levels > 1) formula_parts <- c(formula_parts, "plate")
cat(sprintf("Model formula: TE ~ Gene + %s\n", paste(formula_parts, collapse = " + ")))

# Update expression matrices
exp_PC <- exp_PC[, covariates$sample_id, drop = FALSE]
exp_cancer <- exp_cancer[, covariates$sample_id, drop = FALSE]

# Filter low expression genes (vectorized)
gene_means <- rowMeans(exp_PC)
exp_PC <- exp_PC[gene_means > 1, , drop = FALSE]
exp_PC <- log2(exp_PC + 1) # log2 transformation

# ----------------------------------------------------------------------------
# Summary Statistics
# ----------------------------------------------------------------------------

cat(sprintf("Final samples: %d\n", ncol(exp_PC)))
cat(sprintf("Final TEs: %d\n", nrow(exp_cancer)))
cat(sprintf("Final genes: %d\n", nrow(exp_PC)))

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
    calc_regression_batch_optimized(
      te_expr_matrix = exp_cancer,
      gene_expr_matrix = exp_PC,
      covariates_df = covariates,
      te_chunk = chunk
    )
  },
  mc.cores = min(length(te_chunks), max_cores),
  mc.preschedule = FALSE  # Better load balancing
)

# ----------------------------------------------------------------------------
# Combine & Process Results
# ----------------------------------------------------------------------------

cat("Merging results...\n")
rt_reg_total <- rbindlist(result_list)
setDF(rt_reg_total)  # Convert back to data.frame

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

rt_reg_sig <- rt_reg_total[rt_reg_total$fdr < 0.05 & !is.na(rt_reg_total$fdr), ]
setorder(setDT(rt_reg_sig), fdr)  # Fast sorting with data.table

output_file_sig <- paste0(Cancer, "_regulator_TE_element_regression_sig.txt")
fwrite(rt_reg_sig, output_file_sig, sep = "\t", quote = FALSE)
cat(sprintf("\nSignificant results saved: %s (%d associations)\n", 
            output_file_sig, nrow(rt_reg_sig)))

# Clean memory
rm(exp_cancer, exp_PC, covariates, rt_reg_total, result_list, 
   pan_cli_c, rt_reg_sig, sig_results)
gc()

cat(paste0("\nCompleted: ", Cancer, "\n"))


# ==============================================================================
# End of Analysis
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("All cancer types processed successfully!\n")
cat("==============================================================================\n")