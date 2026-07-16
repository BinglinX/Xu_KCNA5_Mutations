# =============================================================================
# data_cleaning.R
# 
# Purpose: Load, clean and save TPM and raw counts matrices from GDC
# Input:   TCGA-COAD.star_tpm.tsv
#          TCGA-COAD.star_counts.tsv
# Output:  data/cleaned/TCGA-COAD_tpm_cleaned.rds
#          data/cleaned/TCGA-COAD_counts_cleaned.rds
# =============================================================================

library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# -----------------------------------------------------------------------------
# Helper function: clean a GDC expression matrix
# Works for both TPM and counts files (same structure)
# log2_transformed: set TRUE if values are log2(x+1) transformed
# is_counts: set TRUE to round to integers (required for DESeq2)
# -----------------------------------------------------------------------------

clean_gdc_matrix <- function(filepath, 
                              log2_transformed = TRUE,
                              is_counts = FALSE) {
    
    cat("Loading:", filepath, "\n")
    mat <- read.delim(filepath, check.names = FALSE, header = TRUE)
    cat("Dimensions after loading:", dim(mat), "\n")
    
    # Back-transform if log2 transformed
    if (log2_transformed) {
        cat("Back-transforming log2 values...\n")
        expr_mat <- as.data.frame(2^as.matrix(mat[, -1]) - 1)
    } else {
        expr_mat <- as.data.frame(as.matrix(mat[, -1]))
    }
    
    # Round to integers if counts
    if (is_counts) {
        expr_mat <- round(expr_mat)
    }
    
    # Add Ensembl IDs (strip version numbers)
    expr_mat$Ensembl_ID <- sub("\\..*", "", mat$Ensembl_ID)
    
    # Map Ensembl IDs to Hugo symbols
    cat("Mapping Ensembl IDs to Hugo symbols...\n")
    symbols <- mapIds(org.Hs.eg.db,
                      keys      = expr_mat$Ensembl_ID,
                      column    = "SYMBOL",
                      keytype   = "ENSEMBL",
                      multiVals = "first")
    
    expr_mat$Hugo_Symbol <- symbols
    
    # Remove unmapped genes
    n_before <- nrow(expr_mat)
    expr_mat <- expr_mat[!is.na(expr_mat$Hugo_Symbol), ]
    cat("Removed", n_before - nrow(expr_mat), "unmapped genes\n")
    
    # Remove Ensembl ID column
    expr_mat <- expr_mat[, -which(colnames(expr_mat) == "Ensembl_ID")]
    
    # Move Hugo_Symbol to front
    expr_mat <- expr_mat[, c("Hugo_Symbol",
                              setdiff(colnames(expr_mat), "Hugo_Symbol"))]
    
    # Handle duplicates - keep highest mean expressed row
    n_dups <- sum(duplicated(expr_mat$Hugo_Symbol))
    cat("Removing", n_dups, "duplicate gene symbols...\n")
    
    expr_mat <- expr_mat %>%
        mutate(mean_expr = rowMeans(dplyr::select(., -Hugo_Symbol), 
                                    na.rm = TRUE)) %>%
        group_by(Hugo_Symbol) %>%
        slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        dplyr::select(-mean_expr)
    
    # Set row names
    expr_mat <- as.data.frame(expr_mat)
    rownames(expr_mat) <- expr_mat$Hugo_Symbol
    expr_mat <- expr_mat[, -1]
    
    return(expr_mat)
}

# TPM - keep as log2 (already transformed by GDC)
tpm_log2 <- clean_gdc_matrix(
    filepath         = "data/raw/TCGA-COAD.star_tpm.tsv",
    log2_transformed = FALSE,  # don't back-transform
    is_counts        = FALSE
)
saveRDS(tpm_log2, "data/cleaned/TCGA-COAD_tpm_cleaned.rds")


# Raw counts
counts_raw <- clean_gdc_matrix(
    filepath         = "data/raw/TCGA-COAD.star_counts.tsv",
    log2_transformed = TRUE,
    is_counts        = TRUE
)

# Save
saveRDS(counts_raw, "data/cleaned/TCGA-COAD_counts_cleaned.rds")

# -----------------------------------------------------------------------------
# Identify tumour and normal samples
# TCGA barcode convention:
#   -01 = primary solid tumour
#   -11 = solid tissue normal
# -----------------------------------------------------------------------------
tumour_samples <- colnames(counts_raw)[grepl("-01", colnames(counts_raw))]
normal_samples <- colnames(counts_raw)[grepl("-11", colnames(counts_raw))]

# Save sample lists
saveRDS(tumour_samples, "data/cleaned/tumour_samples.rds")
saveRDS(normal_samples, "data/cleaned/normal_samples.rds")
