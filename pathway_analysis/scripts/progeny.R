# =============================================================================
# progeny.R
# KCNA5 Pathway Analysis in TCGA-COAD
#
# Purpose: Run PROGENY pathway scoring and correlate with gene expression
# Input:   data/cleaned/TCGA-COAD_tpm_cleaned.rds
# Output:  results/progeny/*.png
#          results/progeny/progeny_scores.rds
# =============================================================================

library(progeny)
library(dplyr)
library(tibble)
library(ggplot2)

# -----------------------------------------------------------------------------
# Load cleaned TPM data
# -----------------------------------------------------------------------------
tpm <- readRDS("data/cleaned/TCGA-COAD_tpm_cleaned.rds")
cat("TPM matrix loaded:", dim(tpm), "\n")

# -----------------------------------------------------------------------------
# Run PROGENY pathway scoring
# -----------------------------------------------------------------------------
cat("Running PROGENY...\n")
progeny_scores <- progeny(
    as.matrix(tpm[]),
    scale      = TRUE,
    organism   = "Human",
    top        = 500,
    perm       = 1
)

cat("PROGENY scores dimensions:", dim(progeny_scores), "\n")
# Expected: samples x 14 pathways

# Save scores for reuse in other scripts
saveRDS(progeny_scores, "results/progeny/progeny_scores.rds")

# -----------------------------------------------------------------------------
# Helper function: correlate any gene with all PROGENY pathway scores
# Produces a bar plot coloured by correlation direction
# -----------------------------------------------------------------------------
plot_progeny_correlation <- function(gene,
                                      tpm,
                                      progeny_scores,
                                      method      = "pearson",
                                      p_threshold = 0.05,
                                      save        = TRUE,
                                      outdir      = "results/progeny") {
    
    # Check gene exists
    if (!gene %in% rownames(tpm)) {
        stop(paste("Gene", gene, "not found in expression matrix"))
    }
    
    # Extract and align gene expression to progeny sample order
    gene_expr <- as.numeric(tpm[gene, ])
    names(gene_expr) <- colnames(tpm)
    gene_expr <- gene_expr[rownames(progeny_scores)]
    
    # Correlate with each pathway
    cor_results <- apply(progeny_scores, 2, function(pathway) {
        test <- cor.test(gene_expr, pathway, method = method)
        c(correlation = unname(test$estimate),
          pvalue      = test$p.value)
    })
    
    # Build results dataframe
    cor_df <- as.data.frame(t(cor_results)) %>%
        rownames_to_column("Pathway") %>%
        arrange(pvalue)
    
    # Plot
    p <- ggplot(cor_df, aes(x     = reorder(Pathway, -log10(pvalue)),
                             y     = -log10(pvalue),
                             fill  = correlation)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_gradient2(low     = "blue",
                             mid     = "white",
                             high    = "red",
                             midpoint = 0) +
        geom_hline(yintercept = -log10(p_threshold),
                   linetype   = "dashed",
                   colour     = "black") +
        labs(x        = "Pathway",
             y        = "-log10(p-value)",
             fill     = "Pearson r",
             title    = paste(gene, "correlation with PROGENY pathway scores"),
             subtitle = "TCGA-COAD") +
        theme_bw()
    
    print(p)
    
    # Save figure
    if (save) {
        outfile <- file.path(outdir, paste0(gene, "_progeny_correlation_normal.png"))
        ggsave(outfile, plot = p, width = 6, height = 8)
        cat("Saved:", outfile, "\n")
    }
    
    invisible(cor_df)
}

# -----------------------------------------------------------------------------
# Run correlations
# -----------------------------------------------------------------------------

# KCNA5 - gene of interest
kcna5_progeny <- plot_progeny_correlation("KCNA5", tpm, progeny_scores)


# Save correlation tables
write.csv(kcna5_progeny, "results/progeny/KCNA5_progeny_correlations_normal.csv", 
          row.names = FALSE)

