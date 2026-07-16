# =============================================================================
# gsea_nonT1mutant.R
# BIOL0041 - KCNA5 Pathway Analysis in TCGA-COAD
#
# Purpose: GSEA on DESeq2 results for KCNA5 non-T1 domain mutants vs WT
#          Gene sets: Hallmark
# Input:   results/deseq2/DESeq2_KCNA5_nonT1Mutant_vs_WT.csv
# Output:  results/gsea/nonT1mutant/GSEA_hallmark.pdf
#          results/gsea/nonT1mutant/GSEA_hallmark_all.csv
#          results/gsea/nonT1mutant/GSEA_hallmark_significant.csv
# =============================================================================

library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)

# -----------------------------------------------------------------------------
# Load DESeq2 results
# -----------------------------------------------------------------------------
res_df <- read.csv("results/deseq2/DESeq2_KCNA5_nonT1Mutant_vs_WT.csv")

cat("DESeq2 results loaded:", nrow(res_df), "genes\n")
cat("Significant genes (padj < 0.05):",
    sum(res_df$padj < 0.05, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# Build ranked gene list
# -----------------------------------------------------------------------------
gene_ranks <- res_df %>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange)) %>%
    mutate(
        pvalue_adj = ifelse(pvalue == 0, .Machine$double.xmin, pvalue),
        rank       = sign(log2FoldChange) * -log10(pvalue_adj)
    ) %>%
    arrange(desc(rank)) %>%
    { setNames(.$rank, .$gene) }

cat("Ranked gene list:", length(gene_ranks), "genes\n")
cat("Top 5:\n");    print(head(gene_ranks, 5))
cat("Bottom 5:\n"); print(tail(gene_ranks, 5))

# -----------------------------------------------------------------------------
# Get gene sets
# -----------------------------------------------------------------------------
cat("\nDownloading gene sets...\n")

hallmark <- msigdbr(species    = "Homo sapiens",
                    collection = "H") %>%
    split(x = .$gene_symbol, f = .$gs_name)

cat("Hallmark gene sets:", length(hallmark), "\n")

# -----------------------------------------------------------------------------
# Helper function: run GSEA, save and plot results
# (same as 05b but with T1 mutant subtitle)
# -----------------------------------------------------------------------------
run_gsea <- function(gene_ranks,
                      gene_sets,
                      label,
                      outdir,
                      n_top    = 20,
                      subtitle = "TCGA-COAD | Non-T1 mutant (residues < 128 or > 216) vs WT") {
    
    cat("\n=== GSEA:", label, "===\n")
    
    set.seed(42)
    gsea_results <- fgsea(pathways    = gene_sets,
                           stats       = gene_ranks,
                           eps         = 0,
                           nPermSimple = 10000)
    
    gsea_sig <- gsea_results %>%
        filter(padj < 0.05) %>%
        arrange(desc(NES))
    
    cat("Significant pathways (padj < 0.05):", nrow(gsea_sig), "\n")
    
    # Save all results
    write.csv(gsea_results %>%
                  as.data.frame() %>%
                  dplyr::select(-leadingEdge) %>%
                  arrange(pval),
              file.path(outdir, paste0("GSEA_", label, "_all.csv")),
              row.names = FALSE)
    
    # Save significant results
    write.csv(gsea_sig %>%
                  as.data.frame() %>%
                  dplyr::select(-leadingEdge),
              file.path(outdir, paste0("GSEA_", label, "_significant.csv")),
              row.names = FALSE)
    
    # Plot
    if (nrow(gsea_sig) > 0) {
        
        plot_df <- gsea_sig %>%
            as.data.frame() %>%
            arrange(NES) %>%
            head(n_top) %>%
            mutate(
                pathway_clean = gsub(paste0("HALLMARK_"),
                                      "", pathway),
                pathway_clean = gsub("_", " ", pathway_clean),
                direction     = ifelse(NES > 0,
                                       "Up in T1 mutant",
                                       "Down in T1 mutant")
            )
        
        p <- ggplot(plot_df,
                    aes(x    = reorder(pathway_clean, NES),
                        y    = NES,
                        fill = direction)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            scale_fill_manual(values = c("Up in T1 mutant"   = "red",
                                          "Down in T1 mutant" = "blue")) +
            geom_hline(yintercept = 0, linetype = "solid") +
            labs(x        = "",
                 y        = "Normalised Enrichment Score (NES)",
                 fill     = "",
                 title    = paste("GSEA:", label, "pathways"),
                 subtitle = subtitle) +
            theme_bw() +
            theme(axis.text.y     = element_text(size = 8),
                  legend.position = "bottom")
        
        print(p)
        ggsave(file.path(outdir, paste0("GSEA_", label, ".pdf")),
               plot   = p,
               width  = 12,
               height = max(6, nrow(plot_df) * 0.4))
        cat("Saved:", file.path(outdir, paste0("GSEA_", label, ".pdf")), "\n")
        
    } else {
        cat("No significant pathways found for", label, "\n")
    }
    
    invisible(gsea_results)
}

# -----------------------------------------------------------------------------
# Run GSEA for all three gene set collections
# -----------------------------------------------------------------------------
gsea_hallmark <- run_gsea(
    gene_ranks = gene_ranks,
    gene_sets  = hallmark,
    label      = "hallmark",
    outdir     = "results/gsea/nonT1mutant"
)
