# =============================================================================
# deseq2.R
# BIOL0041 - KCNA5 Pathway Analysis in TCGA-COAD
#
# Purpose: DESeq2 differential expression analysis - Tumour vs Normal
#          Check expression of genes of interest (KCNA5, FZD6)
#          Produce volcano plot
# Input:   data/cleaned/TCGA-COAD_counts_cleaned.rds
#          data/cleaned/tumour_samples.rds
#          data/cleaned/normal_samples.rds
# Output:  results/deseq2/DESeq2_tumour_vs_normal.csv
#          results/deseq2/volcano_tumour_vs_normal.pdf
#          results/deseq2/dds.rds (DESeq2 object for GSEA)
# =============================================================================

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
counts_raw     <- readRDS("data/cleaned/TCGA-COAD_counts_cleaned.rds")
tumour_samples <- readRDS("data/cleaned/tumour_samples.rds")
normal_samples <- readRDS("data/cleaned/normal_samples.rds")

cat("Counts matrix loaded:", dim(counts_raw), "\n")
cat("Tumour samples:", length(tumour_samples), "\n")
cat("Normal samples:", length(normal_samples), "\n")

# -----------------------------------------------------------------------------
# Build DESeq2 dataset
# -----------------------------------------------------------------------------

# Subset to tumour and normal samples only
counts_subset <- counts_raw[, c(tumour_samples, normal_samples)]

# Build metadata - Normal is reference level
metadata <- data.frame(
    sample    = c(tumour_samples, normal_samples),
    condition = factor(c(rep("Tumour", length(tumour_samples)),
                         rep("Normal", length(normal_samples))),
                       levels = c("Normal", "Tumour"))
)
rownames(metadata) <- metadata$sample

cat("\nSample condition counts:\n")
print(table(metadata$condition))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData   = metadata,
    design    = ~condition
)

# Filter lowly expressed genes
# Keep genes with at least 10 counts in at least 10 samples
keep <- rowSums(counts(dds) >= 10) >= 10
dds  <- dds[keep, ]
cat("\nGenes retained after filtering:", nrow(dds), "\n")

# -----------------------------------------------------------------------------
# Run DESeq2
# -----------------------------------------------------------------------------
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)

# Save DESeq2 object for use in GSEA script
saveRDS(dds, "results/deseq2/dds.rds")

# -----------------------------------------------------------------------------
# Extract results
# -----------------------------------------------------------------------------
res <- results(dds,
               contrast = c("condition", "Tumour", "Normal"),
               alpha    = 0.05)

cat("\nDESeq2 results summary:\n")
summary(res)

# Convert to dataframe
res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

cat("\nSignificant genes (padj < 0.05):", 
    sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in tumour:", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("Downregulated in tumour:", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE), "\n")

# Save results
write.csv(res_df, "results/deseq2/DESeq2_tumour_vs_normal.csv", 
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Check genes of interest
# -----------------------------------------------------------------------------
genes_of_interest <- c("KCNA5", "FZD6", "CTNNB1", "FZD1", "FZD2", 
                        "WNT5A", "KCNQ3", "KCNQ1","FZD8")

cat("\nGenes of interest:\n")
print(res_df[res_df$gene %in% genes_of_interest, 
             c("gene", "log2FoldChange", "pvalue", "padj")])

# -----------------------------------------------------------------------------
# Volcano plot
# -----------------------------------------------------------------------------
# Genes to highlight
highlight_genes <- c("KCNA5")

volcano_df <- res_df %>%
    mutate(
        significance = case_when(
            padj < 0.05 & log2FoldChange >  1 ~ "Up in tumour",
            padj < 0.05 & log2FoldChange < -1 ~ "Down in tumour",
            TRUE ~ "NS"
        ),
        label = ifelse(gene %in% highlight_genes, gene, "")
    )

p_volcano <- ggplot(volcano_df, 
                    aes(x      = log2FoldChange,
                        y      = -log10(padj),
                        colour = significance,
                        label  = label)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_point(data    = volcano_df %>% filter(gene %in% highlight_genes),
               colour  = "black",
               size    = 3) +
    geom_text_repel(colour            = "black",
                    fontface          = "bold",
                    max.overlaps      = Inf,
                    min.segment.length = 0) +
    scale_colour_manual(values = c("Up in tumour"   = "red",
                                   "Down in tumour" = "blue",
                                   "NS"             = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(x        = "log2 Fold Change (Tumour vs Normal)",
         y        = "-log10(adjusted p-value)",
         title    = "Tumour vs Normal DESeq2 results",
         subtitle = "TCGA-COAD | KCNA5 and FZD6 highlighted",
         colour   = "") +
    theme_bw()

print(p_volcano)
ggsave("results/deseq2/volcano_tumour_vs_normal.pdf", 
       plot  = p_volcano,
       width = 8, 
       height = 6)
cat("Saved: results/deseq2/volcano_tumour_vs_normal.pdf\n")
