# =============================================================================
# gsea.R
# BIOL0041 - KCNA5 Pathway Analysis in TCGA-COAD
#
# Purpose: GSEA on DESeq2 results for tumour vs. normal
#          Gene sets: Hallmark, GO:BP
# Input:   results/deseq2/DESeq2_tumour_vs_normal.csv
# Output:  results/gsea/GSEA_hallmark.pdf
# =============================================================================

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Load DESeq2 results
# -----------------------------------------------------------------------------
res_df <- read.csv("results/deseq2/DESeq2_tumour_vs_normal.csv")

# Build ranked gene list - replace Inf with max finite value
gene_ranks <- res_df %>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange)) %>%
    mutate(
        pvalue_adj = ifelse(pvalue == 0, .Machine$double.xmin, pvalue),
        rank = sign(log2FoldChange) * -log10(pvalue_adj)
    ) %>%
    arrange(desc(rank)) %>%
    { setNames(.$rank, .$gene) }

# Check no infinite values remain
cat("Infinite values:", sum(is.infinite(gene_ranks)), "\n")
cat("Top 5:\n"); print(head(gene_ranks, 5))

# Get gene sets - updated msigdbr syntax
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    split(x = .$gene_symbol, f = .$gs_name)

gobp <- msigdbr(species = "Homo sapiens", 
                collection = "C5", 
                subcollection = "GO:BP") %>%
    split(x = .$gene_symbol, f = .$gs_name)

# Run GSEA - remove nperm to use recommended fgseaMultilevel
set.seed(42)
gsea_hallmark <- fgsea(pathways = hallmark,
                        stats    = gene_ranks)

gsea_hallmark_sig <- gsea_hallmark %>%
    filter(padj < 0.05) %>%
    arrange(desc(NES))

cat("Significant Hallmark pathways:", nrow(gsea_hallmark_sig), "\n")
print(gsea_hallmark_sig[, c("pathway", "NES", "padj")])



# Clean pathway names and plot
gsea_hallmark_sig %>%
    as.data.frame() %>%
    mutate(
        pathway_clean = gsub("HALLMARK_", "", pathway),
        direction = ifelse(NES > 0, "Up in tumour", "Down in tumour")
    ) %>%
    ggplot(aes(x    = reorder(pathway_clean, NES),
               y    = NES,
               fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Up in tumour"   = "red",
                                  "Down in tumour" = "blue")) +
    geom_hline(yintercept = 0, linetype = "solid") +
    labs(x        = "",
         y        = "Normalised Enrichment Score (NES)",
         fill     = "",
         title    = "Hallmark pathways enriched in COAD tumour vs normal",
         subtitle = "TCGA-COAD | GSEA") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          legend.position = "bottom")

ggsave("results/gsea/GSEA_hallmark.pdf", width = 10, height = 8)