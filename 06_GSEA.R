################################################################################
# SCRIPT 06 — GSEA HALLMARK
# Paper 1 — Global transcriptomic landscape
################################################################################

# This script performs Gene Set Enrichment Analysis (GSEA) using the complete
# DESeq2 results from the global longitudinal comparison between diagnosis
# (D0) and day 14 (D14) after treatment.
#
# The analysis uses the DESeq2 Wald statistic (stat) as the ranking metric.
#
# Comparison:
#
#   D14 vs D0
#
# DESeq2 design:
#
#   ~ Patient + Day
#
# Contrast:
#
#   Day D14 vs D0
#
# Gene set collection:
#
#   MSigDB Hallmark
#
# Ranking:
#
#   DESeq2 Wald statistic (stat)
#
# Positive stat values indicate enrichment toward D14, whereas negative
# stat values indicate enrichment toward D0.
#
# No padj filtering is applied before GSEA. All genes with an available
# DESeq2 Wald statistic and valid gene symbol are included in the ranking.
#
# GSEA is performed using the complete ranked gene list and does not apply
# differential expression thresholds before pathway enrichment analysis.
################################################################################


################################################################################
# 1. LOAD PACKAGES
################################################################################

library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(ggplot2)


################################################################################
# 2. DEFINE PATHS
################################################################################

# Change this path to your local project directory
base_dir <- "YOUR_PROJECT_PATH/RNAseq_Analysis"

papers_dir <- file.path(
  base_dir,
  "Papers"
)

paper1_dir <- file.path(
  papers_dir,
  "Paper_1"
)

data_dir <- file.path(
  paper1_dir,
  "data"
)

figures_dir <- file.path(
  paper1_dir,
  "figures"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# 3. DEFINE FILE NAMES
################################################################################

# Change only this name if the final figure name needs to be changed
figure_name <- "NAME_OF_THE_FIGURE"

de_results_file <- file.path(
  data_dir,
  "Paper1_master_DESeq2_D14_vs_D0_all.csv"
)


################################################################################
# 4. LOAD COMPLETE DESeq2 RESULTS
################################################################################

res_df <- read.csv(
  de_results_file,
  stringsAsFactors = FALSE
)


################################################################################
# 5. CHECK DESeq2 RESULTS
################################################################################

cat("\n")
cat("============================================\n")
cat("DESeq2 RESULT CHECK\n")
cat("============================================\n")

cat(
  "Genes in result:",
  nrow(res_df),
  "\n"
)

cat(
  "Genes with gene_name:",
  sum(
    !is.na(res_df$gene_name) &
      res_df$gene_name != ""
  ),
  "\n"
)

cat(
  "Genes with stat:",
  sum(
    !is.na(res_df$stat)
  ),
  "\n"
)

cat("============================================\n")


################################################################################
# 6. PREPARE GSEA RANKING
################################################################################

# The DESeq2 Wald statistic is used as the ranking metric.
#
# stat > 0:
#   enrichment toward D14
#
# stat < 0:
#   enrichment toward D0
#
# No padj filtering is applied.
# All genes with an available statistic are included in the ranking.

gsea_df <- res_df %>%
  filter(
    !is.na(stat),
    !is.na(gene_name),
    gene_name != ""
  ) %>%
  select(
    gene_id,
    gene_name,
    stat
  ) %>%
  distinct(
    gene_name,
    .keep_all = TRUE
  ) %>%
  arrange(
    desc(stat)
  )


################################################################################
# 7. CREATE RANKED GENE LIST
################################################################################

gene_list <- gsea_df$stat

names(gene_list) <- gsea_df$gene_name

gene_list <- gene_list[
  !duplicated(names(gene_list))
]

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)


################################################################################
# 8. CHECK RANKING
################################################################################

cat("\n")
cat("============================================\n")
cat("GSEA RANKING\n")
cat("============================================\n")

cat(
  "Genes in ranking:",
  length(gene_list),
  "\n"
)

cat(
  "Maximum stat:",
  max(gene_list),
  "\n"
)

cat(
  "Minimum stat:",
  min(gene_list),
  "\n"
)

cat("\nFirst genes:\n")

print(
  head(
    names(gene_list),
    10
  )
)

cat("\n============================================\n")


################################################################################
# 9. LOAD MSigDB HALLMARK
################################################################################

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)


################################################################################
# 10. PREPARE TERM2GENE
################################################################################

# The current msigdbr version used in the analysis provides gene symbols
# through the gene_symbol column.

hallmark_t2g <- hallmark %>%
  select(
    gs_name,
    gene_symbol
  ) %>%
  filter(
    !is.na(gs_name),
    gs_name != "",
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  distinct()


################################################################################
# 11. STANDARDIZE GENE SYMBOLS
################################################################################

# Ensure that the gene identifiers used by DESeq2 and MSigDB
# are represented in the same format.

names(gene_list) <- toupper(
  names(gene_list)
)

hallmark_t2g <- hallmark_t2g %>%
  mutate(
    gene_symbol = toupper(gene_symbol)
  ) %>%
  distinct()


################################################################################
# 12. CHECK GENE OVERLAP
################################################################################

overlap_genes <- intersect(
  names(gene_list),
  unique(
    hallmark_t2g$gene_symbol
  )
)

cat("\n")
cat("============================================\n")
cat("GENE ID COMPATIBILITY\n")
cat("============================================\n")

cat(
  "Genes in ranking:",
  length(gene_list),
  "\n"
)

cat(
  "Hallmark genes:",
  length(
    unique(
      hallmark_t2g$gene_symbol
    )
  ),
  "\n"
)

cat(
  "Genes in common:",
  length(overlap_genes),
  "\n"
)

cat(
  "Percentage of ranking with correspondence:",
  round(
    100 *
      length(overlap_genes) /
      length(gene_list),
    2
  ),
  "%\n"
)

cat("\nExamples of overlapping genes:\n")

print(
  head(
    overlap_genes,
    20
  )
)

cat("============================================\n")


################################################################################
# 13. SAFETY CHECK
################################################################################

if (length(overlap_genes) == 0) {
  
  stop(
    paste(
      "\nERROR: no genes from the ranking were found",
      "in the Hallmark gene sets.",
      "\nCheck the identifiers used in gene_name."
    )
  )
}


################################################################################
# 14. RUN GSEA
################################################################################

gsea_result <- GSEA(
  
  geneList = gene_list,
  
  TERM2GENE = hallmark_t2g,
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 1,
  
  pAdjustMethod = "BH",
  
  eps = 0,
  
  verbose = FALSE
)


################################################################################
# 15. CONVERT RESULTS TO DATA FRAME
################################################################################

gsea_df_result <- as.data.frame(
  gsea_result
)


################################################################################
# 16. CHECK GSEA RESULTS
################################################################################

cat("\n")
cat("============================================\n")
cat("GSEA RESULTS\n")
cat("============================================\n")

cat(
  "Pathways tested:",
  nrow(gsea_df_result),
  "\n"
)

cat(
  "Pathways FDR < 0.05:",
  sum(
    gsea_df_result$p.adjust < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Enriched in D14 (NES > 0):",
  sum(
    gsea_df_result$p.adjust < 0.05 &
      gsea_df_result$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Enriched in D0 (NES < 0):",
  sum(
    gsea_df_result$p.adjust < 0.05 &
      gsea_df_result$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)

cat("============================================\n")


################################################################################
# 17. CREATE PATHWAY NAME
################################################################################

gsea_df_result$Pathway <- gsub(
  "^HALLMARK_",
  "",
  gsea_df_result$Description
)


################################################################################
# 18. SAVE COMPLETE RESULTS
################################################################################

write.csv(
  gsea_df_result,
  file.path(
    data_dir,
    "Paper1_master_GSEA_Hallmark_all.csv"
  ),
  row.names = FALSE
)


################################################################################
# 19. SIGNIFICANT PATHWAYS
################################################################################

gsea_sig <- gsea_df_result %>%
  filter(
    p.adjust < 0.05
  ) %>%
  arrange(
    p.adjust
  )


################################################################################
# 20. PREPARE DOTPLOT
################################################################################

# Select the 10 most significant pathways in each direction.

gsea_plot_df <- gsea_df_result %>%
  filter(
    p.adjust < 0.05
  ) %>%
  mutate(
    
    Direction = ifelse(
      NES > 0,
      "D14",
      "D0"
    ),
    
    neg_log10_FDR = -log10(
      p.adjust
    )
    
  ) %>%
  group_by(
    Direction
  ) %>%
  slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()


################################################################################
# 21. ORDER PATHWAYS
################################################################################

gsea_plot_df <- gsea_plot_df %>%
  arrange(
    NES
  )

gsea_plot_df$Pathway <- factor(
  gsea_plot_df$Pathway,
  levels = gsea_plot_df$Pathway
)


################################################################################
# 22. CREATE GSEA DOTPLOT
################################################################################

gsea_plot <- ggplot(
  gsea_plot_df,
  aes(
    x = NES,
    y = Pathway
  )
) +
  
  geom_point(
    aes(
      size = neg_log10_FDR,
      color = NES
    ),
    alpha = 0.9
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey40"
  ) +
  
  scale_color_gradient2(
    low = "#164A8A",
    mid = "white",
    high = "#8B1E3F",
    midpoint = 0,
    name = "NES"
  ) +
  
  scale_size_continuous(
    name = "-log10 (adjusted p-value)",
    range = c(2, 7)
  ) +
  
  labs(
    title = "GSEA - Hallmark pathways",
    subtitle = "D14 vs D0",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  
  theme_bw() +
  
  theme(
    
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 9,
      color = "grey40",
      hjust = 0
    ),
    
    axis.title.x = element_text(
      size = 8
    ),
    
    axis.text.x = element_text(
      size = 8
    ),
    
    axis.text.y = element_text(
      size = 7
    ),
    
    panel.grid.major = element_line(
      color = "white"
    ),
    
    panel.grid.minor = element_line(
      color = "white"
    ),
    
    panel.background = element_rect(
      fill = "grey95",
      color = NA
    ),
    
    panel.border = element_rect(
      color = "white",
      fill = NA,
      linewidth = 0.8
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    legend.position = "right"
  )


################################################################################
# 23. DISPLAY GSEA DOTPLOT
################################################################################

print(gsea_plot)


################################################################################
# 24. SAVE PNG
################################################################################

ggsave(
  filename = file.path(
    figures_dir,
    paste0(
      figure_name,
      ".png"
    )
  ),
  plot = gsea_plot,
  width = 7,
  height = 7,
  dpi = 300
)


################################################################################
# 25. SAVE PDF
################################################################################

ggsave(
  filename = file.path(
    figures_dir,
    paste0(
      figure_name,
      ".pdf"
    )
  ),
  plot = gsea_plot,
  width = 7,
  height = 7
)


################################################################################
# 26. SAVE SIGNIFICANT RESULTS
################################################################################

write.csv(
  gsea_sig,
  file.path(
    data_dir,
    "Paper1_master_GSEA_Hallmark_significant.csv"
  ),
  row.names = FALSE
)


################################################################################
# 27. FINAL SUMMARY
################################################################################

cat("\n")
cat("============================================\n")
cat("SCRIPT 06 — GSEA HALLMARK COMPLETED\n")
cat("============================================\n")

cat(
  "Comparison: D14 vs D0\n"
)

cat(
  "Ranking: DESeq2 Wald statistic (stat)\n"
)

cat(
  "Gene set: MSigDB Hallmark\n"
)

cat(
  "Pathways tested:",
  nrow(gsea_df_result),
  "\n"
)

cat(
  "Significant pathways (FDR < 0.05):",
  nrow(gsea_sig),
  "\n\n"
)

cat(
  "Complete results:\n",
  file.path(
    data_dir,
    "Paper1_master_GSEA_Hallmark_all.csv"
  ),
  "\n\n"
)

cat(
  "Significant results:\n",
  file.path(
    data_dir,
    "Paper1_master_GSEA_Hallmark_significant.csv"
  ),
  "\n\n"
)

cat(
  "PNG:\n",
  file.path(
    figures_dir,
    paste0(
      figure_name,
      ".png"
    )
  ),
  "\n\n"
)

cat(
  "PDF:\n",
  file.path(
    figures_dir,
    paste0(
      figure_name,
      ".pdf"
    )
  ),
  "\n"
)

# ============================================================
# END OF SCRIPT 06
# ============================================================