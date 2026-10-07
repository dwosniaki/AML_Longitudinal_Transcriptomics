# ============================================================
# SCRIPT 04 — TOP 200 DEGs HEATMAP
# Paper 1 — Global transcriptomic landscape
# ============================================================
#
# This script generates a heatmap of the 200 most significant
# differentially expressed genes identified in the global
# D0 versus D14 comparison.

# ============================================================
# 1. PACKAGES
# ============================================================

library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(dplyr)


# ============================================================
# 2. DIRECTORIES
# ============================================================

# Change this path to your local project directory
base_dir <- "YOUR_PROJECT_PATH/RNAseq_Analysis"

papers_dir  <- file.path(base_dir, "Papers")
paper1_dir  <- file.path(papers_dir, "Paper_1")

data_dir    <- file.path(paper1_dir, "data")
objects_dir <- file.path(paper1_dir, "objects")
figures_dir <- file.path(paper1_dir, "figures")

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. FILE NAMES
# ============================================================

# Change only this name if the final figure name needs to be changed
figure_name <- "NAME_OF_THE_FIGURE"

degs_file <- file.path(
  data_dir,
  "Paper1_master_DESeq2_D14_vs_D0_all.csv"
)

dds_file <- file.path(
  objects_dir,
  "Paper1_master_DESeq2_D14_vs_D0.rds"
)

metadata_file <- file.path(
  papers_dir,
  "Metadata.xlsx"
)


# ============================================================
# 4. LOAD DESEQ2 RESULTS
# ============================================================

dds <- readRDS(dds_file)

# Complete DESeq2 results
res_df <- read.csv(
  degs_file,
  stringsAsFactors = FALSE
)


# ============================================================
# 5. IDENTIFY THE TOP 200 DEGs
# ============================================================

# Selection criteria:
# padj < 0.05
# |log2FC| >= 1

DEGs <- res_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  ) %>%
  arrange(padj)


# Select the 200 most significant DEGs
top200 <- DEGs %>%
  slice_head(n = 200)


cat("\n============================================\n")
cat("TOP 200 DEGs\n")
cat("============================================\n")

cat(
  "Total number of DEGs:",
  nrow(DEGs),
  "\n"
)

cat(
  "Number selected for the heatmap:",
  nrow(top200),
  "\n"
)


# ============================================================
# 6. SAVE TOP 200 DEG LIST
# ============================================================

write.csv(
  top200,
  file.path(
    data_dir,
    "Paper1_master_Top200_DEGs_heatmap.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. VST TRANSFORMATION
# ============================================================

# VST is used only for visualization.
# DESeq2 was previously performed using raw counts.

vsd <- vst(
  dds,
  blind = FALSE
)

vsd_matrix <- assay(vsd)


# ============================================================
# 8. SELECT THE TOP 200 GENES
# ============================================================

# The DESeq2 results contain gene_id as the identifier.
# Only genes present in the VST matrix are retained.

gene_ids <- top200$gene_id

gene_ids <- gene_ids[
  gene_ids %in% rownames(vsd_matrix)
]

heatmap_matrix <- vsd_matrix[
  gene_ids,
  ,
  drop = FALSE
]


cat(
  "\nGenes found in the VST:",
  nrow(heatmap_matrix),
  "\n"
)


# ============================================================
# 9. Z-SCORE PER GENE
# ============================================================

# Each gene is standardized across samples.

heatmap_z <- t(
  scale(
    t(heatmap_matrix)
  )
)

# Remove genes that produced NA values
heatmap_z <- heatmap_z[
  complete.cases(heatmap_z),
  ,
  drop = FALSE
]


# ============================================================
# 10. METADATA
# ============================================================

metadata <- read_excel(
  metadata_file
)

# Keep only samples present in the heatmap
metadata <- metadata[
  match(
    colnames(heatmap_z),
    metadata$Sample_ID
  ),
  ,
  drop = FALSE
]

# Check sample correspondence
stopifnot(
  all(
    metadata$Sample_ID ==
      colnames(heatmap_z)
  )
)

metadata$Day <- factor(
  metadata$Day,
  levels = c("D0", "D14")
)

metadata$Patient <- factor(
  metadata$Patient_number
)


# ============================================================
# 11. ORDER SAMPLES
# D0 first -> D14 second
# Patients remain in the same order within both groups
# ============================================================

metadata$Patient_number <- as.numeric(
  metadata$Patient_number
)

metadata$Day <- factor(
  metadata$Day,
  levels = c("D0", "D14")
)

metadata <- metadata %>%
  arrange(
    Day,
    Patient_number
  )

sample_order <- metadata$Sample_ID

heatmap_z <- heatmap_z[
  ,
  sample_order,
  drop = FALSE
]


# ============================================================
# 12. COLUMN LABELS
# ============================================================

column_labels <- paste0(
  "Patient ",
  metadata$Patient_number
)

colnames(heatmap_z) <- column_labels


# ============================================================
# 13. ANNOTATIONS: DAY + PATIENT
# ============================================================

# Colors for the 16 patients
patient_levels <- sort(
  unique(
    metadata$Patient_number
  )
)

patient_colors <- colorRampPalette(
  c(
    "#8B1E3F",
    "#F4B6C2"
  )
)(
  length(patient_levels)
)

names(patient_colors) <- paste0(
  "Patient ",
  patient_levels
)


# Patient categorical variable
metadata$Patient_label <- factor(
  paste0(
    "Patient ",
    metadata$Patient_number
  ),
  levels = paste0(
    "Patient ",
    sort(
      unique(
        metadata$Patient_number
      )
    )
  )
)


# Top annotation
top_annotation <- HeatmapAnnotation(
  
  Day = metadata$Day,
  
  Patient = metadata$Patient_label,
  
  col = list(
    
    # First annotation track
    Day = c(
      D0 = "#164A8A",
      D14 = "#63CFC3"
    ),
    
    # Second annotation track
    Patient = patient_colors
  ),
  
  annotation_name_side = "left",
  
  annotation_name_gp = gpar(
    fontsize = 9
  ),
  
  simple_anno_size = unit(
    5,
    "mm"
  )
)


# ============================================================
# 14. EXPRESSION COLORS
# ============================================================

col_fun <- colorRamp2(
  c(
    -2,
    0,
    2
  ),
  c(
    "#2166AC",
    "white",
    "#8B1E3F"
  )
)


# ============================================================
# 15. GENE NAMES
# ============================================================

# Retrieve gene_name corresponding to each gene_id

gene_annotation <- res_df %>%
  select(
    gene_id,
    gene_name
  ) %>%
  distinct()

gene_names <- gene_annotation$gene_name[
  match(
    rownames(heatmap_z),
    gene_annotation$gene_id
  )
]

# If gene_name is unavailable, use gene_id
gene_names[
  is.na(gene_names) |
    gene_names == ""
] <- rownames(
  heatmap_z
)[
  is.na(gene_names) |
    gene_names == ""
]

rownames(heatmap_z) <- gene_names


# ============================================================
# 16. HEATMAP
# ============================================================

ht <- Heatmap(
  
  heatmap_z,
  
  name = "Z-score",
  
  col = col_fun,
  
  top_annotation = top_annotation,
  
  cluster_rows = TRUE,
  
  cluster_columns = FALSE,
  
  row_dend_gp = gpar(
    lwd = 0.5
  ),
  
  show_row_names = TRUE,
  
  row_names_gp = gpar(
    fontsize = 5,
    fontface = "italic"
  ),
  
  row_names_side = "right",
  
  # Do not show sample names
  show_column_names = FALSE,
  
  # Title
  column_title = "Global gene expression",
  
  column_title_side = "top",
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  border = TRUE,
  
  border_gp = gpar(
    col = "white",
    lwd = 0.8
  ),
  
  use_raster = TRUE
)


# ============================================================
# 17. EXPORT PDF
# ============================================================

pdf(
  file.path(
    figures_dir,
    paste0(
      figure_name,
      ".pdf"
    )
  ),
  width = 12,
  height = 10
)

draw(
  ht,
  
  heatmap_legend_side = "right",
  
  annotation_legend_side = "right",
  
  merge_legends = TRUE
)

dev.off()


# ============================================================
# 18. EXPORT PNG
# ============================================================

png(
  file.path(
    figures_dir,
    paste0(
      figure_name,
      ".png"
    )
  ),
  width = 3600,
  height = 3000,
  res = 300
)

draw(
  ht,
  
  heatmap_legend_side = "right",
  
  annotation_legend_side = "right",
  
  merge_legends = TRUE
)

dev.off()


# ============================================================
# 19. FINAL SUMMARY
# ============================================================

cat("\n============================================\n")
cat("HEATMAP COMPLETED\n")
cat("============================================\n")

cat(
  "Total DEGs:",
  nrow(DEGs),
  "\n"
)

cat(
  "Genes in heatmap:",
  nrow(heatmap_z),
  "\n"
)

cat(
  "Samples:",
  ncol(heatmap_z),
  "\n"
)

cat("\nGenerated files:\n")

cat(
  file.path(
    data_dir,
    "Paper1_master_Top200_DEGs_heatmap.csv"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    paste0(
      figure_name,
      ".pdf"
    )
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    paste0(
      figure_name,
      ".png"
    )
  ),
  "\n"
)

# ============================================================
# END OF SCRIPT 04
# ============================================================