# ============================================================
# PROJECT: RNA-Seq AML
# SCRIPT 01 — INITIAL DATA PREPARATION
# ============================================================
#
# Purpose:
#   - Check and create the required directory structure
#   - Load the master expression matrix
#   - Validate the input data
#   - Generate raw-count matrices
#   - Generate CPM-normalized matrices
#   - Save data and objects
#
# No differential expression analysis is performed here.
# No expression filtering is performed here.
# Original source data are not modified.
# ============================================================


# ============================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================

library(edgeR)


# ============================================================
# 2. DEFINE BASE DIRECTORY
# ============================================================

# Change this path to the location of your RNA-seq analysis project.
#
# Example:
# base_dir <- "/path/to/RNAseq_Analysis"

base_dir <- "YOUR_PROJECT_PATH/RNAseq_Analysis"


# ============================================================
# 3. DEFINE PAPER DIRECTORIES
# ============================================================

paper1_dir <- file.path(
  base_dir,
  "Papers",
  "Paper_1"
)

paper2_dir <- file.path(
  base_dir,
  "Papers",
  "Paper_2"
)


# ============================================================
# 4. CREATE / CHECK DIRECTORY STRUCTURE
# ============================================================

directories <- c(
  
  # -------------------------
  # Paper 1 — mRNA
  # -------------------------
  
  file.path(paper1_dir, "data"),
  file.path(paper1_dir, "scripts"),
  file.path(paper1_dir, "objects"),
  file.path(paper1_dir, "figures"),
  
  # -------------------------
  # Paper 2 — lncRNA
  # -------------------------
  
  file.path(paper2_dir, "data"),
  file.path(paper2_dir, "scripts"),
  file.path(paper2_dir, "objects"),
  file.path(paper2_dir, "figures")
  
)

for (dir in directories) {
  
  if (!dir.exists(dir)) {
    
    dir.create(
      dir,
      recursive = TRUE
    )
    
    message("Created: ", dir)
    
  } else {
    
    message("Already exists: ", dir)
    
  }
  
}


# ============================================================
# 5. DEFINE OUTPUT PATHS
# ============================================================

# -------------------------
# Paper 1
# -------------------------

paper1_data <- file.path(
  paper1_dir,
  "data"
)

paper1_objects <- file.path(
  paper1_dir,
  "objects"
)

paper1_figures <- file.path(
  paper1_dir,
  "figures"
)


# -------------------------
# Paper 2
# -------------------------

paper2_data <- file.path(
  paper2_dir,
  "data"
)

paper2_objects <- file.path(
  paper2_dir,
  "objects"
)

paper2_figures <- file.path(
  paper2_dir,
  "figures"
)


# ============================================================
# 6. DEFINE INPUT DATA
# ============================================================

# Change this path if the master expression matrix is stored
# in a different location within your project.

input_file <- file.path(
  base_dir,
  "data",
  "master.rds"
)


# ============================================================
# 7. CHECK INPUT FILE
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Input file not found:\n",
      input_file
    )
  )
  
}

message("Input file found:")
message(input_file)


# ============================================================
# 8. LOAD MASTER MATRIX
# ============================================================

master <- readRDS(
  input_file
)


# ============================================================
# 9. CHECK OBJECT STRUCTURE
# ============================================================

message("\nObject class:")
print(class(master))

message("\nObject dimensions:")
print(dim(master))

message("\nColumn names:")
print(colnames(master))


if (!is.data.frame(master) &&
    !is.matrix(master)) {
  
  stop(
    "The input object must be a data.frame or matrix."
  )
  
}


master <- as.data.frame(
  master,
  check.names = FALSE
)


# ============================================================
# 10. IDENTIFY ANNOTATION COLUMNS
# ============================================================

annotation_candidates <- c(
  "gene_id",
  "gene_name",
  "gene_biotype",
  "gene_source",
  "seqname",
  "start",
  "end",
  "strand"
)

annotation_columns <- intersect(
  annotation_candidates,
  colnames(master)
)

message("\nAnnotation columns detected:")
print(annotation_columns)


if (!"gene_biotype" %in% annotation_columns) {
  
  stop(
    "gene_biotype column was not detected."
  )
  
}


# ============================================================
# 11. IDENTIFY EXPRESSION / SAMPLE COLUMNS
# ============================================================

expression_columns <- setdiff(
  colnames(master),
  annotation_columns
)

message("\nExpression/sample columns detected:")
print(expression_columns)


if (length(expression_columns) == 0) {
  
  stop(
    "No expression/sample columns were detected."
  )
  
}


# ============================================================
# 12. EXTRACT ANNOTATION AND EXPRESSION DATA
# ============================================================

master_annotation <- master[
  ,
  annotation_columns,
  drop = FALSE
]

master_counts <- master[
  ,
  expression_columns,
  drop = FALSE
]


master_counts <- as.data.frame(
  lapply(
    master_counts,
    function(x) as.numeric(as.character(x))
  ),
  check.names = FALSE
)


# ============================================================
# 13. VALIDATE EXPRESSION DATA
# ============================================================

if (anyNA(master_counts)) {
  
  warning(
    "NA values were detected in the expression matrix."
  )
  
}


if (any(
  master_counts < 0,
  na.rm = TRUE
)) {
  
  stop(
    "Negative expression values were detected."
  )
  
}


# ============================================================
# 14. CHECK GENE BIOTYPE DISTRIBUTION
# ============================================================

message("\nGene biotype distribution:")

biotype_table <- sort(
  table(master$gene_biotype),
  decreasing = TRUE
)

print(biotype_table)


# ============================================================
# 15. CHECK DUPLICATE GENE IDs
# ============================================================

if ("gene_id" %in% colnames(master_annotation)) {
  
  if (anyDuplicated(master_annotation$gene_id) > 0) {
    
    warning(
      "Duplicated gene_id values detected in master."
    )
    
  }
  
}


# ============================================================
# 16. CREATE RAW MATRIX WITH ANNOTATION
# ============================================================

master_raw <- cbind(
  master_annotation,
  master_counts
)


# ============================================================
# 17. CALCULATE CPM
# ============================================================

master_cpm_values <- edgeR::cpm(
  as.matrix(master_counts),
  log = FALSE
)


# ============================================================
# 18. CREATE CPM MATRIX WITH ANNOTATION
# ============================================================

master_cpm <- cbind(
  master_annotation,
  as.data.frame(
    master_cpm_values,
    check.names = FALSE
  )
)


# ============================================================
# 19. CALCULATE LIBRARY SIZES
# ============================================================

master_library_size <- colSums(
  master_counts,
  na.rm = TRUE
)

library_sizes <- data.frame(
  sample = names(master_library_size),
  master_library_size = as.numeric(
    master_library_size
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 20. SAVE PAPER 1 — RAW DATA
# ============================================================

saveRDS(
  master_raw,
  file.path(
    paper1_objects,
    "Paper1_master_raw.rds"
  )
)

write.csv(
  master_raw,
  file.path(
    paper1_data,
    "Paper1_master_raw.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SAVE PAPER 1 — CPM DATA
# ============================================================

saveRDS(
  master_cpm,
  file.path(
    paper1_objects,
    "Paper1_master_CPM.rds"
  )
)

write.csv(
  master_cpm,
  file.path(
    paper1_data,
    "Paper1_master_CPM.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. SAVE LIBRARY SIZE INFORMATION
# ============================================================

write.csv(
  library_sizes,
  file.path(
    paper1_data,
    "library_sizes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. SAVE GENE BIOTYPE SUMMARY
# ============================================================

gene_biotype_summary <- data.frame(
  gene_biotype = names(biotype_table),
  n_genes = as.integer(biotype_table),
  stringsAsFactors = FALSE
)

write.csv(
  gene_biotype_summary,
  file.path(
    paper1_data,
    "gene_biotype_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. SAVE SESSION INFORMATION
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    paper1_objects,
    "sessionInfo.txt"
  )
)


# ============================================================
# 25. FINAL VALIDATION
# ============================================================

message("\n==============================================")
message("FINAL VALIDATION")
message("==============================================")


message(
  "Paper 1 — raw dimensions: ",
  paste(
    dim(master_raw),
    collapse = " x "
  )
)


message(
  "Paper 1 — CPM dimensions: ",
  paste(
    dim(master_cpm),
    collapse = " x "
  )
)


message(
  "Number of samples: ",
  ncol(master_counts)
)


message(
  "Number of genes: ",
  nrow(master_counts)
)


message("\nScript 01 completed successfully.")


# ============================================================
# END OF SCRIPT 01
# ============================================================