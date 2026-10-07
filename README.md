# AML_Longitudinal_Transcriptomics
Reproducible analysis of longitudinal RNA-seq data from paired AML samples collected at diagnosis (D0) and after treatment (D14).

# Longitudinal Transcriptomic Analysis of Therapeutic Response in Acute Myeloid Leukemia

## Overview

This repository contains the R scripts used for the analysis of longitudinal bulk RNA-seq data from patients with acute myeloid leukemia (AML).

The study evaluates changes in global gene expression between diagnosis (D0) and day 14 (D14) after the initiation of treatment using paired samples from the same patients.

The main objective of this analysis is to characterize transcriptional changes associated with the early treatment response in AML.

---

## Study Design

The dataset consists of paired RNA-seq samples collected at two time points:

- **D0:** diagnosis, before treatment
- **D14:** 14 days after the initiation of the first treatment cycle

The global longitudinal analysis includes:

- **16 patients**
- **32 RNA-seq samples**
- **16 paired D0/D14 comparisons**

Samples corresponding to patients treated with ATRA were excluded from the global analysis.

Differential expression was evaluated using a paired DESeq2 design:

```text
~ Patient + Day
```

The main contrast was:

```text
D14 vs D0
```

---

## Analysis Workflow

The analysis is organized into sequential scripts:

```text
Raw RNA-seq data
        │
        ▼
Gene-level count matrix
        │
        ▼
Master expression matrix
        │
        ├── Global PCA
        │
        ├── Differential expression analysis
        │
        ├── Top DEGs heatmap
        │
        ├── Volcano plot
        │
        └── GSEA — Hallmark pathways
```

---

## Repository Structure

```text
.
├── Papers/
│   └── Session_1/
│       ├── scripts/
│       │   ├── Script_01_Initial_Data_Preparation.R
│       │   ├── Script_02_Global_PCA.R
│       │   ├── Script_03_Global_DESeq2.R
│       │   ├── Script_04_Top200_DEGs_Heatmap.R
│       │   ├── Script_05_Volcano_Plot.R
│       │   └── Script_06_GSEA_Hallmark.R
│       │
│       ├── data/
│       ├── objects/
│       └── figures/
│
├── .gitignore
└── README.md
```

The `data`, `objects`, and `figures` directories are excluded from version control.

---

## Scripts

### Script 01 — Initial Data Preparation

Prepares the master expression matrix for downstream analysis.

Main steps include:

- Loading the master expression matrix
- Identifying annotation and expression columns
- Calculating CPM
- Calculating library sizes
- Generating summary information
- Saving processed expression objects

---

### Script 02 — Global PCA

Performs principal component analysis using log2-transformed CPM values.

The analysis:

- Evaluates global expression structure
- Includes paired D0 and D14 samples
- Connects samples from the same patient
- Uses the 32 samples included after exclusion of ATRA-treated patients

---

### Script 03 — Global DESeq2

Performs differential expression analysis between D14 and D0.

Design:

```text
~ Patient + Day
```

Contrast:

```text
D14 vs D0
```

Significant DEGs are defined using:

```text
padj < 0.05
|log2FC| >= 1
```

The complete DESeq2 results are retained for downstream analyses, including GSEA.

---

### Script 04 — Top 200 DEGs Heatmap

Generates a heatmap of the 200 most significant DEGs from the global D14 vs D0 comparison.

Selection criteria:

```text
padj < 0.05
|log2FC| >= 1
```

The 200 most significant genes are selected according to adjusted p-value.

For visualization:

- Variance stabilizing transformation (VST) is applied
- Gene-wise Z-scores are calculated
- Samples are ordered by Day and Patient
- Rows are hierarchically clustered
- Columns are not clustered

---

### Script 05 — Volcano Plot

Generates a volcano plot representing the global differential expression analysis.

The plot displays:

- **X-axis:** log2 Fold Change
- **Y-axis:** -log10(adjusted p-value)

Genes are classified as:

- Upregulated
- Downregulated
- Not significant

using:

```text
padj < 0.05
|log2FC| >= 1
```

---

### Script 06 — GSEA Hallmark

Performs Gene Set Enrichment Analysis using the MSigDB Hallmark collection.

The ranking metric is the DESeq2 Wald statistic:

```text
stat
```

The complete ranked gene list is used without applying differential expression cutoffs before GSEA.

Interpretation of the ranking:

```text
stat > 0  → enrichment toward D14
stat < 0  → enrichment toward D0
```

Pathways are considered significant using:

```text
adjusted p-value < 0.05
```

The analysis generates a Hallmark pathway dot plot based on:

- Normalized Enrichment Score (NES)
- Adjusted p-value

---

## Software and Packages

The analysis was performed using R and the following main packages:

- **DESeq2**
- **edgeR**
- **ComplexHeatmap**
- **circlize**
- **clusterProfiler**
- **msigdbr**
- **ggplot2**
- **tidyverse**
- **readxl**
- **writexl**

Specific package versions and the R session information can be obtained from the corresponding analysis environment.

---

## Data Availability

The original sequencing data, processed expression matrices, metadata containing patient-level information, intermediate analysis objects, and generated results are **not included in this repository**.

The repository contains the analysis scripts required to document and reproduce the computational workflow.

Access to the underlying data is subject to the applicable ethical, institutional, and data-sharing restrictions.

---

## Reproducibility

Before running the scripts:

1. Clone or download this repository.
2. Place the required input data in the corresponding local project directories.
3. Update the project path in each script:

```r
base_dir <- "YOUR_PROJECT_PATH/RNAseq_Analysis"
```

4. Install the required R packages.
5. Run the scripts sequentially.

The scripts are numbered according to the intended analysis workflow.

---

## Analysis Outputs

The workflow generates:

- Processed expression matrices
- Library size summaries
- PCA coordinates and variance information
- DESeq2 results
- Differentially expressed gene lists
- Top 200 DEG heatmap
- Volcano plot
- Hallmark GSEA results
- Hallmark GSEA dot plot
