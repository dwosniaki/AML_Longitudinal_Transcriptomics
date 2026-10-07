#!/bin/bash

# ============================================================
# PROJECT: RNA-Seq AML
# SCRIPT 00 — PROCESSING AND QC
# ============================================================

# ============================================================
# Bulk RNA-seq Processing
# ============================================================
#
# Description:
# This script describes the preprocessing and alignment of
# paired-end bulk RNA-seq data, including the generation of
# STAR chimeric junction files for downstream circRNA analysis.
#
# All commands should be run in a Bash terminal.
#
# Reference genome and annotation:
# Ensembl GRCh38, release 115
#
# Main tools:
#   - STAR
#   - FastQC
#   - MultiQC
#   - Cutadapt
#
# ============================================================
# Pipeline overview
# ============================================================
#
# Phase 0: Reference genome preparation
#   Step 1 - Download and extract FASTA and GTF files
#   Step 2 - Generate STAR genome index
#
# Phase 1: Quality control and trimming
#   Step 1 - Run FastQC on raw reads
#   Step 2 - Run MultiQC on raw reads
#   Step 3 - Perform adapter and quality trimming
#   Step 4 - Run FastQC after trimming
#   Step 5 - Run MultiQC after trimming
#
# Phase 2: Alignment
#   Step 1 - Align trimmed reads with STAR
#   Step 2 - Generate gene-level counts and chimeric junction files
#
# ============================================================

# ============================================================
# Requirements
# ============================================================
#
# The following software should be installed before running
# this pipeline:
#
#   STAR 2.7.11a
#   FastQC
#   MultiQC
#   Cutadapt
#
# FASTA and GTF files should be downloaded from Ensembl.
#
# ============================================================

# ============================================================
# Phase 0: Reference genome preparation
# ============================================================

# ------------------------------------------------------------
# Phase 0 - Step 1
# Download and extract FASTA and GTF files
# ------------------------------------------------------------
#
# The reference FASTA contains the genome sequence.
#
# The GTF file contains genomic annotations for genes, exons,
# transcripts, and other genomic features.
#
# This pipeline uses:
#   - Homo sapiens GRCh38
#   - Ensembl release 115

BASE_DIR="YOUR_PROJECT_PATH"

mkdir -p "$BASE_DIR/Ensembl"
cd "$BASE_DIR/Ensembl"

wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz

# ------------------------------------------------------------
# Phase 0 - Step 2
# Generate the STAR genome index
# ------------------------------------------------------------
#
# The STAR genome index is generated from the reference FASTA
# and GTF annotation files and is subsequently used for
# RNA-seq read alignment.
#
# Change --runThreadN according to the computational resources
# available on your system.

BASE_DIR="YOUR_PROJECT_PATH"

"$BASE_DIR/Software/STAR-2.7.11a/bin/Linux_x86_64_static/STAR" \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir "$BASE_DIR/Ensembl/star_index" \
    --genomeFastaFiles "$BASE_DIR/Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --sjdbGTFfile "$BASE_DIR/Ensembl/Homo_sapiens.GRCh38.115.gtf" \
    --sjdbOverhang 100

# ============================================================
# Phase 1: Quality control and trimming
# ============================================================

# ------------------------------------------------------------
# Phase 1 - Step 1
# Run FastQC on raw sequencing reads
# ------------------------------------------------------------
#
# FastQC evaluates the quality of individual FASTQ files.
#
# The example below runs FastQC for one sample.
# Change SAMPLE to the sample identifier you want to process.
#
# To process all samples, adapt FASTQ_DIR accordingly.

BASE_DIR="YOUR_PROJECT_PATH"
SAMPLE="YOUR_SAMPLE_ID"

FASTQ_DIR="$BASE_DIR/Raw_fastq/$SAMPLE"
RESULTS_DIR="$BASE_DIR/Results/QC"
THREADS=8

mkdir -p "$RESULTS_DIR/$SAMPLE/FastQC"

echo "Running FastQC for $SAMPLE"

fastqc \
    "$FASTQ_DIR"/*.fastq.gz \
    -o "$RESULTS_DIR/$SAMPLE/FastQC" \
    -t "$THREADS"

# ------------------------------------------------------------
# Phase 1 - Step 2
# Run MultiQC on raw sequencing data
# ------------------------------------------------------------
#
# MultiQC aggregates FastQC results from multiple samples
# into a single report.

BASE_DIR="YOUR_PROJECT_PATH"

RESULTS_DIR="$BASE_DIR/Results/QC"

multiqc \
    "$RESULTS_DIR" \
    -o "$RESULTS_DIR/MultiQC"

# ------------------------------------------------------------
# Phase 1 - Step 3
# Perform adapter and quality trimming
# ------------------------------------------------------------
#
# Trimming removes adapter sequences and low-quality bases
# before alignment.
#
# Parameters used in this analysis:
#   - Quality cutoff: 20
#   - Minimum read length: 20 bp
#
# Change the adapter sequences if your library preparation
# uses different adapters.

BASE_DIR="YOUR_PROJECT_PATH"

RAW_DIR="$BASE_DIR/Raw_fastq"
TRIM_DIR="$BASE_DIR/Results/Trimmed"

THREADS=1

ADAPTER_R1="ACTGTCTCTTATACACATCT"
ADAPTER_R2="ACTGTCTCTTATACACATCT"

mkdir -p "$TRIM_DIR"

for sample_dir in "$RAW_DIR"/*; do

    sample=$(basename "$sample_dir")

    R1="$sample_dir/${sample}_R1.fastq.gz"
    R2="$sample_dir/${sample}_R2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Skipping $sample (FASTQ missing)"
        continue
    fi

    echo "Trimming $sample"

    cutadapt \
        -a "$ADAPTER_R1" \
        -a "$ADAPTER_R1" \
        -A "$ADAPTER_R2" \
        -A "$ADAPTER_R2" \
        -q 20 \
        --minimum-length 20 \
        -j "$THREADS" \
        -o "$TRIM_DIR/${sample}_R1.trimmed.fastq.gz" \
        -p "$TRIM_DIR/${sample}_R2.trimmed.fastq.gz" \
        "$R1" "$R2"

    if [ $? -ne 0 ]; then
        echo "ERROR trimming $sample. Aborting."
        exit 1
    fi

done

# ------------------------------------------------------------
# Phase 1 - Step 4
# Run FastQC after trimming
# ------------------------------------------------------------
#
# FastQC is repeated after trimming to evaluate the quality
# of the processed reads.

BASE_DIR="YOUR_PROJECT_PATH"

FASTTRIM_DIR="$BASE_DIR/Results/Trimmed"
RESULTS_DIR="$BASE_DIR/Results/QC_after_trimming"

THREADS=4

mkdir -p "$RESULTS_DIR/FastQC"

fastqc \
    "$FASTTRIM_DIR"/*.fastq.gz \
    -o "$RESULTS_DIR/FastQC" \
    -t "$THREADS"

# ------------------------------------------------------------
# Phase 1 - Step 5
# Run MultiQC after trimming
# ------------------------------------------------------------
#
# MultiQC aggregates the post-trimming FastQC results into
# a single report.

BASE_DIR="YOUR_PROJECT_PATH"

RESULTS_DIR="$BASE_DIR/Results/QC_after_trimming"

multiqc \
    "$RESULTS_DIR" \
    -o "$RESULTS_DIR/MultiQC"

# ============================================================
# Phase 2: Alignment
# ============================================================

# ------------------------------------------------------------
# Phase 2 - Step 1
# STAR alignment
# ------------------------------------------------------------
#
# STAR performs splice-aware alignment of trimmed paired-end
# RNA-seq reads against the reference genome.
#
# The following options were used to:
#   - Perform two-pass alignment
#   - Generate coordinate-sorted BAM files
#   - Generate gene-level read counts
#   - Detect chimeric junctions for downstream circRNA analysis
#
# The STAR parameters below should not be changed unless the
# analysis design is intentionally modified.

set -euo pipefail
shopt -s nullglob

BASE_DIR="YOUR_PROJECT_PATH"

TRIM_DIR="$BASE_DIR/Results/Trimmed"
ALIGN_DIR="$BASE_DIR/Results/Aligned_trimmed"
GENOME_DIR="$BASE_DIR/Ensembl/star_index"
STAR_BIN="$BASE_DIR/Software/STAR-2.7.11a/bin/Linux_x86_64_static/STAR"

THREADS=16


# Check whether required files and directories exist.

[[ -x "$STAR_BIN" ]] || {
    echo "[ERROR] STAR executable not found."
    exit 1
}

[[ -d "$TRIM_DIR" ]] || {
    echo "[ERROR] TRIM_DIR not found."
    exit 1
}

[[ -d "$GENOME_DIR" ]] || {
    echo "[ERROR] GENOME_DIR not found."
    exit 1
}


mkdir -p "$ALIGN_DIR"


# Align all paired-end trimmed samples.

for R1 in "$TRIM_DIR"/*_R1.trimmed.fastq.gz; do

    sample=$(basename "$R1" _R1.trimmed.fastq.gz)

    R2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "[WARN] R2 not found for $sample — skipping"
        continue
    fi

    echo "[INFO] $(date '+%F %T') - Aligning $sample"

    SAMPLE_DIR="$ALIGN_DIR/$sample"

    mkdir -p "$SAMPLE_DIR"

    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$R1" "$R2" \
        --outReadsUnmapped None \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix "$SAMPLE_DIR/${sample}_" \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 10 \
        --alignSJDBoverhangMin 5 \
        --alignIntronMin 20 \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:GRPundef \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --outSAMattributes NH HI AS nM NM MD

done

# ============================================================
# END OF SCRIPT 00
# ============================================================