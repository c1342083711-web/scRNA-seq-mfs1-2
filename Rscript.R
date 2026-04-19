######Indexing, generation of count matrix using simpleaf(piscem)
#using genome.fa and genes.gtf annotation(cellranger, STARsolo also), but not transcript.fa(only require in Alevin-fry, kallisto, salmon)

## Transcriptome-based indexing (Kallisto, Salmon) is faster and uses less memory, but may miss reads that span splice junctions or novel transcripts.
## Genome-based indexing (STAR, Piscem, Cell Ranger) provides more accurate mapping, especially for spliced reads, but is heavier and slower.

setwd("N:/")
# Path to the reference transcriptome FASTA
ref_transcriptome="N:/refdata-gex-GRCh38-2024-A/transcripts.fa"

zcat /mnt/n/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq/MFS-1-1_S9_L003_R1_001.fastq.gz | head -n 20

# Output index directory
index_dir="N:/refdata-gex-GRCh38-2024-A/salmon_index"

REFERENCE_DIR="/mnt/n/refdata-gex-GRCh38-2024-A"
SAMPLE_ID="MFS-1"
READ1="/mnt/n/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq/MFS-1-1_S9_L003_R1_001.fastq.gz"
READ2="/mnt/n/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq/MFS-1-1_S9_L003_R2_001.fastq.gz"

kb-count --help

kb-count --h5 --workflow snrna --transcriptome ${REFERENCE_DIR}/transcripts.fa --genes ${REFERENCE_DIR}/genes.gtf --out ./kb_output_${SAMPLE_ID} --reads "$READ1" "$READ2"

install.packages(c("irlba", "scico"), Ncpus = 2)
library(Matrix)
library(irlba)
library(ggplot2) # Tidyverse is pre-installed, yay!
library(dplyr)
library(scico)
theme_set(theme_bw())

#install cran packages
install.packages(c("tidyverse", "Matrix", "patchwork",
                   "pheatmap", "RColorBrewer", "readxl"))


#install bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "scater",
                       "scran", "DropletUtils", "bluster",
                       "scDblFinder", "AUCell", "PCAtools"))

system.time({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c('multtest', "DropletUtils"), Ncpus = 2)
  install.packages(c("Seurat", "scico", "ggpointdensity"), Ncpus = 2)
})


install kb-python


Python --version
module purge
module load Conda/3
pip install kb-python

module purge
module load Conda/3

FASTQ_DIR="/mnt/n/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq"
REFERENCE_DIR="/mnt/n/refdata-gex-GRCh38-2024-A"
SAMPLE_ID="MFS-1"
transcript fasta = "N:/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
genes.gtf = "N:/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/genes/genes.gtf"


FASTQ_DIR="C:/MFS/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq"
REFERENCE_DIR="C:/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A"
SAMPLE_ID="MFS-1"
genome fasta ="C:\MFS\refdata-gex-GRCh38-2024-A\refdata-gex-GRCh38-2024-A\fasta/genome.fa"
genes.gtf = "C:/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
READ1="C:/MFS/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq/MFS-1-1_S9_L003_R1_001.fastq.gz"
READ2="C:/MFS/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq/MFS-1-1_S9_L003_R2_001.fastq.gz"
transcriptome_output.fa = “C:/MFS/transcriptome_output.fa"
cellranger ="C:\MFS\cellranger_extracted\cellranger-9.0.1\bin/cellranger"

kb-count \
--h5 \
--workflow standard \
--transcriptome "$REFERENCE_DIR/transcripts.fa" \
--genes "$REFERENCE_DIR/refdata-gex-GRCh38-2024-A/genes/genes.gtf" \
--out "$SAMPLE_ID"_kb_output \
--reads "$FASTQ_DIR/MFS-1-1_S9_L003_R1_001.fastq.gz" \
"$FASTQ_DIR/MFS-1-1_S9_L003_R2_001.fastq.gz"



sudo apt-get update
sudo apt-get install kallisto

# Convert Windows paths to Linux paths (assuming /mnt/c/ for Windows C: drive)
FASTQ_DIR="/mnt/c/MFS/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq"
REF_DIR="/mnt/c/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A"
SAMPLE_ID="MFS-1"

TRANSCRIPTOME_FA="/mnt/c/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
# For Kallisto, you need transcriptome FASTA, typically mRNA transcripts, not genome.
# If you only have genome.fa, you need transcriptome sequences. 
# Assuming you have transcript FASTA as transcript.fasta, set its path accordingly.
TRANSCRIPTOME_TXT="/mnt/c/MFS/transcriptome_output.fa"  # Path to transcript FASTA

READ1="${FASTQ_DIR}/${SAMPLE_ID}-1_S9_L003_R1_001.fastq.gz"
READ2="${FASTQ_DIR}/${SAMPLE_ID}-1_S9_L003_R2_001.fastq.gz"

# 1. Build Kallisto index (once)
kallisto index -i kallisto_index.idx "$TRANSCRIPTOME_TXT"

# 2. Quantify transcripts from paired-end reads
kallisto quant -i kallisto_index.idx -o ${SAMPLE_ID}_kallisto_output \
--threads=8 \
"${READ1}" "${READ2}"

# 3. Generate count matrix (transcript-level counts)
# The output is in: ${SAMPLE_ID}_kallisto_output/abundance.tsv
# To get counts:
cut -f1,5 ${SAMPLE_ID}_kallisto_output/abundance.tsv > ${SAMPLE_ID}_counts.tsv





simpleaf index \
  --output "/mnt/c/MFS/index" \
  --fasta "/mnt/c/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/fasta/genome.fa" \
  --gtf "/mnt/c/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/genes/genes.gtf" \
  --decoy-paths "none" \
  --ref-type splici \
  --use-piscem \
  --kmer-length 31 \
  --threads 8



simpleaf index \

  --fasta "/mnt>   --output "/mnt/c/MFS/index" \
--fasta "/mnt/c/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/fasta/genome.fa" \
--gtf "/mnt/c/MFS/refdata-gex-GRCh38-2024-A/refdata-gex-GRCh38-2024-A/genes/genes.gtf" \
--decoy-paths "/mnt/c/MFS/index/empty_decoy.txt" \
--ref-type splic>   --ref-type splici \
--use-piscem \
--kmer-length 31 \
--threads 12

nohup simpleaf index \
--output "/mnt/c/MFS/index" \
--fasta "/mnt/c/MFS/refdata-gex-GRCh38-2024-A/fasta/genome.fa" \
--gtf "/mnt/c/MFS/refdata-gex-GRCh38-2024-A/genes/genes.gtf" \
--decoy-paths "/mnt/c/MFS/index/empty_decoy.txt" \
--ref-type splici \
--use-piscem \
--kmer-length 31 \
--threads 12 \
> simpleaf_index.log 2>&1 &
  
  simpleaf quant \
--input /mnt/c/MFS/fastq/sample1 \
--output /mnt/c/MFS/quant/sample1 \
--index /mnt/c/MFS/index \
--threads 12 \

simpleaf quant \
t/c/MFS/MFS-1/Ho>   --reads1 "/mnt/c/MFS/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq/MFS-1-1_S9_L003_R1_001.fastq.gz" \
>   --reads2 "/mnt/c/MFS/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq/MFS-1-1_S9_L003_R2_001.fastq.gz" \
>   --output "/mnt/c/MFS/quant/sample1" \
>   --index "/mnt/c/MFS/index/index/piscem_idx" \
>   --t2g-map "/mnt/c/MFS/index/ref/t2g_3col.tsv" \
-chemistry 10xv3>   --chemistry 10xv3 \
>   --resolution cr-like \
>   --knee \
ds 12
>   --threads 12


simpleaf quant \
--reads1 "/mnt/c/MFS/MFS-2/HoJWK_10X3RNA_CPOS-230905-ZNS-19640a/primary_seq/MFS-2-1_S3_L002_R1_001.fastq.gz" \
--reads2 "/mnt/c/MFS/MFS-2/HoJWK_10X3RNA_CPOS-230905-ZNS-19640a/primary_seq/MFS-2-1_S3_L002_R2_001.fastq.gz" \
--output "/mnt/d/移动云盘同步盘/Y4S1/FYP/sample2" \
--index "/mnt/c/MFS/index/index/piscem_idx" \
--t2g-map "/mnt/c/MFS/index/ref/t2g_3col.tsv" \
--chemistry 10xv3 \
--resolution cr-like \
--knee \
--threads 12





FASTQ_DIR <- "N:/mnt/n/MFS-1/HoJWK_10X3RNA_CPOS-230829-ZNS-19507a/primary_seq"
REFERENCE_DIR <- "N:/mnt/n/refdata-gex-GRCh38-2024-A"
SAMPLE_ID <- "MFS-1"


kallisto_path <- "N:/Kallisto/kallisto.exe"
bustools_path <- "N:/bustools/bustools.exe"  # 如需用到后续流程

# 输入文件路径
read1 <- file.path(FASTQ_DIR, "MFS-1-1_S9_L003_R1_001.fastq.gz")
read2 <- file.path(FASTQ_DIR, "MFS-1-1_S9_L003_R2_001.fastq.gz")
transcriptome <- file.path(REFERENCE_DIR, "transcripts.fa")
genes_gtf <- file.path(REFERENCE_DIR, "refdata-gex-GRCh38-2024-A/genes/genes.gtf")
out_dir <- paste0(SAMPLE_ID, "_kb_output")

# 构建参数列表
args <- c(
  "count",
  "--h5",
  "--workflow", "snrna",
  "--transcriptome", transcriptome,
  "--genes", genes_gtf,
  "--out", out_dir,
  "--reads", read1,
  read2
)

# 调用系统命令
system2(kallisto_path, args = args)





####QC/normalization(preprocessing), dimesion reduction/clustering, comparative analysis
install.packages("reticulate")
library(reticulate)

# Use a specific Python environment
use_python("/usr/bin/python3")  # or use_virtualenv(), use_condaenv()

# Run Python code inline
py_run_string("import numpy as np; x = np.arange(10)")
py$x  # Access Python variable from R

where python


use_python("D:/download/miniconda/python.exe")
py_run_string("print('Hello from Python')")



packageVersion("cli")  # Should show 3.6.5 or newer
install.packages("devtools")
install.packages("remotes")
install_github("chris-mcginnis-ucsf/DoubletFinder")


# Load libraries
library(reticulate)

library(Seurat)
library(Matrix)
library(dplyr)
library(devtools)
library(DoubletFinder)


# Load MFS1
mfs1_counts <- readMM("C:/MFS/quant/sample1/af_quant/alevin/quants_mat.mtx")
mfs1_features <- readLines("C:/MFS/quant/sample1/af_quant/alevin/quants_mat_cols.txt")
mfs1_barcodes <- readLines("C:/MFS/quant/sample1/af_quant/alevin/quants_mat_rows.txt")
# Transpose the matrix
mfs1_counts <- t(mfs1_counts)
rownames(mfs1_counts) <- mfs1_features
colnames(mfs1_counts) <- mfs1_barcodes
mfs1 <- CreateSeuratObject(counts = mfs1_counts, project = "MFS1")



dim(mfs1_counts)
length(mfs1_features)
length(mfs1_barcodes)

setwd("D:/移动云盘同步盘/Y4S1/FYP")
library(Matrix)

# Assume mfs1_counts is your sparse matrix
dense_mtx <- as.matrix(mfs1_counts)

# Inspect dimensions
dim(dense_mtx)

# Example: look at the first few rows and columns
dense_mtx[1:5, 1:5]

write.csv(dense_mtx, file = "mfs1_counts_dense.csv")


library(Matrix)
library(org.Hs.eg.db)   # Human gene annotation
library(dplyr)

# Build mapping from Ensembl ID to gene symbol
id2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = rownames(mfs1_counts),
                                   keytype = "ENSEMBL",
                                   columns = c("SYMBOL"))

# Create lookup vector
id2symbol_vec <- setNames(id2symbol$SYMBOL, id2symbol$ENSEMBL)

# Build triplet dataframe
triplet_df <- data.frame(
  gene_id   = rownames(mfs1_counts)[mfs1_counts@i + 1],
  gene_name = id2symbol_vec[rownames(mfs1_counts)[mfs1_counts@i + 1]],
  cell      = colnames(mfs1_counts)[mfs1_counts@j + 1],
  count     = mfs1_counts@x
)

# Write to CSV
write.csv(triplet_df, "mfs1_counts_triplet_with_names.csv", row.names = FALSE)


library(Matrix)
library(org.Hs.eg.db)   # Human gene annotation
library(dplyr)

# Build mapping from Ensembl ID to gene symbol
id2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = rownames(mfs2_counts),
                                   keytype = "ENSEMBL",
                                   columns = c("SYMBOL"))

# Create lookup vector
id2symbol_vec <- setNames(id2symbol$SYMBOL, id2symbol$ENSEMBL)

# Build triplet dataframe
triplet_df <- data.frame(
  gene_id   = rownames(mfs2_counts)[mfs2_counts@i + 1],
  gene_name = id2symbol_vec[rownames(mfs2_counts)[mfs2_counts@i + 1]],
  cell      = colnames(mfs2_counts)[mfs2_counts@j + 1],
  count     = mfs2_counts@x
)

# Write to CSV
write.csv(triplet_df, "mfs2_counts_triplet_with_names.csv", row.names = FALSE)





setwd("D:/移动云盘同步盘/Y4S1/FYP“)

library(rhdf5)

library(rhdf5)






# Load MFS2
mfs2_counts <- readMM("D:/移动云盘同步盘/Y4S1/FYP/sample2/af_quant/alevin/quants_mat2.mtx")
mfs2_features <- readLines("D:/移动云盘同步盘/Y4S1/FYP/sample2/af_quant/alevin/quants_mat2_cols.txt")
mfs2_barcodes <- readLines("D:/移动云盘同步盘/Y4S1/FYP/sample2/af_quant/alevin/quants_mat2_rows.txt")
# Transpose the matrix
mfs2_counts <- t(mfs2_counts)
rownames(mfs2_counts) <- mfs2_features
colnames(mfs2_counts) <- mfs2_barcodes
mfs2 <- CreateSeuratObject(counts = mfs2_counts, project = "MFS2")




# Add sample labels
mfs1$sample <- "MFS1"
mfs2$sample <- "MFS2"

# Quality Control: calculate mitochondrial percentage
mfs1[["percent.mt"]] <- PercentageFeatureSet(mfs1, pattern = "^MT-")
mfs2[["percent.mt"]] <- PercentageFeatureSet(mfs2, pattern = "^MT-")

# Filter cells: keep cells with >200 features and <10% mitochondrial reads
mfs1 <- subset(mfs1, subset = nFeature_RNA > 200 & percent.mt < 10)
mfs2 <- subset(mfs2, subset = nFeature_RNA > 200 & percent.mt < 10)

##MFS1: Filtered 5 cells out of 8435 
##MFS2: Filtered 24 cells out of 10794 

# Before filtering
n_before_mfs1 <- ncol(mfs1_counts)
n_before_mfs2 <- ncol(mfs2_counts)

# After filtering
n_after_mfs1 <- ncol(mfs1)
n_after_mfs2 <- ncol(mfs2)

# Number of cells removed
filtered_mfs1 <- n_before_mfs1 - n_after_mfs1
filtered_mfs2 <- n_before_mfs2 - n_after_mfs2

# Print results
cat("MFS1: Filtered", filtered_mfs1, "cells out of", n_before_mfs1, "\n")
cat("MFS2: Filtered", filtered_mfs2, "cells out of", n_before_mfs2, "\n")


# Normalize and find variable features
mfs1 <- NormalizeData(mfs1) %>% FindVariableFeatures()
mfs2 <- NormalizeData(mfs2) %>% FindVariableFeatures()

#How Log Normalization Works in Seurat
#Step 1: For each cell, counts are divided by the total counts for that cell.
#Step 2: The result is multiplied by a scale factor (default = 10,000).
#Step 3: A natural log transformation is applied: log-normalized value = log(1+scaled count)
#log-normalized value=log(1+scaled count)

# MFS1 preprocessing
mfs1 <- NormalizeData(mfs1)
mfs1 <- FindVariableFeatures(mfs1)
mfs1 <- ScaleData(mfs1)
mfs1 <- RunPCA(mfs1)
mfs1 <- FindNeighbors(mfs1, dims = 1:20)
mfs1 <- FindClusters(mfs1, resolution = 0.5)
mfs1 <- RunUMAP(mfs1, dims = 1:20)

# MFS2 preprocessing
mfs2 <- NormalizeData(mfs2)
mfs2 <- FindVariableFeatures(mfs2)
mfs2 <- ScaleData(mfs2)
mfs2 <- RunPCA(mfs2)
mfs2 <- FindNeighbors(mfs2, dims = 1:20)
mfs2 <- FindClusters(mfs2, resolution = 0.5)
mfs2 <- RunUMAP(mfs2, dims = 1:20)

#MFS1 top markers 
#markers_mfs1 <- FindAllMarkers(mfs1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top_genes_mfs1 <- markers_mfs1 %>%
#  group_by(cluster) %>%
#  slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE) %>%
#  select(cluster, gene, avg_log2FC, p_val_adj)

#markers_mfs2 <- FindAllMarkers(mfs2_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top_genes_mfs2 <- markers_mfs2 %>%
#  group_by(cluster) %>%
#  slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE) %>%
#  select(cluster, gene, avg_log2FC, p_val_adj)

#top_genes_mfs1$group <- "MFS1"
#top_genes_mfs2$group <- "MFS2"

#combined_top_genes <- bind_rows(top_genes_mfs1, top_genes_mfs2) %>%
#  arrange(group, cluster)

#DotPlot(mfs1, features = top_genes_mfs1$gene) + RotatedAxis()
#DotPlot(mfs2_clean, features = top_genes_mfs2$gene) + RotatedAxis()

#DoHeatmap(mfs1, features = top_genes_mfs1$gene)
#DoHeatmap(mfs2_clean, features = top_genes_mfs2$gene)


#MFS2 top markers
markers_mfs2 <- FindAllMarkers(mfs2_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_genes_mfs2 <- markers_mfs2 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE)

cluster_summary <- top_genes_mfs2 %>%
  select(cluster, gene, avg_log2FC, p_val_adj) %>%
  arrange(cluster)

DoHeatmap(mfs2_clean, features = top_genes_mfs2$gene)
DotPlot(mfs2_clean, features = top_genes_mfs2$gene) + RotatedAxis()

##new cluster identification
# Load required libraries
library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

install.packages("biomaRt")

# Step 1: Get top 10 marker genes per cluster
mfs1_clean <- NormalizeData(mfs1_clean)
mfs1_clean <- FindVariableFeatures(mfs1_clean)
mfs1_clean <- ScaleData(mfs1_clean)
mfs1_clean <- RunPCA(mfs1_clean)
mfs1_clean <- FindNeighbors(mfs1_clean, dims = 1:20)
mfs1_clean <- FindClusters(mfs1_clean, resolution = 0.5)

mfs2_clean <- NormalizeData(mfs2_clean)
mfs2_clean <- FindVariableFeatures(mfs2_clean)
mfs2_clean <- ScaleData(mfs2_clean)
mfs2_clean <- RunPCA(mfs2_clean)
mfs2_clean <- FindNeighbors(mfs2_clean, dims = 1:20)
mfs2_clean <- FindClusters(mfs2_clean, resolution = 0.5)

markers_mfs1 <- FindAllMarkers(mfs1_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_mfs1 <- markers_mfs1 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)

library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Get top 50 markers per cluster
top50_mfs1 <- markers_mfs1 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE)

top50_mfs1 <- top50_mfs1 %>%
  ungroup() %>%
  dplyr::select(-gene_symbol)
head(top50_mfs1$gene)



# Remove "-U" suffix from Ensembl IDs
top50_mfs1$gene_clean <- gsub("-U$", "", top50_mfs1$gene)

top50_mfs1$gene_symbol <- mapIds(org.Hs.eg.db,
                                 keys = top50_mfs1$gene_clean,  # cleaned IDs
                                 column = "SYMBOL",             # target: gene symbol
                                 keytype = "ENSEMBL",           # input type
                                 multiVals = "first")

top50_mfs1$gene <- ifelse(is.na(top50_mfs1$gene_symbol),
                          top50_mfs1$gene, top50_mfs1$gene_symbol)

head(top50_mfs1$gene)

# Add gene description (function/annotation)
top50_mfs1$gene_function <- mapIds(org.Hs.eg.db,
                                   keys = top50_mfs1$gene_clean,  # cleaned Ensembl IDs
                                   column = "GENENAME",           # target: gene description
                                   keytype = "ENSEMBL",           # input type
                                   multiVals = "first")

markers_mfs2 <- markers_mfs2 %>%
  ungroup() %>%
  dplyr::select(-gene_symbol)

# Step 1: Clean IDs
markers_mfs2$gene_clean <- gsub("-U$", "", markers_mfs2$gene)

# Step 2: Map to gene symbols
markers_mfs2$gene_symbol <- mapIds(org.Hs.eg.db,
                                   keys = markers_mfs2$gene_clean,
                                   column = "SYMBOL",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")

# Step 3: Replace gene column with symbols (keep original if NA)
markers_mfs2$gene <- ifelse(is.na(markers_mfs2$gene_symbol),
                            markers_mfs2$gene, markers_mfs2$gene_symbol)
# Step 4: Drop helper columns
markers_mfs2 <- markers_mfs2 %>% dplyr::select(-gene_clean, -gene_symbol)

head(markers_mfs2$gene)

head(top50_mfs1$gene)

markers_mfs2$gene_clean <- gsub("-U$", "", markers_mfs2$gene)
# Remove -U or -A suffixes
markers_mfs2$gene_clean <- gsub("-[A-Z]$", "", markers_mfs2$gene)

# Add gene description (function/annotation)
markers_mfs2$gene_function <- mapIds(org.Hs.eg.db,
                                   keys = markers_mfs2$gene_clean,  # cleaned Ensembl IDs
                                   column = "GENENAME",           # target: gene description
                                   keytype = "ENSEMBL",           # input type
                                   multiVals = "first")

head(markers_mfs2[, c("gene_clean", "gene_function")])


markers_mfs2 <- FindAllMarkers(mfs2_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top50_mfs2 <- markers_mfs2 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE)


top50_mfs2 <- top50_mfs2 %>%
  ungroup() %>%
  dplyr::select(-gene_symbol)
head(top50_mfs1$gene)


# Remove "-U" suffix from Ensembl IDs
top50_mfs2$gene_clean <- gsub("-U$", "", top50_mfs2$gene)
# Add gene description (function/annotation)
top50_mfs2$gene_function <- mapIds(org.Hs.eg.db,
                                   keys = top50_mfs2$gene_clean,  # cleaned Ensembl IDs
                                   column = "GENENAME",           # target: gene description
                                   keytype = "ENSEMBL",           # input type
                                   multiVals = "first")

write.csv(markers_mfs1, file = "markers_mfs1_new2.csv", row.names = FALSE)
write.csv(markers_mfs2, file = "markers_mfs2_new3.csv", row.names = FALSE)

write.csv(top50_mfs1, file = "top50_mfs1_new.csv", row.names = FALSE)
write.csv(top50_mfs2, file = "top50_mfs2_new.csv", row.names = FALSE)


# Step 2: Clean Ensembl IDs and map to gene symbols
top10_mfs1$ensembl_clean <- gsub("-[A-Z]+$", "", top10_mfs1$gene)
top10_mfs2$ensembl_clean <- gsub("-[A-Z]+$", "", top10_mfs2$gene)

library(AnnotationDbi)
install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

top10_mfs1$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = top10_mfs1$ensembl_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

top10_mfs2$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = top10_mfs2$ensembl_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Step 3: Add cell type annotations (customize as needed)
cell_types_mfs1 <- c("0" = "T cells", "1" = "NK cells", "2" = "Monocytes", "3" = "B cells",
                     "4" = "Dendritic cells", "5" = "CD8+ T cells", "6" = "Tregs", "7" = "Plasma cells",
                     "8" = "Macrophages", "9" = "pDCs", "10" = "cDCs", "11" = "ILC", "12" = "Mast cells",
                     "13" = "Erythrocytes", "14" = "Progenitors")

cell_types_mfs2 <- c("0" = "T cells", "1" = "NK cells", "2" = "Monocytes", "3" = "B cells",
                     "4" = "Dendritic cells", "5" = "CD8+ T cells", "6" = "Tregs", "7" = "Plasma cells",
                     "8" = "Macrophages", "9" = "pDCs", "10" = "cDCs", "11" = "ILC")

top10_mfs1$cell_type <- cell_types_mfs1[as.character(top10_mfs1$cluster)]
top10_mfs2$cell_type <- cell_types_mfs2[as.character(top10_mfs2$cluster)]

# Step 4: Combine and view
top10_mfs1$group <- "MFS1"
top10_mfs2$group <- "MFS2"

combined_markers <- bind_rows(top10_mfs1, top10_mfs2) %>%
  select(group, cluster, gene_symbol, gene, cell_type, avg_log2FC, p_val_adj) %>%
  arrange(group, cluster)

View(combined_markers)

DimPlot(mfs1, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  ggtitle("MFS1 Clusters")

DimPlot(mfs2_clean, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  ggtitle("MFS2 Clusters")



#Estimate expected doublets
nExp_mfs1 <- round(0.05 * ncol(mfs1))  # 5% expected doublets
nExp_mfs2 <- round(0.05 * ncol(mfs2))

#(fail)Run DoubletFinder
# Sweep parameters to find optimal pK
sweep.res <- paramSweep(mfs1, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)

best.pK <- as.character(sweep.stats$pK[which.max(sweep.stats$BCreal)])
best.pK <- as.numeric(best.pK)
best.pK

# fail running doubletfinder


# Step 1: Run paramSweep and summarize
sweep.res <- paramSweep(mfs1, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)

# Step 2: Estimate homotypic doublet proportion
homotypic.prop <- modelHomotypic(mfs1$seurat_clusters)

# Step 3: Adjust expected doublets
nExp_poi <- round(0.05 * ncol(mfs1))  # raw expected
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))  # adjusted for homotypic

#library(DoubletFinder)
#ls("package:DoubletFinder")
# if(!require(DoubletFinder))
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

  require(Matrix)
  require(Seurat)
  require(fields)
  require(KernSmooth)
  require(cowplot)
  require(ggplot2)
  require(RANN)
  


# Find the correct pK name
pk.name <- paste("pANN", best.pK, sep = "_")


pANN <- sweep.res[[1]]$pANN
names(pANN) <- rownames(sweep.res[[1]])

# Initialize all as singlets
classification <- rep("Singlet", length(pANN))
names(classification) <- names(pANN)

# Identify top-scoring cells
doublet.ids <- names(sort(pANN, decreasing = TRUE))[1:nExp_poi.adj]
classification[doublet.ids] <- "Doublet"

# Convert to data frame
classifications <- data.frame(classification)
colnames(classifications) <- paste("DF.classifications", 0.25, best.pK, nExp_poi.adj, sep = "_")

#Add to Seurat object
mfs1 <- AddMetaData(mfs1, metadata = classifications)

table(mfs1[[colnames(classifications)]])

##DF.classifications_0.25__381
#Doublet Singlet 
#381    8049 

mfs1_clean <- subset(mfs1, subset = DF.classifications_0.25_0.005_381 == "Singlet")


#Count matrix for mfs1 
# Load required packages
install.packages("openxlsx")
library(Seurat)
library(openxlsx)   # for Excel export

# Install if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")

library(Seurat)
library(Matrix)
library(rhdf5)

# Subset your Seurat object (already done in your case)
mfs1_clean <- subset(mfs1, subset = DF.classifications_0.25_0.005_381 == "Singlet")

# Extract raw counts from the RNA assay
mfs1_countmatrix <- GetAssayData(mfs1_clean, assay = "RNA", layer = "counts")





# Define output file
h5_file <- "mfs1_countmatrix.h5"

# Create HDF5 file
h5createFile(h5_file)

# Save sparse matrix components
h5write(mfs1_countmatrix@i, h5_file, "row_indices")       # row indices of non-zero entries
h5write(mfs1_countmatrix@p, h5_file, "col_pointers")      # column pointers
h5write(mfs1_countmatrix@x, h5_file, "values")            # non-zero values

# Save dimensions and labels
h5write(dim(mfs1_countmatrix), h5_file, "dimensions")
h5write(rownames(mfs1_countmatrix), h5_file, "genes")
h5write(colnames(mfs1_countmatrix), h5_file, "cells")

# Close file
H5close()

library(Matrix)
library(rhdf5)

# Define file
h5_file <- "mfs1_countmatrix.h5"

# Read components
row_indices <- h5read(h5_file, "row_indices")
col_pointers <- h5read(h5_file, "col_pointers")
values       <- h5read(h5_file, "values")
dims         <- h5read(h5_file, "dimensions")
genes        <- h5read(h5_file, "genes")
cells        <- h5read(h5_file, "cells")

# Reconstruct sparse matrix (dgCMatrix)
mfs1_countmatrix_loaded <- new("dgCMatrix",
                               i = as.integer(row_indices),
                               p = as.integer(col_pointers),
                               x = as.numeric(values),
                               Dim = as.integer(dims),
                               Dimnames = list(genes, cells))

# Check
mfs1_countmatrix_loaded[1:5, 1:5]

h5ls(h5_file)

library(Matrix)

# Assume you already reconstructed the sparse matrix:
# mfs1_countmatrix_loaded <- new("dgCMatrix", ...)

# Assuming your Seurat object is called mfs1_clean
gene_names <- rownames(mfs1_clean)

library(Matrix)

# Triplet form of the sparse matrix
triplet <- summary(mfs1_countmatrix_loaded)

# Build clean data frame using Seurat's gene names
triplet_df <- data.frame(
  Gene  = gene_names[triplet$i],
  Cell  = colnames(mfs1_countmatrix_loaded)[triplet$j],
  Count = triplet$x
)

# Export
write.csv(triplet_df, "mfs1_countmatrix_nonzero1.csv", row.names = FALSE)

head(gene_names, 200)

sum(grepl("^\\.", gene_names))

sum(grepl("^LOC", gene_names))

sum(grepl("^ENSG", gene_names))


##mfs2 count matrix 
# Extract raw counts from the RNA assay
mfs2_countmatrix <- GetAssayData(mfs2_clean, assay = "RNA", layer = "counts")





# Define output file
h5_file2 <- "mfs2_countmatrix.h5"

# Create HDF5 file
h5createFile(h5_file2)

# Save sparse matrix components
h5write(mfs2_countmatrix@i, h5_file2, "row_indices")       # row indices of non-zero entries
h5write(mfs2_countmatrix@p, h5_file2, "col_pointers")      # column pointers
h5write(mfs2_countmatrix@x, h5_file2, "values")            # non-zero values

# Save dimensions and labels
h5write(dim(mfs2_countmatrix), h5_file2, "dimensions")
h5write(rownames(mfs2_countmatrix), h5_file2, "genes")
h5write(colnames(mfs2_countmatrix), h5_file2, "cells")

# Close file
H5close()

library(Matrix)
library(rhdf5)

# Define file
h5_file2 <- "mfs2_countmatrix.h5"

# Read components
row_indices <- h5read(h5_file2, "row_indices")
col_pointers <- h5read(h5_file2, "col_pointers")
values       <- h5read(h5_file2, "values")
dims         <- h5read(h5_file2, "dimensions")
genes        <- h5read(h5_file2, "genes")
cells        <- h5read(h5_file2, "cells")

# Reconstruct sparse matrix (dgCMatrix)
mfs2_countmatrix_loaded <- new("dgCMatrix",
                               i = as.integer(row_indices),
                               p = as.integer(col_pointers),
                               x = as.numeric(values),
                               Dim = as.integer(dims),
                               Dimnames = list(genes, cells))

# Check
mfs1_countmatrix_loaded[1:5, 1:5]

h5ls(h5_file)

library(Matrix)

# Assume you already reconstructed the sparse matrix:
# mfs1_countmatrix_loaded <- new("dgCMatrix", ...)

# Assuming your Seurat object is called mfs1_clean
gene_names2 <- rownames(mfs2_clean)

library(Matrix)

library(AnnotationDbi)
library(org.Hs.eg.db)

# Map Ensembl IDs to HGNC symbols
gene_symbols2 <- mapIds(org.Hs.eg.db,
                        keys = gene_names2,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Triplet form of the sparse matrix
triplet2 <- summary(mfs2_countmatrix_loaded)

# Build triplet with both ID and symbol
triplet_df2 <- data.frame(
  GeneID = gene_names2[triplet2$i],                          # original ID (may be ENSG or LOC)
  Symbol = gene_symbols2[gene_names2[triplet2$i]],           # mapped HGNC symbol
  Cell   = colnames(mfs2_countmatrix_loaded)[triplet2$j],    # cell barcodes
  Count  = triplet2$x                                        # counts
)

write.csv(triplet_df2, "mfs2_countmatrix_nonzero_with_symbols.csv", row.names = FALSE)







library(Seurat)
library(ggplot2)


#MFS2 doubletfinding
sweep.res.mfs2 <- paramSweep(mfs2, PCs = 1:20, sct = FALSE)
sweep.stats.mfs2 <- summarizeSweep(sweep.res.mfs2, GT = FALSE)

#Find best pK
best.pK.mfs2 <- as.numeric(as.character(sweep.stats.mfs2[which.max(sweep.stats.mfs2$BCreal), "pK"]))

#Estimate expected doublets
nExp.mfs2 <- round(0.075 * ncol(mfs2))  # Adjust 0.075 if needed

#Extract pANN scores
pANN.mfs2 <- sweep.res.mfs2[[1]]$pANN
names(pANN.mfs2) <- rownames(sweep.res.mfs2[[1]])

#Classify doublets manually
classification.mfs2 <- rep("Singlet", length(pANN.mfs2))
names(classification.mfs2) <- names(pANN.mfs2)

doublet.ids.mfs2 <- names(sort(pANN.mfs2, decreasing = TRUE))[1:nExp.mfs2]
classification.mfs2[doublet.ids.mfs2] <- "Doublet"

classifications.mfs2 <- data.frame(classification.mfs2)
colnames(classifications.mfs2) <- paste("DF.classifications", 0.25, best.pK.mfs2, nExp.mfs2, sep = "_")

#add to Seurat object 
mfs2 <- AddMetaData(mfs2, metadata = classifications.mfs2)


mfs2_clean <- subset(mfs2, subset = DF.classifications_0.25__808 == "Singlet")

sum(is.na(mfs2@meta.data$DF.classifications_0.25__808))

# Find overlapping cells
common.cells <- intersect(names(pANN.mfs2), colnames(mfs2))

# Reclassify only those cells
classification.mfs2 <- rep("Singlet", length(common.cells))
names(classification.mfs2) <- common.cells

doublet.ids.mfs2 <- names(sort(pANN.mfs2[common.cells], decreasing = TRUE))[1:nExp.mfs2]
classification.mfs2[doublet.ids.mfs2] <- "Doublet"

# Create metadata
classifications.mfs2 <- data.frame(classification.mfs2)
colnames(classifications.mfs2) <- paste("DF.classifications", 0.25, best.pK.mfs2, nExp.mfs2, sep = "_")

# Add to Seurat object
mfs2 <- AddMetaData(mfs2, metadata = classifications.mfs2)

classification_col.mfs2 <- colnames(classifications.mfs2)[1]
mfs2_clean <- subset(mfs2, subset = DF.classifications_0.25_0.01_808 == "Singlet")

table(mfs2@meta.data$DF.classifications_0.25_0.01_808)
#Doublet Singlet 
808    9192 




setwd("D:/移动云盘同步盘/Y4S1/FYP/sample2/af_quant/alevin")

library(Matrix)

# Read the Matrix Market file
mtx <- readMM("quants_mat2.mtx")

# Number of cells (columns)
nrow(mtx)

# Highest cell index actually present
max(which(rowSums(mtx) > 0))


# How many cells in the object?
ncol(mfs2)

# How many cells have a classification?
table(mfs2@meta.data$DF.classifications_0.25_0.01_808, useNA = "ifany")

# How many NA classifications?
sum(is.na(mfs2@meta.data$DF.classifications_0.25_0.01_808))






#UMAP visualization
mfs2 <- RunUMAP(mfs2, dims = 1:20)

classification_col.mfs2 <- grep("^DF\\.classifications", colnames(mfs2@meta.data), value = TRUE)

DimPlot(mfs2, group.by = classification_col.mfs2, reduction = "umap") +
  ggtitle("UMAP: DoubletFinder Classification (mfs2)") +
  theme_minimal()





# Extract the classification column name
classification_col <- grep("^DF\\.classifications", colnames(mfs1@meta.data), value = TRUE)

# Plot UMAP colored by doublet status
DimPlot(mfs1, group.by = classification_col, reduction = "umap") +
  ggtitle("UMAP: DoubletFinder Classification") +
  theme_minimal()











#Use clean singlet objects
anchors <- FindIntegrationAnchors(object.list = list(mfs1_clean, mfs2_clean), dims = 1:30)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)

#Downstream analysis
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined)


DimPlot(combined, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP: Integrated Clusters") +
  theme_minimal()


#Compare cluster composition across samples
table(combined$seurat_clusters, combined$orig.ident)



#Identify marker genes

DefaultAssay(combined) <- "RNA"  # Switch to RNA for expression

combined <- JoinLayers(combined)

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

head(markers)

top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


#Compare clusters across samples
DimPlot(combined, reduction = "umap", group.by = "orig.ident")

setwd("D:/移动云盘同步盘/Y4S1/FYP")
write.csv(markers, file = "markers.csv", row.names = FALSE)

top_genes <- head(markers$gene, 10)  # Adjust number as needed

library(Seurat)
library(dplyr)


# Visualize top marker genes across samples
library(dplyr)



top_genes_by_cluster <- markers %>%
  mutate(score = avg_log2FC / (p_val_adj + 1e-10)) %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  arrange(cluster) %>%
  pull(gene)


#Calculate average expression per sample
avg_exp <- AverageExpression(combined, group.by = "orig.ident", return.seurat = TRUE)
avg_df <- as.data.frame(avg_exp$RNA[top_genes_by_cluster, ])

colnames(avg_df)

head(rownames(avg_exp$RNA))

print(top_genes_by_cluster)

# Check which genes have non-zero expression in either sample
nonzero_genes <- rownames(avg_exp$RNA)[rowSums(avg_exp$RNA) > 0]

# Filter valid genes again
valid_genes_filtered <- intersect(valid_genes, nonzero_genes)

# Extract expression
avg_df <- as.data.frame(avg_exp$RNA[valid_genes_filtered, ])



DefaultAssay(mfs1_clean) <- "RNA"
DefaultAssay(mfs2_clean) <- "RNA"



avg_mfs1 <- AverageExpression(mfs1_clean, return.seurat = TRUE)
avg_mfs2 <- AverageExpression(mfs2_clean, return.seurat = TRUE)

# Extract log-normalized expression matrix from the "data" layer
expr_mfs1 <- GetAssayData(avg_mfs1, assay = "RNA", slot = "data")
expr_mfs2 <- GetAssayData(avg_mfs2, assay = "RNA", slot = "data")

# Clean gene IDs
genes_to_check <- gsub("-[A-Z]+$", "", top_genes_by_cluster)

# Filter valid genes
valid_genes <- intersect(genes_to_check, rownames(expr_mfs1))
valid_genes <- intersect(valid_genes, rownames(expr_mfs2))

# Extract expression
df_mfs1 <- expr_mfs1[valid_genes, , drop = FALSE]
df_mfs2 <- expr_mfs2[valid_genes, , drop = FALSE]

# If there's only one column (e.g., one group), use it directly
avg_df <- data.frame(
  gene = valid_genes,
  MFS1 = df_mfs1[, 1],
  MFS2 = df_mfs2[, 1]
)

avg_df_filtered <- avg_df %>%
  filter(MFS1 > 0 | MFS2 > 0)

only_in_mfs1 <- avg_df_filtered %>%
  filter(MFS1 > 0.25 & MFS2 <= 0.01)

VlnPlot(combined, features = only_in_mfs1$gene, group.by = "orig.ident", pt.size = 0)

length(valid_genes)
head(avg_df_filtered)




# Clean suffixes if needed
genes_to_check <- gsub("-[A-Z]+$", "", genes_to_check)

# Filter valid genes
valid_genes <- intersect(genes_to_check, rownames(avg_mfs1$RNA))
valid_genes <- intersect(valid_genes, rownames(avg_mfs2$RNA))

# Extract expression
df_mfs1 <- avg_mfs1$RNA[valid_genes, , drop = FALSE]
df_mfs2 <- avg_mfs2$RNA[valid_genes, , drop = FALSE]

# Combine
avg_df <- data.frame(
  gene = valid_genes,
  MFS1 = df_mfs1[, 1],
  MFS2 = df_mfs2[, 1]
)

str(avg_mfs1$RNA)
str(avg_mfs2$RNA)

avg_df_filtered <- avg_df %>%
  filter(MFS1 > 0 | MFS2 > 0)

only_in_mfs1 <- avg_df_filtered %>%
  filter(MFS1 > 0.25 & MFS2 <= 0.25)

VlnPlot(combined, features = only_in_mfs1$gene, group.by = "orig.ident", pt.size = 0)




#4. export results
write.csv(avg_exp_df, "MFS1_vs_MFS2_gene_expression.csv")

#5.Identify sample-specific markers (optional)
Idents(combined) <- "orig.ident"
sample_markers <- FindMarkers(combined, ident.1 = "MFS1", ident.2 = "MFS2", only.pos = TRUE)
head(sample_markers)






# Install if not already installed
install.packages("BiocManager")
BiocManager::install(c("SingleR", "celldex", "Seurat", "SingleCellExperiment"))

# Load libraries
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

install.packages("digest", force = TRUE)
install.packages("reshape2", force = TRUE)
BiocManager::install("BiocParallel", force = TRUE)


# Load broader immune reference
ref <- celldex::BlueprintEncodeData()


# For mfs1_clean
sce1 <- as.SingleCellExperiment(mfs1_clean)

# For mfs2_clean
sce2 <- as.SingleCellExperiment(mfs2_clean)

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(AnnotationDbi)

# Convert Ensembl to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(sce1),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Apply gene symbols to rownames
rownames(sce1) <- gene_symbols

# Remove rows with NA gene names
sce1 <- sce1[!is.na(rownames(sce1)), ]

# Remove duplicated gene names (SingleR requires unique rownames)
sce1 <- sce1[!duplicated(rownames(sce1)), ]

# Annotate mfs1_clean
pred1 <- SingleR(test = sce1, ref = ref, labels = ref$label.main)
mfs1_clean$SingleR.labels <- pred1$labels

# Convert Ensembl IDs to gene symbols for sce2
gene_symbols2 <- mapIds(org.Hs.eg.db,
                        keys = rownames(sce2),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Apply gene symbols
rownames(sce2) <- gene_symbols2

# Clean up: remove NA and duplicated gene names
sce2 <- sce2[!is.na(rownames(sce2)), ]
sce2 <- sce2[!duplicated(rownames(sce2)), ]

# Annotate mfs2_clean
pred2 <- SingleR(test = sce2, ref = ref, labels = ref$label.main)
mfs2_clean$SingleR.labels <- pred2$labels

# UMAP with SingleR labels
DimPlot(mfs1_clean, group.by = "SingleR.labels", label = TRUE)
DimPlot(mfs2_clean, group.by = "SingleR.labels", label = TRUE)

DimPlot(mfs1_clean, group.by = "SingleR.labels", label = FALSE)
DimPlot(mfs2_clean, group.by = "SingleR.labels", label = FALSE)

library(ggplot2)

# Plot colored by SingleR cell type (legend = cell types) 
p <- DimPlot(mfs1_clean, group.by = "SingleR.labels", label = FALSE) + theme(legend.position = "right") 
# Overlay cluster numbers using Seurat Idents 
p <- LabelClusters(DimPlot(mfs1_clean, group.by = "ident", label = FALSE), id = "ident", repel = TRUE) + theme(legend.position = "right") 
p

library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)

# 1. Create a combined cluster + cell type label (unique per cluster)
mfs1_clean$cluster_celltype <- paste0(Idents(mfs1_clean), ": ", mfs1_clean$SingleR.labels)

# 2. Plot colored by cluster-celltype combo (legend = cluster: cell type)
p <- DimPlot(mfs1_clean,
             group.by = "cluster_celltype",
             label = FALSE) + 
  theme(legend.position = "right")

# 3. Extract UMAP coordinates and cluster IDs
umap_df <- as.data.frame(mfs1_clean@reductions$umap@cell.embeddings)
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- Idents(mfs1_clean)

# 4. Compute cluster centers
centers <- umap_df %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))

# 5. Add cluster numbers manually at cluster centers
p <- p + geom_text(data = centers,
                   aes(x = UMAP_1, y = UMAP_2, label = cluster),
                   color = "black", size = 5)

# 6. Display the plot
p

library(dplyr)
library(ggplot2)

# 1. Create a combined cluster + cell type label (unique per cluster)
mfs2_clean$cluster_celltype <- paste0(Idents(mfs2_clean), ": ", mfs2_clean$SingleR.labels)

# 2. Plot colored by cluster-celltype combo (legend = cluster: cell type)
p2 <- DimPlot(mfs2_clean,
              group.by = "cluster_celltype",
              label = FALSE) + 
  theme(legend.position = "right")

# 3. Extract UMAP coordinates and cluster IDs
umap_df2 <- as.data.frame(mfs2_clean@reductions$umap@cell.embeddings)
colnames(umap_df2)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df2$cluster <- Idents(mfs2_clean)

# 4. Compute cluster centers
centers2 <- umap_df2 %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))

# 5. Add cluster numbers manually at cluster centers
p2 <- p2 + geom_text(data = centers2,
                     aes(x = UMAP_1, y = UMAP_2, label = cluster),
                     color = "black", size = 5)

# 6. Display the plot
p2






# Save enlarged plot
p1 <- DimPlot(mfs1_clean, group.by = "SingleR.labels", label = TRUE, repel = TRUE) +
  theme(text = element_text(size = 14))  # Increase font size
ggsave("mfs1_clean_plot.png", plot = p1, width = 10, height = 8, dpi = 300)

p2 <- DimPlot(mfs2_clean, group.by = "SingleR.labels", label = TRUE, repel = TRUE) +
  theme(text = element_text(size = 14))
ggsave("mfs2_clean_plot.png", plot = p2, width = 10, height = 8, dpi = 300)

DimPlot(mfs1_clean, group.by = "SingleR.labels", label = TRUE, repel = TRUE)
DimPlot(mfs2_clean, group.by = "SingleR.labels", label = TRUE, repel = TRUE)


library(ggplot2)
# Extract UMAP coordinates and cluster info
umap_data <- Embeddings(mfs1_clean, "umap")
cluster_ids <- mfs1_clean$seurat_clusters

# Compute cluster centers
centers <- aggregate(umap_data, by = list(cluster = cluster_ids), FUN = mean)

# Rename columns to match UMAP axes
colnames(centers)[2:3] <- c("UMAP_1", "UMAP_2")

# Create UMAP plot colored by SingleR labels
p <- DimPlot(mfs1_clean, group.by = "SingleR.labels", label = FALSE)



# Add cluster number labels to the plot
p + geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = cluster), size = 5)


# Set identities to Seurat clusters
Idents(mfs1_clean) <- "seurat_clusters"

# Plot colored by Seurat clusters, labels repel
p <- DimPlot(mfs1_clean, group.by = "seurat_clusters", label = FALSE)
p <- LabelClusters(p, id = "seurat_clusters", repel = TRUE, size = 6)


library(Seurat)
library(ggplot2)

# 1. Keep Seurat cluster numbers as identities
Idents(mfs1_clean) <- "seurat_clusters"

# 2. Plot colored by SingleR labels (legend shows annotations)
p <- DimPlot(mfs1_clean, group.by = "SingleR.labels", label = FALSE)

# 3. Compute cluster centers directly from Seurat embeddings
umap_data <- Embeddings(mfs1_clean, "umap")
cluster_ids <- Idents(mfs1_clean)
centers <- aggregate(umap_data, by = list(cluster = cluster_ids), FUN = mean)
colnames(centers)[2:3] <- c("UMAP_1", "UMAP_2")

# 4. Overlay cluster numbers (0–14) at centers
p <- p + geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = cluster),
                   size = 6, fontface = "bold", color = "black")

# 5. Save high‑resolution figure
ggsave("umap_clusters_highres.png", plot = p,
       width = 12, height = 10, dpi = 600)


# Create UMAP plot colored by SingleR labels
p2 <- DimPlot(mfs2_clean, group.by = "SingleR.labels", label = FALSE)

# Extract UMAP coordinates and cluster info
umap_data2 <- Embeddings(mfs2_clean, "umap")
cluster_ids2 <- mfs2_clean$seurat_clusters

# Compute cluster centers
centers2 <- aggregate(umap_data2, by = list(cluster = cluster_ids2), FUN = mean)

# Rename columns to match UMAP axes
colnames(centers2)[2:3] <- c("UMAP_1", "UMAP_2")

# Add cluster number labels to the plot
p2 + geom_text(data = centers2, aes(x = UMAP_1, y = UMAP_2, label = cluster), size = 5)



# Optional: Compare to original clusters
table(mfs1_clean$seurat_clusters, mfs1_clean$SingleR.labels)

# Cluster-level annotation for mfs1_clean
cluster_pred1 <- SingleR(test = sce1, ref = ref, labels = ref$label.main, clusters = mfs1_clean$seurat_clusters)
mfs1_clean$cluster_annotation <- cluster_pred1$labels[mfs1_clean$seurat_clusters]

# Same for mfs2_clean
cluster_pred2 <- SingleR(test = sce2, ref = ref, labels = ref$label.main, clusters = mfs2_clean$seurat_clusters)
mfs2_clean$cluster_annotation <- cluster_pred2$labels[mfs2_clean$seurat_clusters]


#Blueprint Encode reference dataset (part of celldex), which provides curated immune and stromal cell type profiles derived from the Blueprint and ENCODE projects.

#This reference is widely used for immune cell annotation because it contains transcriptomic signatures for major immune lineages (T cells, B cells, NK cells, monocytes/macrophages, dendritic cells, etc.).


#subpopulation mfs1 
library(Seurat)
library(igraph)
set.seed(42)

# Subset
cl13 <- subset(mfs1_clean, idents = "13")

# Option A: SCTransform (robust for small n)
cl13 <- SCTransform(cl13, verbose = FALSE)
cl13 <- RunPCA(cl13, verbose = FALSE)

# Choose PCs (inspect ElbowPlot if needed)
ElbowPlot(cl13, ndims = 50)
use_dims <- 1:15

# Build graph with smaller k
cl13 <- FindNeighbors(cl13, dims = use_dims, k.param = 15)

# Try higher resolution; if too many tiny clusters, dial back
cl13 <- FindClusters(cl13, resolution = 1.2)  # 0.8–2.0 typical

# UMAP with tighter min.dist to encourage separation
cl13 <- RunUMAP(
  cl13, dims = use_dims,
  n.neighbors = 15, min.dist = 0.1,
  spread = 1.5, verbose = FALSE
)

DimPlot(cl13, reduction = "umap", group.by = "seurat_clusters",
        label = TRUE, repel = TRUE)

# Genes of interest
deg_genes <- c("ICOS","CD3G","CD2","CTLA4","PDCD1","CD28","TNFRSF9","CD96")

# Keep only genes that exist in the object
deg_genes_present <- deg_genes[deg_genes %in% rownames(cl13)]

# Extract expression values
expr_mat <- FetchData(cl13, vars = deg_genes_present)

# Compute average expression per cluster
avg_expr <- expr_mat %>%
  mutate(cluster = cl13$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(all_of(deg_genes_present), mean, na.rm = TRUE))

print(avg_expr)

library(dplyr)

# Your avg_expr tibble (already computed)
avg_expr <- tibble::tibble(
  cluster = factor(c(0,1,2,3,4)),
  ICOS    = c(0.490,0.910,0.449,0.480,0.149),
  CD2     = c(0.193,0.157,0.153,0.116,0.116),
  CTLA4   = c(0.0792,0,0.114,0.0578,0),
  PDCD1   = c(0.0396,0,0.0495,0.0578,0),
  CD28    = c(0.0396,0.120,0.0495,0.0289,0.0578),
  TNFRSF9 = c(0.0908,0.606,0.0248,0.202,0.231)
)

# 1) Compute thresholds (median per gene across clusters)
thresholds <- apply(avg_expr[,-1], 2, median)

# 2) Add High/Low calls as new columns (preserve numeric values)
expr_calls <- avg_expr
for (gene in names(thresholds)) {
  expr_calls[[paste0(gene,"_call")]] <- ifelse(avg_expr[[gene]] > thresholds[gene], "High", "Low")
}

# 3) Annotate clusters using a scoring system
cluster_annotations <- sapply(avg_expr$cluster, function(cl){
  row <- expr_calls[expr_calls$cluster==cl,]
  
  cd4_score <- sum(row[,c("ICOS_call","PDCD1_call","CTLA4_call")]=="High")
  cd8_score <- sum(row[,c("CD2_call","CD28_call")]=="High")
  nk_score  <- sum(row[,c("TNFRSF9_call")]=="High")
  
  if(cd4_score >= max(cd8_score, nk_score)) {
    return("CD4-like")
  } else if(cd8_score >= max(cd4_score, nk_score)) {
    return("CD8-like")
  } else if(nk_score >= max(cd4_score, cd8_score)) {
    return("NK-like")
  } else {
    return("Unassigned")
  }
})

expr_calls$Annotation <- cluster_annotations

# 4) Print final table
print(expr_calls)





#2b subpopulation cluster 13
library(Seurat)
library(dplyr)

# 1) Subset cluster 13
cl13 <- subset(mfs1_clean, idents = "13")

# 2) NK-oriented subset from your list
nk_genes <- unique(c(
  "CD96", "CST7", "TNFRSF9", "CLEC2D", "HCST", "KLRB1",
  "SYTL3", "FYB1", "SRGN", "FYN"
))

# 3) Filter to genes present in cl13
nk_genes_present <- nk_genes[nk_genes %in% rownames(cl13)]
# Check variance
gene_vars <- apply(FetchData(cl13, vars = nk_genes_present), 2, var)
nk_genes_filtered <- names(gene_vars[gene_vars > 0])

# Rerun PCA with filtered genes
cl13 <- ScaleData(cl13, features = nk_genes_filtered)
cl13 <- RunPCA(cl13, features = nk_genes_filtered)
cl13 <- FindNeighbors(cl13, dims = 1:3)
cl13 <- FindClusters(cl13, resolution = 0.6)
cl13 <- RunUMAP(cl13, dims = 1:3)

DimPlot(cl13, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# 5) Generate average expression table per cluster
expr_mat <- FetchData(cl13, vars = nk_genes_present)

avg_expr <- expr_mat %>%
  mutate(cluster = cl13$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(avg_expr)

#CD8/CD4/NK cluster 13
cd4_genes <- c(
  "IL7R","ICOS","CTLA4","PDCD1","TIGIT","TOX","BCL11B","GATA3",
  "CD6","CD69","CXCR4","IKZF2","BATF","STAT4","CD48","IL21R",
  "CCR7","LTB","PRKCQ","IL2RG","PTPN22","RUNX3","IKZF1","ARHGDIB",
  "SAMSN1","TNFRSF18","TNFRSF4","NFATC2"
)

cd8_genes <- c(
  "CD2","CD3D","CD3E","CD3G","CD247","CD28","LCK","TRBC1","TRBC2",
  "RUNX3","SYTL3","ITK","SKAP1","FYB1","PRKCQ","STAT5B","STAT5A",
  "IKZF1","IKZF2","PTPN22","DOCK2","RAC2","GZMB","PRF1","TAGAP",
  "TNFRSF9","TNFRSF18"
)

nk_genes <- c(
  "CD96","CST7","TNFRSF9","CLEC2D","HCST","KLRB1","SYTL3","FYB1",
  "SRGN","FYN","DOCK2","RAC2","GBP5","GNG2","TAGAP","TNFRSF4",
  "ISG20","CXCR4","ARHGAP15","STAT4","STAT5B","STAT5A","IKZF2",
  "NFATC2","ADGRE5","SERPINB9"
)


library(Seurat)
library(dplyr)

# Subset cluster 13
cl13 <- subset(mfs1_clean, idents = "13")

# Filter to genes present
cd4_genes_present <- cd4_genes[cd4_genes %in% rownames(cl13)]

# Remove zero‑variance genes
gene_vars_cd4 <- apply(FetchData(cl13, vars = cd4_genes_present), 2, var)
cd4_genes_filtered <- names(gene_vars_cd4[gene_vars_cd4 > 0])

# PCA/UMAP clustering
cl13_cd4 <- ScaleData(cl13, features = cd4_genes_filtered)
cl13_cd4 <- RunPCA(cl13_cd4, features = cd4_genes_filtered)
dims_use <- 1:min(10, ncol(cl13_cd4[["pca"]]@cell.embeddings))
cl13_cd4 <- FindNeighbors(cl13_cd4, dims = dims_use)
cl13_cd4 <- FindClusters(cl13_cd4, resolution = 0.6)
cl13_cd4 <- RunUMAP(cl13_cd4, dims = dims_use)

DimPlot(cl13_cd4, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# Average expression table
expr_mat_cd4 <- FetchData(cl13_cd4, vars = cd4_genes_present)
avg_expr_cd4 <- expr_mat_cd4 %>%
  mutate(cluster = cl13_cd4$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(avg_expr_cd4)
cl13_cd4@commands$NormalizeData

library(Seurat)
library(dplyr)

# --- Step 1: Define CD4 gene panel ---
cd4_genes <- c(
  "IL7R","ICOS","CTLA4","PDCD1","TIGIT","TOX","BCL11B","GATA3",
  "CD6","CD69","CXCR4","IKZF2","BATF","STAT4","CD48","IL21R",
  "CCR7","LTB","PRKCQ","IL2RG","PTPN22","RUNX3","IKZF1","ARHGDIB",
  "SAMSN1","TNFRSF18","TNFRSF4","NFATC2"
)

# --- Step 2: Subset cluster 13 ---
cl13 <- subset(mfs1_clean, idents = "13")

# --- Step 3: Filter to genes present ---
cd4_genes_present <- cd4_genes[cd4_genes %in% rownames(cl13)]

# --- Step 4: Keep only cells with >2 CD4 markers expressed ---
expr_mat_cd4 <- FetchData(cl13, vars = cd4_genes_present, slot = "data")
cells_keep <- rownames(expr_mat_cd4)[rowSums(expr_mat_cd4 > 0) > 2]

cl13_cd4 <- subset(cl13, cells = cells_keep)

# --- Step 5: Remove zero‑variance genes ---
gene_vars_cd4 <- apply(FetchData(cl13_cd4, vars = cd4_genes_present, slot = "data"), 2, var)
cd4_genes_filtered <- names(gene_vars_cd4[gene_vars_cd4 > 0])

# --- Step 6: Scale, PCA, clustering, UMAP ---
cl13_cd4 <- ScaleData(cl13_cd4, features = cd4_genes_filtered)
cl13_cd4 <- RunPCA(cl13_cd4, features = cd4_genes_filtered)
dims_use <- 1:min(10, ncol(cl13_cd4[["pca"]]@cell.embeddings))
cl13_cd4 <- FindNeighbors(cl13_cd4, dims = dims_use)
cl13_cd4 <- FindClusters(cl13_cd4, resolution = 0.8)
cl13_cd4 <- RunUMAP(cl13_cd4, dims = dims_use, n.neighbors = min(15, ncol(cl13_cd4) - 1))

DimPlot(cl13_cd4, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 7: Average expression table (ln(1+UMI)) ---
expr_mat_cd4 <- FetchData(cl13_cd4, vars = cd4_genes_present, slot = "data")
avg_expr_cd4 <- expr_mat_cd4 %>%
  mutate(cluster = cl13_cd4$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(avg_expr_cd4)

# --- Step 8: Check how many cells were retained ---
cat("Original cells in cluster 13:", ncol(cl13), "\n")
cat("Filtered cells (>2 CD4 markers expressed):", ncol(cl13_cd4), "\n")

library(pheatmap)

# --- Step 9: Convert avg_expr_cd4 to matrix ---
avg_expr_mat_cd4 <- as.data.frame(avg_expr_cd4)
rownames(avg_expr_mat_cd4) <- avg_expr_mat_cd4$cluster
avg_expr_mat_cd4 <- avg_expr_mat_cd4[ , -1]  # remove cluster column

# --- Step 10: Heatmap by Z-score (row scaling) ---
pheatmap(avg_expr_mat_cd4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",   # <-- Z-score per gene across clusters
         color = colorRampPalette(c("navy","white","firebrick"))(50))





library(Seurat)
library(dplyr)

# Subset cluster 13
cl13 <- subset(mfs1_clean, idents = "13")

# --- CD8 markers ---
cd8_genes <- c(
  "CD2","CD3D","CD3E","CD3G","CD247","CD28","LCK","TRBC1","TRBC2","RUNX3",
  "SYTL3","ITK","SKAP1","FYB1","PRKCQ","STAT5B","STAT5A","IKZF1","IKZF2",
  "PTPN22","DOCK2","RAC2","TAGAP","TNFRSF9","TNFRSF18"
)

cd8_genes_present <- cd8_genes[cd8_genes %in% rownames(cl13)]
gene_vars_cd8 <- apply(FetchData(cl13, vars = cd8_genes_present), 2, var)
cd8_genes_filtered <- names(gene_vars_cd8[gene_vars_cd8 > 0])

cl13_cd8 <- ScaleData(cl13, features = cd8_genes_filtered)
cl13_cd8 <- RunPCA(cl13_cd8, features = cd8_genes_filtered)
dims_use_cd8 <- 1:min(10, ncol(cl13_cd8[["pca"]]@cell.embeddings))
cl13_cd8 <- FindNeighbors(cl13_cd8, dims = dims_use_cd8)
cl13_cd8 <- FindClusters(cl13_cd8, resolution = 0.6)
cl13_cd8 <- RunUMAP(cl13_cd8, dims = dims_use_cd8)

DimPlot(cl13_cd8, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

expr_mat_cd8 <- FetchData(cl13_cd8, vars = cd8_genes_present)
avg_expr_cd8 <- expr_mat_cd8 %>%
  mutate(cluster = cl13_cd8$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
print(avg_expr_cd8)

library(Seurat)
library(dplyr)

# --- Step 1: Define CD8 gene panel ---
cd8_genes <- c(
  "CD2","CD3D","CD3E","CD3G","CD247","CD28","LCK","TRBC1","TRBC2","RUNX3",
  "SYTL3","ITK","SKAP1","FYB1","PRKCQ","STAT5B","STAT5A","IKZF1","IKZF2",
  "PTPN22","DOCK2","RAC2","TAGAP","TNFRSF9","TNFRSF18"
)

cd8_genes_present <- cd8_genes[cd8_genes %in% rownames(cl13)]

# --- Step 2: Keep only cells with >2 CD8 genes expressed ---
expr_mat_cd8 <- FetchData(cl13, vars = cd8_genes_present, slot = "data")

# Filter cells: must express more than 2 markers (non-zero ln(1+UMI))
cells_keep <- rownames(expr_mat_cd8)[rowSums(expr_mat_cd8 > 0) > 2]

cl13_cd8 <- subset(cl13, cells = cells_keep)

# --- Step 3: Filter CD8 genes with non-zero variance ---
gene_vars_cd8 <- apply(FetchData(cl13_cd8, vars = cd8_genes_present, slot = "data"), 2, var)
cd8_genes_filtered <- names(gene_vars_cd8[gene_vars_cd8 > 0])

# --- Step 4: Scale, PCA, clustering, UMAP ---
cl13_cd8 <- ScaleData(cl13_cd8, features = cd8_genes_filtered)
cl13_cd8 <- RunPCA(cl13_cd8, features = cd8_genes_filtered)
dims_use_cd8 <- 1:min(10, ncol(cl13_cd8[["pca"]]@cell.embeddings))
cl13_cd8 <- FindNeighbors(cl13_cd8, dims = dims_use_cd8)
cl13_cd8 <- FindClusters(cl13_cd8, resolution = 1.0)
cl13_cd8 <- RunUMAP(cl13_cd8, dims = dims_use_cd8, n.neighbors = 10)

DimPlot(cl13_cd8, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 5: Average expression per subcluster (ln(1+UMI)) ---
expr_mat_cd8 <- FetchData(cl13_cd8, vars = cd8_genes_present, slot = "data")
avg_expr_cd8 <- expr_mat_cd8 %>%
  mutate(cluster = cl13_cd8$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(avg_expr_cd8)

# Extract mean CD28 expression for each cluster
cd28_means <- avg_expr_cd8 %>%
  select(cluster, CD28)

print(cd28_means)

# If you only want clusters 0 and 1:
cd28_means %>% filter(cluster %in% c("0","1"))

# --- Step 6: Check how many cells were retained ---
cat("Original cells:", ncol(cl13), "\n")
cat("Filtered cells (>2 CD8 markers expressed):", ncol(cl13_cd8), "\n")

library(pheatmap)

# --- Step 7: Convert avg_expr_cd8 to matrix ---
avg_expr_mat_cd8 <- as.data.frame(avg_expr_cd8)
rownames(avg_expr_mat_cd8) <- avg_expr_mat_cd8$cluster
avg_expr_mat_cd8 <- avg_expr_mat_cd8[ , -1]  # remove cluster column

# --- Step 8: Heatmap by Z-score (row scaling) ---
pheatmap(avg_expr_mat_cd8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",   # <-- Z-score per gene across clusters
         color = colorRampPalette(c("navy","white","firebrick"))(50))


# --- NK markers ---
nk_genes <- c(
  "CD96","CST7","TNFRSF9","CLEC2D","HCST","KLRB1","SYTL3","FYB1","SRGN","FYN",
  "DOCK2","RAC2","GBP5","GNG2","TAGAP","TNFRSF4","ISG20","CXCR4","ARHGAP15",
  "STAT4","STAT5B","STAT5A","IKZF2","NFATC2","ADGRE5","SERPINB9"
)

nk_genes_present <- nk_genes[nk_genes %in% rownames(cl13)]
gene_vars_nk <- apply(FetchData(cl13, vars = nk_genes_present), 2, var)
nk_genes_filtered <- names(gene_vars_nk[gene_vars_nk > 0])

cl13_nk <- ScaleData(cl13, features = nk_genes_filtered)
cl13_nk <- RunPCA(cl13_nk, features = nk_genes_filtered)
dims_use_nk <- 1:min(10, ncol(cl13_nk[["pca"]]@cell.embeddings))
cl13_nk <- FindNeighbors(cl13_nk, dims = dims_use_nk)
cl13_nk <- FindClusters(cl13_nk, resolution = 0.6)
cl13_nk <- RunUMAP(cl13_nk, dims = dims_use_nk)

DimPlot(cl13_nk, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

expr_mat_nk <- FetchData(cl13_nk, vars = nk_genes_present)
avg_expr_nk <- expr_mat_nk %>%
  mutate(cluster = cl13_nk$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
print(avg_expr_nk)



# Cells you intended to keep
length(cells_keep)

# Cells actually present
ncol(cl13_nk)

# Compare overlap
setdiff(colnames(cl13_nk), cells_keep)


library(Seurat)
library(dplyr)

# --- Step 1: Define NK gene panel ---
library(Seurat)
library(dplyr)

# --- Step 1: Define NK gene panel ---
nk_genes <- c(
  "CD96","CST7","TNFRSF9","CLEC2D","HCST","KLRB1","SYTL3","FYB1","SRGN","FYN",
  "DOCK2","RAC2","GBP5","GNG2","TAGAP","TNFRSF4","ISG20","CXCR4","ARHGAP15",
  "STAT4","STAT5B","STAT5A","IKZF2","NFATC2","ADGRE5","SERPINB9"
)

nk_genes_present <- nk_genes[nk_genes %in% rownames(cl13)]

# --- Step 2: Keep only cells with ≥N NK genes expressed ---
expr_mat_nk <- FetchData(cl13, vars = nk_genes_present, layer = "data")

# Set threshold here (e.g., 2, 3, or 4 genes expressed)
threshold <- 3
cells_keep <- rownames(expr_mat_nk)[rowSums(expr_mat_nk > 0) >= threshold]

# Force exact subsetting by cell names
cl13_nk <- subset(cl13, cells = cells_keep)

# Confirm cell count
ncol(cl13_nk)



# --- Step 3: Filter NK genes with non-zero variance ---
gene_vars_nk <- apply(FetchData(cl13_nk, vars = nk_genes_present, slot = "data"), 2, var)
nk_genes_filtered <- names(gene_vars_nk[gene_vars_nk > 0])

# --- Step 4: Scale, PCA, clustering, UMAP ---
cl13_nk <- ScaleData(cl13_nk, features = nk_genes_filtered)
cl13_nk <- RunPCA(cl13_nk, features = nk_genes_filtered)
dims_use_nk <- 1:min(10, ncol(cl13_nk[["pca"]]@cell.embeddings))
cl13_nk <- FindNeighbors(cl13_nk, dims = dims_use_nk)
cl13_nk <- FindClusters(cl13_nk, resolution = 1.0)
cl13_nk <- RunUMAP(cl13_nk, dims = dims_use_nk, n.neighbors = 10)

DimPlot(cl13_nk, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 5: Average expression per subcluster (ln(1+UMI)) ---

expr_mat_nk <- FetchData(cl13_nk, vars = nk_genes_present, slot = "data")

avg_expr_nk <- expr_mat_nk %>%
  mutate(cluster = droplevels(cl13_nk$seurat_clusters)) %>%  # drop empty levels
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(avg_expr_nk)

library(pheatmap)

# --- Step 6: Convert avg_expr_nk to matrix ---
avg_expr_mat_nk <- as.data.frame(avg_expr_nk)
rownames(avg_expr_mat_nk) <- avg_expr_mat_nk$cluster
avg_expr_mat_nk <- avg_expr_mat_nk[ , -1]  # remove cluster column

# --- Step 7: Heatmap by Z-score (row scaling) ---
pheatmap(avg_expr_mat_nk,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",   # <-- Z-score per gene across clusters
         color = colorRampPalette(c("navy","white","firebrick"))(50))

# Number of cells kept
length(cells_keep)

# Compare before vs after filtering
ncol(cl13)       # total cells in original object
length(cells_keep)  # cells retained with ≥ threshold NK genes
ncol(cl13_nk)    # should match length(cells_keep)



#Cluster 15 DC/B cells
bcell_genes <- c(
  "MS4A1","CD79A","CD37","CD40","CD72","BANK1","SPIB","BLK","BCL11A",
  "IGHM","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DPB1","HLA-DMB",
  "HLA-DPA1","CD74","CD83","CD55","CD48","LTB","IKZF3","CXCR4","CCR7"
)

dc_genes <- c(
  "HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPB1","HLA-DMB","HLA-DPA1",
  "CD74","CD83","CD40","CCR7","SPIB","IRF8","BATF","CD48","CD69","CD52",
  "CD37","LAPTM5","PTPRC","CYTIP","CARD11","RUNX3","MEF2C","REL","STAT4"
)

library(Seurat)
library(dplyr)

# 1) Subset cluster 15
cl15 <- subset(mfs1_clean, idents = "15")

# --- B cell pipeline ---
bcell_genes_present <- bcell_genes[bcell_genes %in% rownames(cl15)]
gene_vars_b <- apply(FetchData(cl15, vars = bcell_genes_present), 2, var)
bcell_genes_filtered <- names(gene_vars_b[gene_vars_b > 0])

cl15_b <- ScaleData(cl15, features = bcell_genes_filtered)
cl15_b <- RunPCA(cl15_b, features = bcell_genes_filtered)
dims_use_b <- 1:min(10, ncol(cl15_b[["pca"]]@cell.embeddings))
cl15_b <- FindNeighbors(cl15_b, dims = dims_use_b)
cl15_b <- FindClusters(cl15_b, resolution = 0.6)
cl15_b <- RunUMAP(cl15_b, dims = dims_use_b)

DimPlot(cl15_b, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

expr_mat_b <- FetchData(cl15_b, vars = bcell_genes_present)
avg_expr_b <- expr_mat_b %>%
  mutate(cluster = cl15_b$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
print(avg_expr_b)


# --- DC pipeline ---
dc_genes_present <- dc_genes[dc_genes %in% rownames(cl15)]
gene_vars_dc <- apply(FetchData(cl15, vars = dc_genes_present), 2, var)
dc_genes_filtered <- names(gene_vars_dc[gene_vars_dc > 0])

cl15_dc <- ScaleData(cl15, features = dc_genes_filtered)
cl15_dc <- RunPCA(cl15_dc, features = dc_genes_filtered)
dims_use_dc <- 1:min(10, ncol(cl15_dc[["pca"]]@cell.embeddings))
cl15_dc <- FindNeighbors(cl15_dc, dims = dims_use_dc)
cl15_dc <- FindClusters(cl15_dc, resolution = 0.6)
cl15_dc <- RunUMAP(cl15_dc, dims = dims_use_dc)

DimPlot(cl15_dc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

expr_mat_dc <- FetchData(cl15_dc, vars = dc_genes_present)
avg_expr_dc <- expr_mat_dc %>%
  mutate(cluster = cl15_dc$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
print(avg_expr_dc)

# Percent expressing per subcluster for key DC markers
dc_core <- c("HLA-DRB1","HLA-DQA1","CCR7","CD83","IRF8")
df <- FetchData(cl15_dc, vars = dc_core)
df$cluster <- cl15_dc$seurat_clusters

pct_expr <- df %>%
  group_by(cluster) %>%
  summarise(across(all_of(dc_core), ~ mean(. > 0), .names = "pct_{col}"))

print(pct_expr)

##No clear B‑cell cluster: The canonical B‑cell markers (MS4A1, CD79A, IGHM, BLK, BANK1) are either low or absent.

#Antigen‑presentation genes dominate: HLA‑class II and CD74/CD83 are higher, which suggests DC‑like or activated APCs, not pure B cells.

#Cluster 15 may be mixed: It could represent DCs, monocytes, or activated APCs with some low B‑cell contamination, but not a distinct B‑cell population.


#cluster 16 B cells
library(Seurat)
library(dplyr)

# 1) Subset cluster 16
cl16 <- subset(mfs1_clean, idents = "16")

# 2) Define B cell marker set (expanded from your list)
bcell_genes <- c(
  "MZB1","JCHAIN","DERL3","POU2AF1","CD79A","CD79B","FCRL5","TNFRSF17",
  "IGHM","IGHG1","IGHG3","IGHG4","IGHA1","IGKC","IGLC2","IGLC3",
  "BLNK","SLAMF7","CD27","SDC1","CD38","PRDM1","IRF4","XBP1",
  "CD37","CXCR4","ST6GAL1","PLCG2","DAPP1","LAPTM5","BANK1"
)

# 3) Filter to genes present in cluster 16
bcell_genes_present <- bcell_genes[bcell_genes %in% rownames(cl16)]

# 4) Remove zero‑variance genes
gene_vars_b <- apply(FetchData(cl16, vars = bcell_genes_present), 2, var)
bcell_genes_filtered <- names(gene_vars_b[gene_vars_b > 0])

# 5) PCA/UMAP clustering
cl16_b <- ScaleData(cl16, features = bcell_genes_filtered)
cl16_b <- RunPCA(cl16_b, features = bcell_genes_filtered)
dims_use_b <- 1:min(10, ncol(cl16_b[["pca"]]@cell.embeddings))
cl16_b <- FindNeighbors(cl16_b, dims = dims_use_b)
cl16_b <- FindClusters(cl16_b, resolution = 0.6)
cl16_b <- RunUMAP(cl16_b, dims = dims_use_b)

DimPlot(cl16_b, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# 6) Average expression table per cluster
expr_mat_b <- FetchData(cl16_b, vars = bcell_genes_present)
avg_expr_b <- expr_mat_b %>%
  mutate(cluster = cl16_b$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(avg_expr_b)


##Cluster 16 is a small, homogeneous B‑cell group. The PCA loadings suggest a mix of plasma‑like and naïve/memory‑like signatures, 
#but with only 45 cells, it’s not enough to form distinct subclusters.

#Cluster 5
library(Seurat)
library(dplyr)
library(pheatmap)

macrophage_genes <- c(
  # Core macrophage markers
  "CD68","CD163","MRC1","MSR1","CSF1R","AIF1","FOLR2","VSIG4","STAB1",
  "HLA-DRA","HLA-DRB1","HLA-DPB1","C1QA","C1QB","C1QC",
  # M1-like (pro-inflammatory)
  "IL1B","TNF","TLR2","TLR4","NLRP3","CXCL8","CCR1","CD86",
  # M2-like (anti-inflammatory/tissue repair)
  "IL10RA","HMOX1","MERTK","TGFBI","CD163L1","MSR1","VSIG4"
)
# Subset cluster 5
c5 <- subset(mfs1_clean, idents = "5")

# Keep only genes present
macrophage_genes_present <- macrophage_genes[macrophage_genes %in% rownames(c5)]

# Filter cells with >2 macrophage markers expressed
expr_mat_mac <- FetchData(c5, vars = macrophage_genes_present, layer = "data")
cells_keep <- rownames(expr_mat_mac)[rowSums(expr_mat_mac > 0) > 2]

c5_mac <- subset(c5, cells = cells_keep)

cat("Original cells in cluster 5:", ncol(c5), "\n")
cat("Filtered macrophage-like cells:", ncol(c5_mac), "\n")



gene_vars_mac <- apply(FetchData(c5_mac, vars = macrophage_genes_present, slot = "data"), 2, var)
macrophage_genes_filtered <- names(gene_vars_mac[gene_vars_mac > 0])

c5_mac <- ScaleData(c5_mac, features = macrophage_genes_filtered)
c5_mac <- RunPCA(c5_mac, features = macrophage_genes_filtered)
dims_use_mac <- 1:min(10, ncol(c5_mac[["pca"]]@cell.embeddings))
c5_mac <- FindNeighbors(c5_mac, dims = dims_use_mac)
c5_mac <- FindClusters(c5_mac, resolution = 0.6)
c5_mac <- RunUMAP(c5_mac, dims = dims_use_mac, n.neighbors = min(15, ncol(c5_mac) - 1))

DimPlot(c5_mac, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

expr_mat_mac <- FetchData(c5_mac, vars = macrophage_genes_present, slot = "data")
avg_expr_mac <- expr_mat_mac %>%
  mutate(cluster = droplevels(c5_mac$seurat_clusters)) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

avg_expr_mat_mac <- as.data.frame(avg_expr_mac)
rownames(avg_expr_mat_mac) <- avg_expr_mat_mac$cluster
avg_expr_mat_mac <- avg_expr_mat_mac[ , -1]

# Lineage annotation
m1_markers <- c("IL1B","TNF","TLR2","TLR4","NLRP3","CXCL8","CCR1","CD86")
m2_markers <- c("IL10RA","HMOX1","MERTK","TGFBI","CD163L1","MSR1","VSIG4")

marker_lineage <- c(
  setNames(rep("M1_like", length(m1_markers)), m1_markers),
  setNames(rep("M2_like", length(m2_markers)), m2_markers),
  setNames(rep("Macrophage_core", length(macrophage_genes)), macrophage_genes)
)

genes_present <- intersect(names(marker_lineage), colnames(avg_expr_mat_mac))
annotation_row <- data.frame(Lineage = marker_lineage[genes_present])
rownames(annotation_row) <- genes_present

pheatmap::pheatmap(
  t(avg_expr_mat_mac[, genes_present]),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_row = annotation_row,
  color = colorRampPalette(c("navy","white","firebrick"))(50)
)





#Subpopulation mfs2
#C5
library(Seurat)
library(dplyr)
library(pheatmap)
library(org.Hs.eg.db)

# --- Step 0: Subset cluster 5 ---
c5 <- subset(mfs2_clean, idents = "5")
cat("Total cells in C5 (cluster 5):", ncol(c5), "\n")

# --- Step 1: Map monocyte markers (symbols) to Ensembl IDs ---
mono_symbols <- c(
  "CD14","LYZ","CST3","FCGR3A","CD68","SERPINA1","C1QA","C1QB","C1QC",
  "IL1B","TNFAIP2","CD86","PLAUR","CXCL8","TFEC","CSF1R","CD302",
  "FCGR1A","ITGAX","MRC1","NRP2","SIGLEC1","LILRB1","LILRB2","CXCL16","MAFB"
)

map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = mono_symbols,
                             keytype = "SYMBOL",
                             columns = "ENSEMBL")

mono_genes_present <- map$ENSEMBL[map$ENSEMBL %in% rownames(c5)]
cat("Number of monocyte markers found in dataset:", length(mono_genes_present), "\n")

# --- Step 2: Filter monocyte-like cells (>7 markers expressed) ---
expr_mat_mono <- FetchData(c5, vars = mono_genes_present, slot = "counts")
cells_keep <- rownames(expr_mat_mono)[rowSums(expr_mat_mono > 1) > 20]
c5_mono <- subset(c5, cells = cells_keep)
cat("Cells retained in C5 after monocyte marker filter (>7 markers):", ncol(c5_mono), "\n")

# --- Step 3: Reclustering ---
c5_mono <- ScaleData(c5_mono, features = mono_genes_present)
c5_mono <- RunPCA(c5_mono, features = mono_genes_present)
c5_mono <- FindNeighbors(c5_mono, dims = 1:10)
c5_mono <- FindClusters(c5_mono, resolution = 1.0)
c5_mono <- RunUMAP(c5_mono, dims = 1:10)

DimPlot(c5_mono, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 4: Average expression per subcluster ---
expr_mat_mono <- FetchData(c5_mono, vars = mono_genes_present, slot = "data")
avg_expr_mono <- expr_mat_mono %>%
  mutate(cluster = droplevels(c5_mono$seurat_clusters)) %>%
  group_by(cluster) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

avg_expr_mat <- as.data.frame(avg_expr_mono)
rownames(avg_expr_mat) <- avg_expr_mat$cluster
avg_expr_mat <- avg_expr_mat[, -1, drop = FALSE]

# --- Step 5: Map Ensembl IDs back to gene symbols for readability ---
map_back <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = colnames(avg_expr_mat),
                                  keytype = "ENSEMBL",
                                  columns = "SYMBOL")
id2symbol <- setNames(map_back$SYMBOL, map_back$ENSEMBL)

mapped_cols <- id2symbol[colnames(avg_expr_mat)]
keep <- nzchar(mapped_cols) & !is.na(mapped_cols)
avg_expr_mat <- avg_expr_mat[, keep, drop = FALSE]
colnames(avg_expr_mat) <- mapped_cols[keep]

# --- Step 6: Build lineage annotation ---
classical_mono    <- c("CD14","LYZ","CST3","FCGR3A","CD68","SERPINA1","C1QA","C1QB","C1QC")
intermediate_mono <- c("IL1B","TNFAIP2","CD86","PLAUR","CXCL8","TFEC","CSF1R","CD302")
nonclassical_mono <- c("FCGR1A","ITGAX","MRC1","NRP2","SIGLEC1","LILRB1","LILRB2","CXCL16","MAFB")

marker_lineage <- c(
  setNames(rep("Classical", length(classical_mono)), classical_mono),
  setNames(rep("Intermediate", length(intermediate_mono)), intermediate_mono),
  setNames(rep("NonClassical", length(nonclassical_mono)), nonclassical_mono)
)

genes_present <- intersect(names(marker_lineage), colnames(avg_expr_mat))
annotation_row <- data.frame(Lineage = marker_lineage[genes_present])
rownames(annotation_row) <- genes_present

# --- Step 7: Heatmap ---
mat <- avg_expr_mat[, genes_present, drop = FALSE]
vars <- apply(mat, 2, function(x) var(x, na.rm = TRUE)); vars[is.na(vars)] <- 0
scale_opt <- if (all(vars == 0)) "none" else "row"

pheatmap(
  t(mat),
  cluster_rows = FALSE,
  cluster_cols = (nrow(mat) > 1),
  scale = scale_opt,
  annotation_row = annotation_row,
  color = colorRampPalette(c("navy","white","firebrick"))(50),
  main = "Monocyte subclusters (Classical / Intermediate / NonClassical)"
)




library(Seurat)
library(dplyr)
library(pheatmap)
library(biomaRt)


# ============================
# 1. Macrophages
# ============================


library(Seurat)
library(dplyr)
library(pheatmap)
library(biomaRt)

# --- Ensure the active identities are seurat_clusters before subsetting ---
# If Idents are not seurat_clusters, set them explicitly:
if (!identical(colnames(mfs2_clean@meta.data)[1], "seurat_clusters")) {
  # If seurat_clusters exists in metadata, set Idents to it
  if ("seurat_clusters" %in% colnames(mfs2_clean@meta.data)) {
    Idents(mfs2_clean) <- "seurat_clusters"
  } else {
    stop("seurat_clusters not found in metadata; set Idents(mfs2_clean) appropriately or subset with a logical filter.")
  }
}

# --- Step 0: Subset cluster 5 from mfs2_clean ---
c5 <- subset(mfs2_clean, idents = "5")

# --- Marker sets (single definition) ---
macrophage_core <- c("CD68","CD163","MRC1","MSR1","CSF1R","AIF1","HLA-DRA","HLA-DRB1","HLA-DPB1","C1QA","C1QB","C1QC")
m1_like <- c("IL1B","TNF","TLR2","TLR4","NLRP3","CXCL8","CXCL2","CD86")
m2_like <- c("IL10RA","HMOX1","MERTK","TGFBI","VSIG4","SIGLEC1","FOLR2","CCR1","STAB1")
macrophage_genes <- unique(c(macrophage_core, m1_like, m2_like))

# --- Step 1: Prepare Biomart once and convert symbols to dataset IDs if needed ---
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

convert_to_dataset_ids <- function(gene_list, seurat_obj, mart) {
  rn <- rownames(seurat_obj)
  if (length(rn) == 0) stop("Object has no feature names.")
  uses_ens <- grepl("^ENS", rn[1])
  if (!uses_ens) {
    return(gene_list[gene_list %in% rn])
  } else {
    map <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                 filters = "hgnc_symbol", values = gene_list, mart = mart)
    ids <- unique(map$ensembl_gene_id)
    ids[ids %in% rn]
  }
}

macrophage_genes_present <- convert_to_dataset_ids(macrophage_genes, c5, mart)
if (length(macrophage_genes_present) == 0) stop("No macrophage markers found in object rownames.")

# --- Step 2: Filter cells with >2 markers expressed (ln(1+UMI) > 0) ---
expr_mat_mac <- FetchData(c5, vars = macrophage_genes_present, slot = "data")
cells_keep <- rownames(expr_mat_mac)[rowSums(expr_mat_mac > 0) > 20]
c5_mac <- subset(c5, cells = cells_keep)

cat("Original cells in cluster 5:", ncol(c5), "\n")
cat("Macrophage-like cells retained (>2 markers):", ncol(c5_mac), "\n")

# --- Step 3: Filter features by non-zero variance within retained cells ---
gene_vars_mac <- apply(FetchData(c5_mac, vars = macrophage_genes_present, slot = "data"), 2, var)
macrophage_genes_filtered <- names(gene_vars_mac[gene_vars_mac > 0])
if (length(macrophage_genes_filtered) < 5) warning("Few variable macrophage markers after filtering.")

# --- Step 4: Scale, PCA, neighbors, clustering, UMAP ---
c5_mac <- ScaleData(c5_mac, features = macrophage_genes_filtered)
c5_mac <- RunPCA(c5_mac, features = macrophage_genes_filtered)
dims_use_mac <- 1:min(10, ncol(c5_mac[["pca"]]@cell.embeddings))
c5_mac <- FindNeighbors(c5_mac, dims = dims_use_mac)
c5_mac <- FindClusters(c5_mac, resolution = 1.0)
c5_mac <- RunUMAP(c5_mac, dims = dims_use_mac, n.neighbors = max(5, min(15, ncol(c5_mac) - 1)))

DimPlot(c5_mac, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 5: Average expression per macrophage subcluster (ln(1+UMI)) ---
expr_mat_mac <- FetchData(c5_mac, vars = macrophage_genes_present, slot = "data")
avg_expr_mac <- expr_mat_mac %>%
  mutate(cluster = droplevels(c5_mac$seurat_clusters)) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(avg_expr_mac)

# --- Step 6: Optional symbol conversion for readability ---
avg_expr_mac_sym <- avg_expr_mac
if (grepl("^ENS", colnames(avg_expr_mac_sym)[2])) {
  ensembl_ids <- colnames(avg_expr_mac_sym)[-1]
  mapping <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                   filters = "ensembl_gene_id", values = ensembl_ids, mart = mart)
  id2symbol <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
  colnames(avg_expr_mac_sym)[-1] <- ifelse(
    nzchar(id2symbol[colnames(avg_expr_mac_sym)[-1]]),
    id2symbol[colnames(avg_expr_mac_sym)[-1]],
    colnames(avg_expr_mac_sym)[-1]
  )
}

# --- Step 7: Z-score heatmap with core/M1/M2 annotation ---
avg_expr_mat <- as.data.frame(avg_expr_mac_sym)
rownames(avg_expr_mat) <- avg_expr_mat$cluster
avg_expr_mat <- avg_expr_mat[, -1, drop = FALSE]

ann_labels <- c(
  setNames(rep("Macrophage_core", length(macrophage_core)), macrophage_core),
  setNames(rep("M1_like", length(m1_like)), m1_like),
  setNames(rep("M2_like", length(m2_like)), m2_like)
)
ann_labels <- ann_labels[intersect(names(ann_labels), colnames(avg_expr_mat))]
annotation_row <- data.frame(Lineage = ann_labels)
rownames(annotation_row) <- names(ann_labels)

pheatmap(
  t(avg_expr_mat[, names(ann_labels), drop = FALSE]),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_row = annotation_row,
  color = colorRampPalette(c("navy","white","firebrick"))(50),
  main = "Macrophage subclusters (Z-score, core/M1/M2)"
)

# Distribution of marker counts per cell
marker_counts <- rowSums(expr_mat_mac > 0)
summary(marker_counts)
table(marker_counts)

# Compare to raw counts
expr_mat_mac_counts <- FetchData(c5, vars = macrophage_genes_present, slot = "counts")
marker_counts_raw <- rowSums(expr_mat_mac_counts > 0)
summary(marker_counts_raw)
table(marker_counts_raw)

sum(marker_counts_raw > 2)   # cells with >2 markers
sum(marker_counts_raw > 5)   # cells with >5 markers
sum(marker_counts_raw > 7)   # cells with >7 markers


# ============================
# 2. Dendritic Cells
# ============================
dc_genes <- c(
  "HLA-DRA","HLA-DRB1","HLA-DMB","HLA-DQA1","HLA-DQB1","HLA-DPB1",
  "CD74","CD83","IRF8","SPI1","CIITA","ITGAX","FCER1G","CLEC7A",
  "LILRB4","GPR183","TFEC","CXCL8","IL18","CD40","CD86","CD80",
  "CST3","IFI30","MPEG1","CYBB","TYROBP"
)

dc_present <- convert_to_dataset_ids(dc_genes, mfs2_c5)
if (length(dc_present) == 0) stop("No DC markers found in dataset")

gene_vars_dc <- apply(FetchData(mfs2_c5, vars = dc_present), 2, var)
dc_filtered <- names(gene_vars_dc[gene_vars_dc > 0])

mfs2_c5_dc <- ScaleData(mfs2_c5, features = dc_filtered)
mfs2_c5_dc <- RunPCA(mfs2_c5_dc, features = dc_filtered)
dims_use_dc <- 1:min(10, ncol(mfs2_c5_dc[["pca"]]@cell.embeddings))
mfs2_c5_dc <- FindNeighbors(mfs2_c5_dc, dims = dims_use_dc)
mfs2_c5_dc <- FindClusters(mfs2_c5_dc, resolution = 0.6)
mfs2_c5_dc <- RunUMAP(mfs2_c5_dc, dims = dims_use_dc)

DimPlot(mfs2_c5_dc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

expr_mat_dc <- FetchData(mfs2_c5_dc, vars = dc_present)
avg_expr_dc <- expr_mat_dc %>%
  mutate(cluster = mfs2_c5_dc$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

avg_expr_dc <- replace_ids_with_symbols(avg_expr_dc)
print(head(avg_expr_dc))

mat <- as.matrix(avg_expr_dc[,-1]); rownames(mat) <- avg_expr_dc$cluster

# --- Define DC subtype marker sets ---
cDC1 <- c("CLEC9A","XCR1","BATF3","IRF8")
cDC2 <- c("CD1C","FCER1A","IRF4","CLEC10A")
pDC  <- c("IL3RA","TCF4","LILRA4","GZMB")

intersect(cDC1, colnames(avg_expr_dc))
intersect(cDC2, colnames(avg_expr_dc))
intersect(pDC,  colnames(avg_expr_dc))


cDC2 %in% rownames(mfs2_c5)
pDC  %in% rownames(mfs2_c5)

all_markers <- colnames(avg_expr_dc)[-1]
marker_lineage <- c(
  setNames(rep("cDC1", length(cDC1)), cDC1),
  setNames(rep("cDC2", length(cDC2)), cDC2),
  setNames(rep("pDC", length(pDC)), pDC)
)
marker_lineage <- marker_lineage[all_markers]

annotation_row <- data.frame(Subtype = marker_lineage)
rownames(annotation_row) <- all_markers

pheatmap(t(mat),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_row = annotation_row,
         main = "Dendritic cell subclusters with subtype annotation")

#C0 CD8 
library(Seurat)
library(dplyr)
library(pheatmap)
library(biomaRt)

# --- Subset cluster 0 (CD8 cells) ---
mfs2_c0 <- subset(mfs2_clean, idents = "0")
head(rownames(mfs2_c0))

# --- Define CD8 marker set (select 20–30 from your long list) ---
cd8_markers <- c(
  "CD8A","CD3D","CD3E","CD3G","GZMK","GZMA","GZMM","GZMH",
  "NKG7","CCL5","CST7","CTSW","CXCR4","CD69","CRTAM","RUNX3",
  "LAG3","CXCR3","IL7R","CD2","CD44","CD48","CD52","TRAC",
  "TRBC1","TRBC2","PTPRC","HCST"
)


# --- Step 1: Convert to dataset IDs if needed ---
library(org.Hs.eg.db)

# Map symbols to Ensembl IDs
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = cd8_markers,
                             keytype = "SYMBOL",
                             columns = c("ENSEMBL","SYMBOL"))

# Keep only those present in your dataset
cd8_present <- map$ENSEMBL[map$ENSEMBL %in% rownames(mfs2_c0)]

if (length(cd8_present) == 0) stop("No CD8 markers found in dataset")

# --- Step 2: Select cells with >2 markers expressed ---
expr_mat_cd8_counts <- FetchData(mfs2_c0, vars = cd8_present, slot = "counts")
cells_keep <- rownames(expr_mat_cd8_counts)[rowSums(expr_mat_cd8_counts > 1) > 15]
mfs2_c0 <- subset(mfs2_c0, cells = cells_keep)

# --- Step 3: Filter genes by variance ---
gene_vars_cd8 <- apply(FetchData(mfs2_c0, vars = cd8_present), 2, var)
cd8_filtered <- names(gene_vars_cd8[gene_vars_cd8 > 0])

# --- Step 4: Reclustering ---
mfs2_c0_cd8 <- ScaleData(mfs2_c0, features = cd8_filtered)
mfs2_c0_cd8 <- RunPCA(mfs2_c0_cd8, features = cd8_filtered)
dims_use_cd8 <- 1:min(10, ncol(mfs2_c0_cd8[["pca"]]@cell.embeddings))
mfs2_c0_cd8 <- FindNeighbors(mfs2_c0_cd8, dims = dims_use_cd8)
mfs2_c0_cd8 <- FindClusters(mfs2_c0_cd8, resolution = 0.6)
mfs2_c0_cd8 <- RunUMAP(mfs2_c0_cd8, dims = dims_use_cd8)

DimPlot(mfs2_c0_cd8, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 5: Average expression per subcluster ---
expr_mat_cd8 <- FetchData(mfs2_c0_cd8, vars = cd8_present)
avg_expr_cd8 <- expr_mat_cd8 %>%
  mutate(cluster = mfs2_c0_cd8$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# --- Step 6: Map Ensembl IDs back to gene symbols ---
ensembl_ids <- colnames(avg_expr_cd8)[-1]
mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                 filters="ensembl_gene_id",
                 values=ensembl_ids,
                 mart=mart)
id2symbol <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
colnames(avg_expr_cd8)[-1] <- ifelse(nzchar(id2symbol[colnames(avg_expr_cd8)[-1]]),
                                     id2symbol[colnames(avg_expr_cd8)[-1]],
                                     colnames(avg_expr_cd8)[-1])

# --- Step 7: Annotate CD8 subtypes ---
Effector_CD8 <- c("GZMK","GZMA","GZMM","NKG7","CCL5","CST7","CTSW")
Memory_CD8   <- c("IL7R","CXCR3","CD44","CD69","CRTAM","RUNX3")
Naive_CD8    <- c("CCR7","LEF1","TCF7","SELL")
Exhausted_CD8<- c("LAG3","PDCD1","TIGIT","CTLA4")

all_markers <- colnames(avg_expr_cd8)[-1]
marker_lineage <- c(
  setNames(rep("Effector", length(Effector_CD8)), Effector_CD8),
  setNames(rep("Memory",   length(Memory_CD8)),   Memory_CD8),
  setNames(rep("Naive",    length(Naive_CD8)),    Naive_CD8),
  setNames(rep("Exhausted",length(Exhausted_CD8)),Exhausted_CD8)
)
marker_lineage <- marker_lineage[all_markers]

annotation_row <- data.frame(Subtype = marker_lineage)
rownames(annotation_row) <- all_markers

# --- Step 8: Heatmap ---
mat <- as.matrix(avg_expr_cd8[,-1]); rownames(mat) <- avg_expr_cd8$cluster
vars <- apply(mat, 2, function(x) var(x, na.rm = TRUE)); vars[is.na(vars)] <- 0
scale_opt <- if (all(vars == 0)) "none" else "row"

pheatmap(t(mat),
         cluster_rows = FALSE,
         cluster_cols = (nrow(mat) > 1),
         scale = scale_opt,
         annotation_row = annotation_row,
         main = "CD8 T cell subclusters with subtype annotation",
         cellheight = 12, cellwidth = 12, fontsize = 10)





library(Seurat)
library(dplyr)
library(pheatmap)
library(org.Hs.eg.db)   # stable offline mapping

# --- Step 0: Subset cluster 1 (candidate CD4 cells) ---
mfs2_c1 <- subset(mfs2_clean, idents = "1")

# --- Step 1: Define marker panel (20–30 from your long list) ---
cd4_markers <- c(
  "IL7R","CCR7","CD28","ICOS","TNFRSF4","TNFRSF18","CTLA4","TIGIT",
  "CD2","CD5","CD44","CD52","CD3D","CD3E","CD3G","CD27","CD247","IL2RA",
  "FOXP1","SATB1","BCL11B","GATA3","STAT3","BATF","MAF","TOX"
)

# --- Step 2: Map symbols to Ensembl IDs if needed ---
if (grepl("^ENSG", rownames(mfs2_c1)[1])) {
  map <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = cd4_markers,
                               keytype = "SYMBOL",
                               columns = "ENSEMBL")
  cd4_present <- map$ENSEMBL[map$ENSEMBL %in% rownames(mfs2_c1)]
} else {
  cd4_present <- cd4_markers[cd4_markers %in% rownames(mfs2_c1)]
}
if (length(cd4_present) == 0) stop("No CD4 markers found in dataset")

# --- Step 3: Select cells with >2 markers expressed ---
expr_mat_cd4_counts <- FetchData(mfs2_c1, vars = cd4_present, slot = "counts")
cells_keep <- rownames(expr_mat_cd4_counts)[rowSums(expr_mat_cd4_counts > 1) > 10]
mfs2_c1 <- subset(mfs2_c1, cells = cells_keep)
cat("Cells retained in cluster 1 after CD4 marker filter (>2 markers):", ncol(mfs2_c1), "\n")

# --- Step 4: Filter genes by variance ---
gene_vars_cd4 <- apply(FetchData(mfs2_c1, vars = cd4_present), 2, var)
cd4_filtered <- names(gene_vars_cd4[gene_vars_cd4 > 0])

# --- Step 5: Reclustering ---
mfs2_c1_cd4 <- ScaleData(mfs2_c1, features = cd4_filtered)
mfs2_c1_cd4 <- RunPCA(mfs2_c1_cd4, features = cd4_filtered)
dims_use_cd4 <- 1:min(10, ncol(mfs2_c1_cd4[["pca"]]@cell.embeddings))
mfs2_c1_cd4 <- FindNeighbors(mfs2_c1_cd4, dims = dims_use_cd4)
mfs2_c1_cd4 <- FindClusters(mfs2_c1_cd4, resolution = 1.2)
mfs2_c1_cd4 <- RunUMAP(mfs2_c1_cd4, dims = dims_use_cd4)

DimPlot(mfs2_c1_cd4, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 6: Average expression per subcluster ---
expr_mat_cd4 <- FetchData(mfs2_c1_cd4, vars = cd4_filtered)
avg_expr_cd4 <- expr_mat_cd4 %>%
  mutate(cluster = mfs2_c1_cd4$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# --- Step 7: Map Ensembl IDs back to symbols if needed ---
ensembl_ids <- colnames(avg_expr_cd4)[-1]
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = ensembl_ids,
                             keytype = "ENSEMBL",
                             columns = "SYMBOL")
id2symbol <- setNames(map$SYMBOL, map$ENSEMBL)
colnames(avg_expr_cd4)[-1] <- ifelse(nzchar(id2symbol[colnames(avg_expr_cd4)[-1]]),
                                     id2symbol[colnames(avg_expr_cd4)[-1]],
                                     colnames(avg_expr_cd4)[-1])

print(head(avg_expr_cd4))

# --- Step 8: Annotate CD4 subtypes ---
Naive_CD4     <- c("CCR7","IL7R","CD27")
Memory_CD4    <- c("CD44","CD28","FOXP1","SATB1")
Effector_CD4  <- c("ICOS","TNFRSF4","TNFRSF18","CD52","CD2","CD5")
Regulatory_CD4<- c("CTLA4","IL2RA","TIGIT","TOX","BCL11B","GATA3")

all_markers <- colnames(avg_expr_cd4)[-1]
marker_lineage <- c(
  setNames(rep("Naive", length(Naive_CD4)), Naive_CD4),
  setNames(rep("Memory", length(Memory_CD4)), Memory_CD4),
  setNames(rep("Effector", length(Effector_CD4)), Effector_CD4),
  setNames(rep("Regulatory", length(Regulatory_CD4)), Regulatory_CD4)
)
marker_lineage <- marker_lineage[all_markers]

annotation_row <- data.frame(Subtype = marker_lineage)
rownames(annotation_row) <- all_markers

# --- Step 9: Heatmap ---
mat <- as.matrix(avg_expr_cd4[,-1]); rownames(mat) <- avg_expr_cd4$cluster
vars <- apply(mat, 2, function(x) var(x, na.rm = TRUE)); vars[is.na(vars)] <- 0
scale_opt <- if (all(vars == 0)) "none" else "row"

pheatmap(t(mat),
         cluster_rows = FALSE,
         cluster_cols = (nrow(mat) > 1),
         scale = scale_opt,
         annotation_row = annotation_row,
         main = "CD4 T cell subclusters with subtype annotation",
         cellheight = 12, cellwidth = 12, fontsize = 10)



library(Seurat)
library(dplyr)
library(pheatmap)
library(org.Hs.eg.db)   # stable offline mapping

# --- Step 0: Subset cluster 2 (NK cells) ---
mfs2_c2 <- subset(mfs2_clean, idents = "2")

# --- Step 1: Define NK marker set (20–30 genes) ---
nk_markers <- c(
  "NKG7","GZMB","GZMA","GZMK","GZMH","PRF1","IFNG","CCL3","CCL4","CCL5",
  "XCL1","XCL2","CD69","CD8A","CD3D","CD3E","CD3G","CRTAM","CST7","CTSW",
  "KLRK1","SLAMF7","FASLG","RUNX3","ZEB2","TNFRSF9"
)

# --- Step 2: Map symbols to Ensembl IDs if needed ---
if (grepl("^ENSG", rownames(mfs2_c2)[1])) {
  map <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = nk_markers,
                               keytype = "SYMBOL",
                               columns = "ENSEMBL")
  nk_present <- map$ENSEMBL[map$ENSEMBL %in% rownames(mfs2_c2)]
} else {
  nk_present <- nk_markers[nk_markers %in% rownames(mfs2_c2)]
}
if (length(nk_present) == 0) stop("No NK markers found in dataset")

# --- Step 3: Select cells with >2 NK markers expressed ---
expr_mat_nk_counts <- FetchData(mfs2_c2, vars = nk_present, slot = "counts")
cells_keep <- rownames(expr_mat_nk_counts)[rowSums(expr_mat_nk_counts > 1) > 15]
mfs2_c2 <- subset(mfs2_c2, cells = cells_keep)
cat("Cells retained in cluster 2 after NK marker filter (>2 markers):", ncol(mfs2_c2), "\n")

# --- Step 4: Filter genes by variance ---
gene_vars_nk <- apply(FetchData(mfs2_c2, vars = nk_present), 2, var)
nk_filtered <- names(gene_vars_nk[gene_vars_nk > 0])

# --- Step 5: Reclustering ---
mfs2_c2_nk <- ScaleData(mfs2_c2, features = nk_filtered)
mfs2_c2_nk <- RunPCA(mfs2_c2_nk, features = nk_filtered)
dims_use_nk <- 1:min(10, ncol(mfs2_c2_nk[["pca"]]@cell.embeddings))
mfs2_c2_nk <- FindNeighbors(mfs2_c2_nk, dims = dims_use_nk)
mfs2_c2_nk <- FindClusters(mfs2_c2_nk, resolution = 1.0)
mfs2_c2_nk <- RunUMAP(mfs2_c2_nk, dims = dims_use_nk)

DimPlot(mfs2_c2_nk, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 6: Average expression per subcluster ---
expr_mat_nk <- FetchData(mfs2_c2_nk, vars = nk_filtered)
avg_expr_nk <- expr_mat_nk %>%
  mutate(cluster = mfs2_c2_nk$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# --- Step 7: Map Ensembl IDs back to gene symbols if needed ---
ensembl_ids <- colnames(avg_expr_nk)[-1]
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = ensembl_ids,
                             keytype = "ENSEMBL",
                             columns = "SYMBOL")
id2symbol <- setNames(map$SYMBOL, map$ENSEMBL)
colnames(avg_expr_nk)[-1] <- ifelse(nzchar(id2symbol[colnames(avg_expr_nk)[-1]]),
                                    id2symbol[colnames(avg_expr_nk)[-1]],
                                    colnames(avg_expr_nk)[-1])

print(head(avg_expr_nk))

# --- Step 8: Annotate NK subtypes ---
Cytotoxic_NK <- c("GZMB","GZMA","GZMK","GZMH","PRF1","NKG7","CST7","CTSW")
Activated_NK <- c("CD69","IFNG","TNFRSF9","FASLG","XCL1","XCL2","CCL3","CCL4","CCL5")
Memory_NK    <- c("RUNX3","ZEB2","SLAMF7","CRTAM")

all_markers <- colnames(avg_expr_nk)[-1]
marker_lineage <- c(
  setNames(rep("Cytotoxic", length(Cytotoxic_NK)), Cytotoxic_NK),
  setNames(rep("Activated", length(Activated_NK)), Activated_NK),
  setNames(rep("Memory", length(Memory_NK)), Memory_NK)
)
marker_lineage <- marker_lineage[all_markers]

annotation_row <- data.frame(Subtype = marker_lineage)
rownames(annotation_row) <- all_markers

# --- Step 9: Heatmap ---
mat <- as.matrix(avg_expr_nk[,-1]); rownames(mat) <- avg_expr_nk$cluster
vars <- apply(mat, 2, function(x) var(x, na.rm = TRUE)); vars[is.na(vars)] <- 0
scale_opt <- if (all(vars == 0)) "none" else "row"

pheatmap(t(mat),
         cluster_rows = FALSE,
         cluster_cols = (nrow(mat) > 1),
         scale = scale_opt,
         annotation_row = annotation_row,
         main = "NK cell subclusters with subtype annotation",
         cellheight = 12, cellwidth = 12, fontsize = 10)


#C4 NK cells 


library(Seurat)
library(dplyr)
library(pheatmap)
library(org.Hs.eg.db)   # stable offline mapping

# --- Step 0: Subset cluster 4 (NK cells) ---
mfs2_c4 <- subset(mfs2_clean, idents = "4")

# --- Step 1: Define NK marker panel (20–30 genes from your list) ---
nk_markers <- c(
  "CD69","TNFRSF9","CRTAM","CLEC2D","RUNX3","ZEB2","STAT4","REL","PTPN22",
  "DOCK2","ARHGAP25","CD247","CD96","KLRK1","IKZF1","IKZF3","PRKCQ","ITK",
  "VAV1","SYTL3","SKAP1","CYTIP","TRAF1","NR3C1","TOX"
)


# --- Step 2: Map symbols to Ensembl IDs if needed ---
if (grepl("^ENSG", rownames(mfs2_c4)[1])) {
  map <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = nk_markers,
                               keytype = "SYMBOL",
                               columns = "ENSEMBL")
  nk_present <- map$ENSEMBL[map$ENSEMBL %in% rownames(mfs2_c4)]
} else {
  nk_present <- nk_markers[nk_markers %in% rownames(mfs2_c4)]
}
if (length(nk_present) == 0) stop("No NK markers found in dataset")

# --- Step 3: Select cells with >2 NK markers expressed ---
expr_mat_nk_counts <- FetchData(mfs2_c4, vars = nk_present, slot = "counts")
cells_keep <- rownames(expr_mat_nk_counts)[rowSums(expr_mat_nk_counts > 1) > 5]
mfs2_c4 <- subset(mfs2_c4, cells = cells_keep)
cat("Cells retained in cluster 2 after NK marker filter (>2 markers):", ncol(mfs2_c2), "\n")

# --- Step 4: Filter genes by variance ---
gene_vars_nk <- apply(FetchData(mfs2_c4, vars = nk_present), 2, var)
nk_filtered <- names(gene_vars_nk[gene_vars_nk > 0])

# --- Step 5: Reclustering ---
mfs2_c4_nk <- ScaleData(mfs2_c4, features = nk_filtered)
mfs2_c4_nk <- RunPCA(mfs2_c4_nk, features = nk_filtered)
dims_use_nk <- 1:min(10, ncol(mfs2_c4_nk[["pca"]]@cell.embeddings))
mfs2_c4_nk <- FindNeighbors(mfs2_c4_nk, dims = dims_use_nk)
mfs2_c4_nk <- FindClusters(mfs2_c4_nk, resolution = 1.0)
mfs2_c4_nk <- RunUMAP(mfs2_c4_nk, dims = dims_use_nk)

DimPlot(mfs2_c4_nk, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 6: Average expression per subcluster ---
expr_mat_nk <- FetchData(mfs2_c4_nk, vars = nk_filtered)
avg_expr_nk <- expr_mat_nk %>%
  mutate(cluster = mfs2_c4_nk$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# --- Step 7: Map Ensembl IDs back to gene symbols if needed ---
ensembl_ids <- colnames(avg_expr_nk)[-1]
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = ensembl_ids,
                             keytype = "ENSEMBL",
                             columns = "SYMBOL")
id2symbol <- setNames(map$SYMBOL, map$ENSEMBL)
colnames(avg_expr_nk)[-1] <- ifelse(nzchar(id2symbol[colnames(avg_expr_nk)[-1]]),
                                    id2symbol[colnames(avg_expr_nk)[-1]],
                                    colnames(avg_expr_nk)[-1])

print(head(avg_expr_nk))

# --- Step 8: Annotate NK subtypes ---
# --- NK subtype annotation (aligned with nk_markers) ---

Activated_NK <- c("CD69","TNFRSF9","CRTAM","CLEC2D")

Memory_NK    <- c("RUNX3","ZEB2")

Signaling_NK <- c("STAT4","REL","PTPN22","DOCK2","ARHGAP25","CD247","CD96","KLRK1",
                  "IKZF1","IKZF3","PRKCQ","ITK","VAV1","SYTL3","SKAP1","CYTIP",
                  "TRAF1","NR3C1","TOX")

# Build lineage mapping
all_markers <- colnames(avg_expr_nk)[-1]
marker_lineage <- c(
  setNames(rep("Activated", length(Activated_NK)), Activated_NK),
  setNames(rep("Memory",    length(Memory_NK)),    Memory_NK),
  setNames(rep("Signaling", length(Signaling_NK)), Signaling_NK)
)
marker_lineage <- marker_lineage[all_markers]

annotation_row <- data.frame(Subtype = marker_lineage)
rownames(annotation_row) <- all_markers


# --- Step 9: Heatmap ---
mat <- as.matrix(avg_expr_nk[,-1]); rownames(mat) <- avg_expr_nk$cluster
vars <- apply(mat, 2, function(x) var(x, na.rm = TRUE)); vars[is.na(vars)] <- 0
scale_opt <- if (all(vars == 0)) "none" else "row"

pheatmap(t(mat),
         cluster_rows = FALSE,
         cluster_cols = (nrow(mat) > 1),
         scale = scale_opt,
         annotation_row = annotation_row,
         main = "NK cell subclusters with subtype annotation",
         cellheight = 12, cellwidth = 12, fontsize = 10)



library(Seurat)
library(dplyr)
library(pheatmap)
library(biomaRt)

# --- Step 0: Subset cluster 6 (cycling CD8 cells) ---
mfs2_c6 <- subset(mfs2_clean, idents = "6")
cat("Total cells in C6 (cluster 6):", ncol(mfs2_c6), "\n")

# --- Step 1: Define cycling CD8 marker set (30 genes) ---
cycling_cd8_markers <- c(
  "MKI67","TOP2A","BIRC5","CCNB1","CCNB2","CDK1","UBE2C","CENPF","CENPE",
  "CDC20","CDC6","CDC45","CDCA8","CDCA7","CDKN3","TYMS","TK1","RRM2","PRC1",
  "STMN1","ASF1B","FOXM1","AURKB","KIF11","KIF2C","KIF23","KIFC1","TPX2",
  "HELLS","DTL"
)

# --- Step 2: Map gene symbols to dataset IDs if needed ---
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")

convert_to_dataset_ids <- function(gene_list, dataset) {
  if (grepl("^ENSG", rownames(dataset)[1])) {
    mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                     filters="hgnc_symbol", values=gene_list, mart=mart)
    ids <- mapping$ensembl_gene_id[mapping$ensembl_gene_id %in% rownames(dataset)]
    return(ids)
  } else {
    return(gene_list[gene_list %in% rownames(dataset)])
  }
}

cycling_present <- convert_to_dataset_ids(cycling_cd8_markers, mfs2_c6)
if (length(cycling_present) == 0) stop("No cycling CD8 markers found in dataset")

# --- Step 3: Filter CD8 cycling cells (>10 markers expressed) ---
expr_mat_cycling <- FetchData(mfs2_c6, vars = cycling_present, slot = "counts")
cells_keep <- rownames(expr_mat_cycling)[rowSums(expr_mat_cycling > 1) > 10]
mfs2_c6 <- subset(mfs2_c6, cells = cells_keep)
cat("Cells retained in C6 after CD8 cycling marker filter (>10 markers):", ncol(mfs2_c6), "\n")

# --- Step 4: Variance filter ---
gene_vars_cycling <- apply(FetchData(mfs2_c6, vars = cycling_present), 2, var)
cycling_filtered <- names(gene_vars_cycling[gene_vars_cycling > 0])
if (length(cycling_filtered) == 0) cycling_filtered <- cycling_present

# --- Step 5: PCA/UMAP reclustering ---
mfs2_c6_cycling <- ScaleData(mfs2_c6, features = cycling_filtered)
mfs2_c6_cycling <- RunPCA(mfs2_c6_cycling, features = cycling_filtered)
dims_use_cycling <- 1:min(10, ncol(mfs2_c6_cycling[["pca"]]@cell.embeddings))
mfs2_c6_cycling <- FindNeighbors(mfs2_c6_cycling, dims = dims_use_cycling)
mfs2_c6_cycling <- FindClusters(mfs2_c6_cycling, resolution = 1.0)
mfs2_c6_cycling <- RunUMAP(mfs2_c6_cycling, dims = dims_use_cycling)

DimPlot(mfs2_c6_cycling, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 6: Average expression per reclustered subcluster ---
expr_mat_cycling <- FetchData(mfs2_c6_cycling, vars = cycling_present)
avg_expr_cycling <- expr_mat_cycling %>%
  mutate(cluster = mfs2_c6_cycling$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# --- Step 7: Convert Ensembl IDs back to gene symbols if needed ---
ensembl_ids <- colnames(avg_expr_cycling)[-1]
if (grepl("^ENSG", ensembl_ids[1])) {
  mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                   filters="ensembl_gene_id",
                   values=ensembl_ids,
                   mart=mart)
  id2symbol <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
  colnames(avg_expr_cycling)[-1] <- id2symbol[colnames(avg_expr_cycling)[-1]]
}

print(head(avg_expr_cycling))

# --- Step 8: Annotate cycling phases ---
G2M_phase <- c("CCNB1","CCNB2","CDK1","CDC20","UBE2C","CENPF","CENPE")
S_phase   <- c("MKI67","TOP2A","TYMS","RRM2","TK1","DTL","HELLS")
Mitotic   <- c("BIRC5","AURKB","KIF11","KIF2C","KIF23","TPX2","FOXM1")

all_markers <- colnames(avg_expr_cycling)[-1]
marker_lineage <- c(
  setNames(rep("G2M", length(G2M_phase)), G2M_phase),
  setNames(rep("S", length(S_phase)), S_phase),
  setNames(rep("Mitotic", length(Mitotic)), Mitotic)
)
marker_lineage <- marker_lineage[all_markers]

annotation_row <- data.frame(Subtype = marker_lineage)
rownames(annotation_row) <- all_markers

# --- Step 9: Heatmap ---
mat <- as.matrix(avg_expr_cycling[,-1]); rownames(mat) <- avg_expr_cycling$cluster
pheatmap(t(mat),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_row = annotation_row,
         main = "Cycling CD8 T cell subclusters with phase annotation",
         cellheight = 12, cellwidth = 12, fontsize = 10)


#C4 CD8
library(Seurat)
library(dplyr)
library(pheatmap)
library(biomaRt)

# --- Step 0: Check what IDs your dataset uses ---
head(rownames(mfs2_clean), 20)   # Are these ENSG IDs or gene symbols?

# --- Step 1: Define CD8 marker list (symbols, curated) ---
c4_cd8_markers <- c(
  "RUNX3","ZEB2","CRTAM","KLRK1","CD96","SLAMF7","TNFRSF9","CD247","CD69",
  "IFNG","REL","STAT4","STAT5B","NFATC2","FOXP1","BCL11B","IKZF3","PTPN22",
  "PRKCQ","ITK","CYTIP","CD44","CD6","CLEC2B","SP140"
)

# --- Step 2: Map gene symbols to dataset IDs if needed ---
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")

if (grepl("^ENSG", rownames(mfs2_clean)[1])) {
  mapping <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = c4_cd8_markers,
                   mart = mart)
  c4_present <- mapping$ensembl_gene_id[mapping$ensembl_gene_id %in% rownames(mfs2_clean)]
} else {
  c4_present <- c4_cd8_markers[c4_cd8_markers %in% rownames(mfs2_clean)]
}

if (length(c4_present) == 0) stop("No CD8 markers found in dataset")

# --- Step 3: Subset cluster 4 (CD8 cells) ---
mfs2_c4 <- subset(mfs2_clean, idents = "4")
cat("Total cells in C4 (cluster 4):", ncol(mfs2_c4), "\n")

# --- Step 4: Filter CD8-like cells (>10 markers expressed) ---
expr_mat_c4 <- FetchData(mfs2_c4, vars = c4_present, slot = "counts")
cells_keep <- rownames(expr_mat_c4)[rowSums(expr_mat_c4 > 1) > 5]
mfs2_c4 <- subset(mfs2_c4, cells = cells_keep)
cat("Cells retained in C4 after CD8 marker filter (>10 markers):", ncol(mfs2_c4), "\n")

# --- Step 5: Variance filter ---
gene_vars_c4 <- apply(FetchData(mfs2_c4, vars = c4_present), 2, var)
c4_filtered <- names(gene_vars_c4[gene_vars_c4 > 0])
if (length(c4_filtered) == 0) c4_filtered <- c4_present

# --- Step 6: PCA/UMAP clustering ---
mfs2_c4_cd8 <- ScaleData(mfs2_c4, features = c4_filtered)
mfs2_c4_cd8 <- RunPCA(mfs2_c4_cd8, features = c4_filtered)
dims_use_c4 <- 1:min(10, ncol(mfs2_c4_cd8[["pca"]]@cell.embeddings))
mfs2_c4_cd8 <- FindNeighbors(mfs2_c4_cd8, dims = dims_use_c4)
mfs2_c4_cd8 <- FindClusters(mfs2_c4_cd8, resolution = 1.0)
mfs2_c4_cd8 <- RunUMAP(mfs2_c4_cd8, dims = dims_use_c4)

DimPlot(mfs2_c4_cd8, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# --- Step 7: Average expression per subcluster ---
expr_mat_c4 <- FetchData(mfs2_c4_cd8, vars = c4_present)
avg_expr_c4 <- expr_mat_c4 %>%
  mutate(cluster = mfs2_c4_cd8$seurat_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# --- Step 8: Convert Ensembl IDs back to gene symbols if needed ---
ensembl_ids <- colnames(avg_expr_c4)[-1]
if (grepl("^ENSG", ensembl_ids[1])) {
  mapping_back <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                        filters = "ensembl_gene_id",
                        values = ensembl_ids,
                        mart = mart)
  id2symbol <- setNames(mapping_back$hgnc_symbol, mapping_back$ensembl_gene_id)
  colnames(avg_expr_c4)[-1] <- id2symbol[colnames(avg_expr_c4)[-1]]
}

# --- Step 9: Heatmap ---
mat <- as.matrix(avg_expr_c4[,-1])
rownames(mat) <- avg_expr_c4$cluster
pheatmap(t(mat),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "row",
         main = "CD8 T cell subclusters (Cluster 4, curated panel)",
         cellheight = 12, cellwidth = 12, fontsize = 10)

# --- Step 10: Heatmap with subtype annotation ---
gene_subtypes <- c(
  "RUNX3"="Effector","ZEB2"="Effector","CRTAM"="Effector","KLRK1"="Effector",
  "CD96"="Effector","SLAMF7"="Effector","TNFRSF9"="Effector","CD247"="Effector",
  "IFNG"="Effector","PRKCQ"="Effector","ITK"="Effector","CYTIP"="Effector",
  "CD69"="Activation","CD44"="Activation","CD6"="Activation","CLEC2B"="Activation",
  "FOXP1"="Regulatory","BCL11B"="Regulatory","IKZF3"="Regulatory","STAT4"="Regulatory",
  "STAT5B"="Regulatory","NFATC2"="Regulatory","REL"="Regulatory","SP140"="Regulatory",
  "PTPN22"="Regulatory"
)

genes_in_mat <- colnames(mat)
valid_genes <- genes_in_mat[genes_in_mat %in% names(gene_subtypes)]

annotation_row <- data.frame(Subtype = gene_subtypes[valid_genes])
rownames(annotation_row) <- valid_genes

pheatmap(t(mat[, valid_genes]),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_row = annotation_row,
         main = "CD8 T cell subclusters (Cluster 4, subtype-labeled)",
         cellheight = 12, cellwidth = 12, fontsize = 10)





#single cell GSEA
# Install ESCAPE if not already
install.packages("rlang")
packageVersion("rlang")

library(rlang)
devtools::install_github("BorchLab/escape")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("escape")
library(escape)
library(Seurat)
library(escape)
# Extract normalized expression matrix from Seurat
expr_matrix <- GetAssayData(mfs1_clean, assay = "RNA", layer = "data")


install.packages("vctrs")
install.packages("msigdbr")
library(msigdbr)


# Load the immune signature gene sets
GS.immune <- msigdbr(species = "Homo sapiens", category = "C7")
GS.immune.list   <- split(GS.immune$gene_symbol,   GS.immune$gs_name)

# Check what’s inside
names(GS.immune.list)[1:10]   # show first 10 gene set names

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")
library(GSVA)

# Run ESCAPE enrichment
mfs1_clean <- runEscape(
  obj = mfs1_clean,
  input.data = mfs1_clean, # required expression matrix
  method = "ssGSEA",
  gene.sets = GS.immune.list,           # immune signature database
  groups = 1000,                     # number of permutations
  min.size = 3,                      # minimum gene set size
  new.assay.name = "escape.immune"   # save enrichment scores as new assay
)

library(escape)
library(future)
plan("multicore", workers = 8)

mfs1_clean <- runEscape(
    obj = mfs1_clean,
    input.data = mfs1_clean,
    method = "ssGSEA",
    gene.sets = GS.immune.list,
    groups = 250,              # fewer permutations
    min.size = 10,             # skip tiny sets
    new.assay.name = "escape.immune"
)

# Check that enrichment scores are now stored
Assays(mfs1_clean)

# Check that enrichment scores are now stored
Assays(mfs1_clean)


# Add results back into Seurat metadata
mfs1_clean <- AddMetaData(mfs1_clean, immune_results)






C0_cd8_markers <- c(
  "CD8A","CD3D","CD3E","CD3G","GZMK","GZMA","GZMM","GZMH",
  "NKG7","CCL5","CST7","CTSW","CXCR4","CD69","CRTAM","RUNX3",
  "LAG3","CXCR3","IL7R","CD2","CD44","CD48","CD52","TRAC",
  "TRBC1","TRBC2","PTPRC","HCST"
)
Effector_CD8 <- c("GZMK","GZMA","GZMM","NKG7","CCL5","CST7","CTSW")
Memory_CD8   <- c("IL7R","CXCR3","CD44","CD69","CRTAM","RUNX3")
Naive_CD8    <- c("CCR7","LEF1","TCF7","SELL")
Exhausted_CD8<- c("LAG3","PDCD1","TIGIT","CTLA4")

C1_cd4_markers <- c(
  "IL7R","CCR7","CD28","ICOS","TNFRSF4","TNFRSF18","CTLA4","TIGIT",
  "CD2","CD5","CD44","CD52","CD3D","CD3E","CD3G","CD27","CD247","IL2RA",
  "FOXP1","SATB1","BCL11B","GATA3","STAT3","BATF","MAF","TOX"
)
Naive_CD4     <- c("CCR7","IL7R","CD27")
Memory_CD4    <- c("CD44","CD28","FOXP1","SATB1")
Effector_CD4  <- c("ICOS","TNFRSF4","TNFRSF18","CD52","CD2","CD5")
Regulatory_CD4<- c("CTLA4","IL2RA","TIGIT","TOX","BCL11B","GATA3")

C2_nk_markers <- c(
  "NKG7","GZMB","GZMA","GZMK","GZMH","PRF1","IFNG","CCL3","CCL4","CCL5",
  "XCL1","XCL2","CD69","CD8A","CD3D","CD3E","CD3G","CRTAM","CST7","CTSW",
  "KLRK1","SLAMF7","FASLG","RUNX3","ZEB2","TNFRSF9"
)
Cytotoxic_NK <- c("GZMB","GZMA","GZMK","GZMH","PRF1","NKG7","CST7","CTSW")
Activated_NK <- c("CD69","IFNG","TNFRSF9","FASLG","XCL1","XCL2","CCL3","CCL4","CCL5")
Memory_NK    <- c("RUNX3","ZEB2","SLAMF7","CRTAM")


C4_nk_markers <- c(
  "CD69","TNFRSF9","CRTAM","CLEC2D","RUNX3","ZEB2","STAT4","REL","PTPN22",
  "DOCK2","ARHGAP25","CD247","CD96","KLRK1","IKZF1","IKZF3","PRKCQ","ITK",
  "VAV1","SYTL3","SKAP1","CYTIP","TRAF1","NR3C1","TOX"
)
Activated_NK <- c("CD69","TNFRSF9","CRTAM","CLEC2D")

Memory_NK    <- c("RUNX3","ZEB2")

Signaling_NK <- c("STAT4","REL","PTPN22","DOCK2","ARHGAP25","CD247","CD96","KLRK1",
                  "IKZF1","IKZF3","PRKCQ","ITK","VAV1","SYTL3","SKAP1","CYTIP",
                  "TRAF1","NR3C1","TOX")

c4_cd8_markers <- c(
  "RUNX3","ZEB2","CRTAM","KLRK1","CD96","SLAMF7","TNFRSF9","CD247","CD69",
  "IFNG","REL","STAT4","STAT5B","NFATC2","FOXP1","BCL11B","IKZF3","PTPN22",
  "PRKCQ","ITK","CYTIP","CD44","CD6","CLEC2B","SP140"
)
gene_subtypes <- c(
  "RUNX3"="Effector","ZEB2"="Effector","CRTAM"="Effector","KLRK1"="Effector",
  "CD96"="Effector","SLAMF7"="Effector","TNFRSF9"="Effector","CD247"="Effector",
  "IFNG"="Effector","PRKCQ"="Effector","ITK"="Effector","CYTIP"="Effector",
  "CD69"="Activation","CD44"="Activation","CD6"="Activation","CLEC2B"="Activation",
  "FOXP1"="Regulatory","BCL11B"="Regulatory","IKZF3"="Regulatory","STAT4"="Regulatory",
  "STAT5B"="Regulatory","NFATC2"="Regulatory","REL"="Regulatory","SP140"="Regulatory",
  "PTPN22"="Regulatory"
)


C5_macrophage_core <- c("CD68","CD163","MRC1","MSR1","CSF1R","AIF1","HLA-DRA","HLA-DRB1","HLA-DPB1","C1QA","C1QB","C1QC")
m1_like <- c("IL1B","TNF","TLR2","TLR4","NLRP3","CXCL8","CXCL2","CD86")
m2_like <- c("IL10RA","HMOX1","MERTK","TGFBI","VSIG4","SIGLEC1","FOLR2","CCR1","STAB1")


C5_classical_mono    <- c("CD14","LYZ","CST3","FCGR3A","CD68","SERPINA1","C1QA","C1QB","C1QC")
intermediate_mono <- c("IL1B","TNFAIP2","CD86","PLAUR","CXCL8","TFEC","CSF1R","CD302")
nonclassical_mono <- c("FCGR1A","ITGAX","MRC1","NRP2","SIGLEC1","LILRB1","LILRB2","CXCL16","MAFB")


C6_cycling_cd8_markers <- c(
  "MKI67","TOP2A","BIRC5","CCNB1","CCNB2","CDK1","UBE2C","CENPF","CENPE",
  "CDC20","CDC6","CDC45","CDCA8","CDCA7","CDKN3","TYMS","TK1","RRM2","PRC1",
  "STMN1","ASF1B","FOXM1","AURKB","KIF11","KIF2C","KIF23","KIFC1","TPX2",
  "HELLS","DTL"
)
G2M_phase <- c("CCNB1","CCNB2","CDK1","CDC20","UBE2C","CENPF","CENPE")
S_phase   <- c("MKI67","TOP2A","TYMS","RRM2","TK1","DTL","HELLS")
Mitotic   <- c("BIRC5","AURKB","KIF11","KIF2C","KIF23","TPX2","FOXM1")





                             #single cell GSEA


library(Seurat)
library(fgsea)
library(dplyr)

# 1. Subset Seurat object to clusters 5 and 13
# Subset cluster 5 only
c5_subset <- subset(mfs1_clean, idents = 5)

# Subset cluster 13 only
c13_subset <- subset(mfs1_clean, idents = 13)

# 2. Reclustering within immune subset
# Reclustering C5
c5_subset <- FindVariableFeatures(c5_subset)
c5_subset <- ScaleData(c5_subset)
c5_subset <- RunPCA(c5_subset)
c5_subset <- FindNeighbors(c5_subset, dims = 1:20)
c5_subset <- FindClusters(c5_subset, resolution = 0.8)
c5_subset <- RunUMAP(c5_subset, dims = 1:20)

# Subset cluster 13 into CD4 and CD8 populations
cl13_cd4 <- subset(mfs1_clean, idents = 13, features = cd4_genes_filtered)
cl13_cd8 <- subset(mfs1_clean, idents = 13, features = cd8_genes_filtered)

# CD4 reclustering
cl13_cd4 <- ScaleData(cl13_cd4, features = cd4_genes_filtered)
cl13_cd4 <- RunPCA(cl13_cd4, features = cd4_genes_filtered)
dims_use_cd4 <- 1:min(10, ncol(cl13_cd4[["pca"]]@cell.embeddings))
cl13_cd4 <- FindNeighbors(cl13_cd4, dims = dims_use_cd4)
cl13_cd4 <- FindClusters(cl13_cd4, resolution = 0.8)
cl13_cd4 <- RunUMAP(cl13_cd4, dims = dims_use_cd4, n.neighbors = min(15, ncol(cl13_cd4) - 1))

# CD8 reclustering
cl13_cd8 <- ScaleData(cl13_cd8, features = cd8_genes_filtered)
cl13_cd8 <- RunPCA(cl13_cd8, features = cd8_genes_filtered)
dims_use_cd8 <- 1:min(10, ncol(cl13_cd8[["pca"]]@cell.embeddings))
cl13_cd8 <- FindNeighbors(cl13_cd8, dims = dims_use_cd8)
cl13_cd8 <- FindClusters(cl13_cd8, resolution = 0.8)
cl13_cd8 <- RunUMAP(cl13_cd8, dims = dims_use_cd8, n.neighbors = min(15, ncol(cl13_cd8) - 1))

# NK reclustering
cl13_nk <- ScaleData(cl13_nk, features = nk_genes_filtered)
cl13_nk <- RunPCA(cl13_nk, features = nk_genes_filtered)
dims_use_nk <- 1:min(10, ncol(cl13_nk[["pca"]]@cell.embeddings))
cl13_nk <- FindNeighbors(cl13_nk, dims = dims_use_nk)
cl13_nk <- FindClusters(cl13_nk, resolution = 1.0)
cl13_nk <- RunUMAP(cl13_nk, dims = dims_use_nk, n.neighbors = 10)

c5_mac <- ScaleData(c5_mac, features = macrophage_genes_filtered)
c5_mac <- RunPCA(c5_mac, features = macrophage_genes_filtered)
dims_use_mac <- 1:min(10, ncol(c5_mac[["pca"]]@cell.embeddings))
c5_mac <- FindNeighbors(c5_mac, dims = dims_use_mac)
c5_mac <- FindClusters(c5_mac, resolution = 0.6)
c5_mac <- RunUMAP(c5_mac, dims = dims_use_mac, n.neighbors = min(15, ncol(c5_mac) - 1))

c5_mac
c5mac_fgsea

install.packages('devtools')
devtools::install_github('immunogenomics/presto')
pathways <- gmtPathways("all_gene_lists.gmt")

run_fgsea_subclusters <- function(obj) {
  subcluster_fgsea <- list()
  for (cl in levels(Idents(obj))) {
    markers_cl <- FindMarkers(obj, ident.1 = cl, min.pct = 0.1)
    gene_stats <- markers_cl$avg_log2FC
    names(gene_stats) <- rownames(markers_cl)
    gene_stats <- tapply(gene_stats, names(gene_stats), max)
    ranks <- sort(gene_stats, decreasing = TRUE)
    fgseaRes_cl <- fgseaMultilevel(pathways = pathways, stats = ranks)
    subcluster_fgsea[[cl]] <- fgseaRes_cl
  }
  return(subcluster_fgsea)
}

cd4_fgsea <- run_fgsea_subclusters(cl13_cd4)
cd8_fgsea <- run_fgsea_subclusters(cl13_cd8)
c5mono_fgsea <- run_fgsea_subclusters(c5_mono)
c5mac_fgsea <- run_fgsea_subclusters(c5_mac)

cd8_fgsea

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnsDb.Hsapiens.v86")

library(EnsDb.Hsapiens.v86)
keep <- !is.na(rownames(c5_mac))
c5_mac <- c5_mac[keep, ]
clean_ids <- clean_ids[keep]
symbol_vec <- symbol_vec[keep]

# Clean IDs first (remove suffixes like "-U")
rownames(c5_mac) <- make.unique(rownames(c5_mac))

clean_ids <- gsub("-.*", "", rownames(c5_mac))
clean_ids <- clean_ids[!is.na(clean_ids)]

clean_ids <- gsub("-.*", "", rownames(c5_mac))

# Map Ensembl IDs to gene symbols
clean_ids <- gsub("-.*", "", rownames(c5_mac))

# Map only non-NA IDs
valid_ids <- clean_ids[!is.na(clean_ids)]
gene_map <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                  keys = valid_ids,
                                  keytype = "GENEID",
                                  columns = c("SYMBOL"))

# Align symbols back to all features
symbol_vec <- gene_map$SYMBOL[match(clean_ids, gene_map$GENEID)]



anno_df <- data.frame(
  ensembl = clean_ids,
  symbol  = symbol_vec,
  row.names = rownames(c5_mac)
)

Annotation(c5_mac, assay = "RNA") <- anno_df


# Attach to Seurat Assay5
Annotation(c5_mac[["RNA"]]) <- anno_df

# Force unique feature names directly
rownames(c5_mac[["RNA"]]) <- make.unique(rownames(c5_mac[["RNA"]]))

# Verify uniqueness
anyDuplicated(rownames(c5_mac[["RNA"]]))   # should return 0

# Extract raw counts
# Extract raw counts from the RNA assay
raw_counts <- LayerData(c5_mac[["RNA"]], layer = "counts")

# Check dimensions
dim(raw_counts)
  
# Force unique gene names
rownames(raw_counts) <- make.unique(rownames(raw_counts))

# Recreate Seurat object with metadata
c5_mac <- CreateSeuratObject(counts = raw_counts, meta.data = c5_mac@meta.data)

# Proceed with standard workflow
c5_mac <- NormalizeData(c5_mac)
c5_mac <- FindVariableFeatures(c5_mac)
c5_mac <- ScaleData(c5_mac)
c5_mac <- RunPCA(c5_mac)
c5_mac <- FindNeighbors(c5_mac, dims = 1:20)
c5_mac <- FindClusters(c5_mac, resolution = 0.5)

# Check identities
table(Idents(c5_mac))

# Run marker detection
markers_cl <- FindMarkers(c5_mac, ident.1 = 0, ident.2 = 1)





pathways <- gmtPathways("all_gene_lists.gmt")

run_fgsea_subclusters <- function(obj) {
  subcluster_fgsea <- list()
  for (cl in levels(Idents(obj))) {
    markers_cl <- FindMarkers(obj, ident.1 = cl, min.pct = 0.1)
    gene_stats <- markers_cl$avg_log2FC
    names(gene_stats) <- rownames(markers_cl)
    gene_stats <- tapply(gene_stats, names(gene_stats), max)
    ranks <- sort(gene_stats, decreasing = TRUE)
    fgseaRes_cl <- fgseaMultilevel(pathways = pathways, stats = ranks)
    subcluster_fgsea[[cl]] <- fgseaRes_cl
  }
  return(subcluster_fgsea)
}

nk_fgsea <- run_fgsea_subclusters(cl13_nk)




# CD4
pathways_cd4 <- unique(unlist(lapply(cd4_fgsea, function(x) x$pathway)))
cd4_NES <- sapply(cd4_fgsea, function(res) {
  nes <- setNames(res$NES, res$pathway)
  nes[pathways_cd4]
})
rownames(cd4_NES) <- pathways_cd4

# CD8
pathways_cd8 <- unique(unlist(lapply(cd8_fgsea, function(x) x$pathway)))
cd8_NES <- sapply(cd8_fgsea, function(res) {
  nes <- setNames(res$NES, res$pathway)
  nes[pathways_cd8]
})
rownames(cd8_NES) <- pathways_cd8

#NK
pathways_nk <- unique(unlist(lapply(nk_fgsea, function(x) x$pathway)))
nk_NES <- sapply(nk_fgsea, function(res) {
  nes <- setNames(res$NES, res$pathway)
  nes[pathways_nk]
})
rownames(nk_NES) <- pathways_nk


#mono
pathways_c5mono <- unique(unlist(lapply(c5mono_fgsea, function(x) x$pathway)))
c5mono_NES <- sapply(c5mono_fgsea, function(res) {
  nes <- setNames(res$NES, res$pathway)
  nes[pathways_c5mono]
})
rownames(c5mono_NES) <- pathways_c5mono



#mac


# Collect NES values per cluster
c5mac_NES_list <- lapply(c5mac_fgsea, function(res) {
  nes <- setNames(res$NES, res$pathway)
  nes[pathways_c5mac]   # align to full pathway list
})

# Combine into a matrix (fill missing with NA)
c5mac_NES <- do.call(cbind, c5mac_NES_list)

# Add rownames and colnames
rownames(c5mac_NES) <- pathways_c5mac
colnames(c5mac_NES) <- names(c5mac_fgsea)

pathways_c5mac <- unique(unlist(lapply(c5mac_fgsea, function(x) x$pathway)))
c5mac_NES <- sapply(c5mac_fgsea, function(res) {
  nes <- setNames(res$NES, res$pathway)
  nes[pathways_c5mac]
})
rownames(c5mac_NES) <- pathways_c5mac





# Clean NA values
nk_NES_clean <- nk_NES
nk_NES_clean[is.na(nk_NES_clean)] <- 0


c5mac_NES_top <- c5mac_NES_clean[top_pathways_c5mac, ]

c5mac_NES_clean <- c5mac_NES
c5mac_NES_clean[is.na(c5mac_NES_clean)] <- 0


c5mono_NES_clean <- c5mono_NES
c5mono_NES_clean[is.na(c5mono_NES_clean)] <- 0

c13_NES_clean <- c13_NES
c13_NES_clean[is.na(c13_NES_clean)] <- 0


# Clean NA values
nk_NES_clean <- nk_NES
nk_NES_clean[is.na(nk_NES_clean)] <- 0

library(pheatmap)
library(pheatmap)





# Subset NES matrix
c5_NES_top <- c5_NES_clean[top_pathways_c5, ]

pheatmap(c5_NES_top, scale = "row",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "C5 Subcluster Top 20 Pathways")

# Combine all fgsea results for CD4 subclusters
all_cd4 <- do.call(rbind, cd4_fgsea)

top_pathways_cd4 <- all_cd4 %>%
  arrange(padj) %>%
  head(20) %>%
  pull(pathway)

cd4_NES_top <- cd4_NES[top_pathways_cd4, ]

pheatmap(cd4_NES_top, scale = "row",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "C13 CD4 Subcluster Top 20 Pathways")

# Combine all fgsea results for CD8 subclusters
all_cd8 <- do.call(rbind, cd8_fgsea)

top_pathways_cd8 <- all_cd8 %>%
  arrange(padj) %>%
  head(20) %>%
  pull(pathway)

cd8_NES_top <- cd8_NES[top_pathways_cd8, ]

pheatmap(cd8_NES_top, scale = "row",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "C13 CD8 Subcluster Top 20 Pathways")
         
         
library(dplyr)

all_nk <- do.call(rbind, nk_fgsea)
top_pathways_nk <- all_nk %>%
  arrange(padj) %>%
  head(20) %>%
  pull(pathway)

nk_NES_top <- nk_NES_clean[top_pathways_nk, ]

library(pheatmap)

pheatmap(nk_NES_top, scale = "row",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "C13 NK Subcluster Top 20 Pathways")



c5mono_NES_top <- c5mono_NES_clean[top_pathways_c5mono, ]

all_c5mono <- do.call(rbind, c5mono_fgsea)
top_pathways_c5mono <- all_c5mono %>%
  arrange(padj) %>%
  head(20) %>%
  pull(pathway)


library(pheatmap)

pheatmap(c5mono_NES_top, scale = "row",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "C5 monocytes Subcluster Top 20 Pathways")




# Combine fgsea results into one data frame
all_c5mac <- do.call(rbind, c5mac_fgsea)

# Select top 20 pathways by adjusted p-value
top_pathways_c5mac <- all_c5mac %>%
  arrange(padj) %>%
  head(20) %>%
  pull(pathway)

# Build NES matrix
c5mac_NES_list <- lapply(c5mac_fgsea, function(res) {
  nes <- setNames(res$NES, res$pathway)
  nes[top_pathways_c5mac]   # align only to top pathways
})

c5mac_NES <- do.call(cbind, c5mac_NES_list)
rownames(c5mac_NES) <- top_pathways_c5mac
colnames(c5mac_NES) <- names(c5mac_fgsea)

# Replace NA with 0
c5mac_NES[is.na(c5mac_NES)] <- 0

# Plot heatmap
library(pheatmap)
pheatmap(c5mac_NES, scale = "row",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "C5 macrophages Subcluster Top 20 Pathways")






    avg_expr_complement,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = "Complement Activation Genes in C5 Macrophages",
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    width = 8, height = 10,
    color = colorRampPalette(c("navy","white","firebrick"))(50)
  )
} else {
  message("No complement activation genes with non-zero expression found.")
}



library(Seurat)
library(pheatmap)

# Define macrophage marker genes
macrophage_genes <- c(
  # Core macrophage markers
  "CD68","CD163","MRC1","MSR1","CSF1R","AIF1","FOLR2","VSIG4","STAB1",
  "HLA-DRA","HLA-DRB1","HLA-DPB1","C1QA","C1QB","C1QC",
  # M1-like (pro-inflammatory)
  "IL1B","TNF","TLR2","TLR4","NLRP3","CXCL8","CCR1","CD86",
  # M2-like (anti-inflammatory/tissue repair)
  "IL10RA","HMOX1","MERTK","TGFBI","CD163L1","MSR1","VSIG4"
)

# Combine pathway genes (top 30 or selected) with macrophage markers
genes_to_plot <- unique(c(macrophage_genes, top_genes))  # top_genes_pw from earlier selection

# Keep only genes present in Seurat object
genes_present <- intersect(genes_to_plot, rownames(c5mac_obj))

# Average expression per cluster
avg_expr <- AverageExpression(c5mac_obj, features = genes_present,
                              assays = "RNA", slot = "data")$RNA

# Clean matrix
avg_expr <- avg_expr[rowSums(avg_expr, na.rm = TRUE) > 0, , drop = FALSE]

# Row annotation: pathway membership
pathway_annotation <- data.frame(
  Pathway = rep(NA, nrow(avg_expr)),
  row.names = rownames(avg_expr)
)
for (pw in selected_pathways) {
  genes_pw <- pathways[[pw]]
  pathway_annotation$Pathway[rownames(avg_expr) %in% genes_pw] <- pw
}
# Mark macrophage genes
pathway_annotation$Pathway[rownames(avg_expr) %in% macrophage_genes] <- "Macrophage markers"

# Column annotation: cluster IDs
annotation_col <- data.frame(
  Cluster = colnames(avg_expr),
  row.names = colnames(avg_expr)
)

# Plot heatmap
pheatmap::pheatmap(
  avg_expr,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = pathway_annotation,
  annotation_col = annotation_col,
  main = "Pathway + Macrophage Marker Genes in C5 Macrophage Clusters",
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 45,
  width = 12, height = 14,
  color = colorRampPalette(c("navy","white","firebrick"))(50)
)




#####                          mfs1
library(Seurat)
library(fgsea)
library(dplyr)
library(ggplot2)
library(ggpubr)   # for stat_compare_means

plot_top_genes_per_pathway <- function(seurat_obj, fgsea_res, celltype_name, n_pathways = 5, n_genes = 5) {
  # Get top pathways by adjusted p-value
  top_pathways <- fgsea_res %>%
    arrange(padj) %>%
    head(n_pathways) %>%
    pull(pathway)
  
  # Open PDF to save plots
  pdf(paste0(celltype_name, "_top", n_pathways, "_pathways_top", n_genes, "genes.pdf"), width = 12, height = 8)
  
  for (pw in top_pathways) {
    # Get genes in pathway
    genes_pw <- pathways[[pw]]
    genes_present <- genes_pw[genes_pw %in% rownames(seurat_obj)]
    
    if (length(genes_present) > 0) {
      # Rank genes by average expression
      avg_expr <- rowMeans(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[genes_present, , drop = FALSE])
      top_genes <- names(sort(avg_expr, decreasing = TRUE))[1:min(n_genes, length(avg_expr))]
      
      # Violin plot with statistical tests
      p1 <- VlnPlot(seurat_obj, features = top_genes, group.by = "seurat_clusters", pt.size = 0) +
        ggtitle(paste(celltype_name, "-", pw, "- Violin")) +
        stat_compare_means(method = "wilcox.test", label = "p.signif")  # adds Wilcoxon test p-values
      
      print(p1)
      
      # UMAP plot
      p2 <- FeaturePlot(seurat_obj, features = top_genes) +
        ggtitle(paste(celltype_name, "-", pw, "- UMAP"))
      print(p2)
    }
  }
  
  dev.off()
}

# Example calls
plot_top_genes_per_pathway(cl13_cd4, all_cd4, "CD4")
plot_top_genes_per_pathway(cl13_cd8, all_cd8, "CD8")
plot_top_genes_per_pathway(cl13_nk, all_nk, "NK")
plot_top_genes_per_pathway(c5_mac, all_c5mac, "Macrophage")
plot_top_genes_per_pathway(c5_mono, all_c5mono, "Monocyte")

# Create results folder explicitly
dir.create("results", showWarnings = FALSE)

# Test PDF writing
pdf("results/test.pdf")
plot(1:10)
dev.off()

plot_top_genes_per_pathway <- function(seurat_obj, fgsea_res, celltype_name,
                                       n_pathways = 5, n_genes = 5, outdir = "results") {
  # Ensure output directory exists
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Get top pathways by adjusted p-value
  top_pathways <- fgsea_res %>%
    arrange(padj) %>%
    head(n_pathways) %>%
    pull(pathway)
  
  # Build output filename
  outfile <- file.path(outdir, paste0(celltype_name, "_top", n_pathways,
                                      "_pathways_top", n_genes, "genes.pdf"))
  
  # Open PDF
  pdf(outfile, width = 12, height = 8)
  
  for (pw in top_pathways) {
    genes_pw <- pathways[[pw]]
    genes_present <- genes_pw[genes_pw %in% rownames(seurat_obj)]
    
    if (length(genes_present) > 0) {
      avg_expr <- rowMeans(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[genes_present, , drop = FALSE])
      top_genes <- names(sort(avg_expr, decreasing = TRUE))[1:min(n_genes, length(avg_expr))]
      
      p1 <- VlnPlot(seurat_obj, features = top_genes, group.by = "seurat_clusters", pt.size = 0) +
        ggtitle(paste(celltype_name, "-", pw, "- Violin"))
      print(p1)
      
      p2 <- FeaturePlot(seurat_obj, features = top_genes, reduction = "umap") +
        ggtitle(paste(celltype_name, "-", pw, "- UMAP"))
      print(p2)
    }
  }
  
  dev.off()
  message("Saved plots to: ", outfile)
}


  
  
  
  
library(Seurat)
library(pheatmap)

# ---- Create output folder ----
outdir <- "Celltype_Heatmaps"
if (!dir.exists(outdir)) dir.create(outdir)

# ---- Define marker sets ----
cd4_genes <- c("IL7R","ICOS","CTLA4","PDCD1","TIGIT","TOX","BCL11B","GATA3",
               "CD6","CD69","CXCR4","IKZF2","BATF","STAT4","CD48","IL21R",
               "CCR7","LTB","PRKCQ","IL2RG","PTPN22","RUNX3","IKZF1","ARHGDIB",
               "SAMSN1","TNFRSF18","TNFRSF4","NFATC2")

cd8_genes <- c("CD2","CD3D","CD3E","CD3G","CD247","CD28","LCK","TRBC1","TRBC2",
               "RUNX3","SYTL3","ITK","SKAP1","FYB1","PRKCQ","STAT5B","STAT5A",
               "IKZF1","IKZF2","PTPN22","DOCK2","RAC2","GZMB","PRF1","TAGAP",
               "TNFRSF9","TNFRSF18")

nk_genes <- c("CD96","CST7","TNFRSF9","CLEC2D","HCST","KLRB1","SYTL3","FYB1",
              "SRGN","FYN","DOCK2","RAC2","GBP5","GNG2","TAGAP","TNFRSF4",
              "ISG20","CXCR4","ARHGAP15","STAT4","STAT5B","STAT5A","IKZF2",
              "NFATC2","ADGRE5","SERPINB9")

classical_mono <- c("CD14","LYZ","CST3","FCGR3A","CD68","SERPINA1","C1QA","C1QB","C1QC")
intermediate_mono <- c("IL1B","TNFAIP2","CD86","PLAUR","CXCL8","TFEC","CSF1R","CD302")
nonclassical_mono <- c("FCGR1A","ITGAX","MRC1","NRP2","SIGLEC1","LILRB1","LILRB2","CXCL16","MAFB")

# ---- Pathway sets ----
pathways_cd4 <- c("Alpha beta T cell activation",
                  "Leukocyte chemotaxis",
                  "B cell activation",
                  "Interleukin-4 and Interleukin-13 signaling",
                  "Interleukin-2 family signaling")

pathways_cd8 <- c("Gamma-delta T cell activation",
                  "TCR signaling",
                  "T cell proliferation", 
                  "Costimulation by the CD28 family")

pathways_nk <- c("Complement activation",
                 "Leukocyte migration involved in inflammatory response",
                 "Leukocyte mediated immunity",
                 "Leukocyte chemotaxis",
                 "Antiviral innate immune response")

pathways_mono <- c("MHC class II antigen presentation",
                   "Antigen processing and presentation",
                   "Activation of innate immune response",
                   "Interferon gamma signaling",
                   "Costimulation by the CD28 family")
                   
                   
library(dplyr)

unique(cd8_fgsea_res$pathway)[1:20]   # show first 20 names
grep("Immunoregulatory", cd8_fgsea_res$pathway, value = TRUE)

grep("Gamma-delta", cd8_fgsea_res$pathway, value = TRUE)






library(Seurat)
library(pheatmap)
library(gridExtra)

# ---- Create output folder ----
outdir <- "Celltype_Heatmaps"
if (!dir.exists(outdir)) dir.create(outdir)

# ---- Function to plot and save heatmap (restricted to ~30 genes) ----
plot_celltype_heatmap <- function(seurat_obj, marker_genes, pathways, selected_pathways, filename, title) {
  
  # Collect pathway genes
  genes_selected <- unique(unlist(pathways[selected_pathways]))
  
  # Combine with marker genes
  genes_to_plot <- unique(c(marker_genes, genes_selected))
  
  # Keep only genes present
  genes_present <- intersect(genes_to_plot, rownames(seurat_obj))
  
  # Average expression per cluster
  avg_expr <- AverageExpression(seurat_obj, features = genes_present,
                                assays = "RNA", slot = "data")$RNA
  
  # Clean matrix
  avg_expr <- avg_expr[rowSums(avg_expr, na.rm = TRUE) > 0, , drop = FALSE]
  
  # ---- Restrict to top 30 variable genes ----
  gene_var <- apply(avg_expr, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(30, length(gene_var))]
  avg_expr <- avg_expr[top_genes, , drop = FALSE]
  
  # Row annotation
  pathway_annotation <- data.frame(
    Pathway = rep(NA, nrow(avg_expr)),
    row.names = rownames(avg_expr)
  )
  for (pw in selected_pathways) {
    genes_pw <- pathways[[pw]]
    pathway_annotation$Pathway[rownames(avg_expr) %in% genes_pw] <- pw
  }
  pathway_annotation$Pathway[rownames(avg_expr) %in% marker_genes] <- "Cell-type markers"
  
  # Column annotation
  annotation_col <- data.frame(
    Cluster = colnames(avg_expr),
    row.names = colnames(avg_expr)
  )
  
  # Generate heatmap object
  ht <- pheatmap::pheatmap(
    avg_expr,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_row = pathway_annotation,
    annotation_col = annotation_col,
    main = title,
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45,
    color = colorRampPalette(c("navy","white","firebrick"))(50),
    silent = TRUE
  )
  
  # Save as PNG
  png(file.path(outdir, filename), width = 1600, height = 1800, res = 200)
  grid::grid.draw(ht$gtable)
  dev.off()
}

plot_celltype_heatmap(cl13_cd4, cd4_genes, pathways, pathways_cd4,
                      "CD4_heatmap.png", "Pathway + CD4 Marker Genes in cl13_cd4 Clusters")

plot_celltype_heatmap(cl13_cd8, cd8_genes, pathways, pathways_cd8,
                      "CD8_heatmap.png", "Pathway + CD8 Marker Genes in cl13_cd8 Clusters")

plot_celltype_heatmap(cl13_nk, nk_genes, pathways, pathways_nk,
                      "NK_heatmap.png", "Pathway + NK Marker Genes in cl13_nk Clusters")

plot_celltype_heatmap(c5_mono, classical_mono, pathways, pathways_mono,
                      "Monocyte_heatmap.png", "Pathway + Monocyte Marker Genes in c5_mono Clusters")





library(Seurat)
library(fgsea)
library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ---- Force fgsea to use a stable local temp directory ----
tempdir_path <- "D:/Rtmp_fgsea_local"   # plain local folder, not synced
dir.create(tempdir_path, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = tempdir_path)
options(tempdir = tempdir_path)

# ---- Load pathways once ----
pathways <- gmtPathways("all_gene_lists.gmt")

# ---- Function: run fgsea per subcluster ----
run_fgsea_subclusters <- function(obj, resolution = 0.5) {
  raw_counts <- LayerData(obj[["RNA"]], layer = "counts")
  rownames(raw_counts) <- make.unique(rownames(raw_counts))
  
  # Map Ensembl IDs to gene symbols
  ensembl_ids <- rownames(raw_counts)
  symbol_map <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  gene_symbols <- ifelse(is.na(symbol_map), ensembl_ids, symbol_map)
  
  # Collapse duplicates
  raw_counts <- rowsum(raw_counts, group = gene_symbols)
  rownames(raw_counts) <- make.unique(rownames(raw_counts))
  
  # Rebuild Seurat object
  obj <- CreateSeuratObject(counts = raw_counts, meta.data = obj@meta.data)
  
  # Standard workflow
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 30)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = 1:20, umap.method = "uwot",
                 metric = "euclidean", nn.method = "fnn",
                 n.neighbors = 15, n.epochs = 200)
  
  # Run fgsea per subcluster (simple version only, single-thread)
  subcluster_fgsea <- list()
  for (cl in levels(Idents(obj))) {
    markers_cl <- tryCatch(
      FindMarkers(obj, ident.1 = cl, min.pct = 0.1, logfc.threshold = 0),
      error = function(e) NULL
    )
    if (is.null(markers_cl) || nrow(markers_cl) == 0) next
    
    gene_stats <- markers_cl$avg_log2FC
    names(gene_stats) <- rownames(markers_cl)
    gene_stats <- tapply(gene_stats, names(gene_stats), max)
    ranks <- sort(gene_stats, decreasing = TRUE)
    if (length(ranks) == 0) next
    
    fgseaRes_cl <- fgseaSimple(pathways = pathways, stats = ranks, nperm = 10000, nproc = 1)
    subcluster_fgsea[[cl]] <- fgseaRes_cl
  }
  return(list(seurat_obj = obj, fgsea_res = subcluster_fgsea))
}

# ---- Function: plot top 20 pathways ----
plot_top20_pathways <- function(results, celltype_name, outdir = "D:/Rtmp_fgsea_local/output") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  fgsea_res_list <- results$fgsea_res
  if (length(fgsea_res_list) == 0) {
    message("No fgsea results for ", celltype_name)
    return(NULL)
  }
  
  all_res <- do.call(rbind, fgsea_res_list)
  if (is.null(all_res) || nrow(all_res) == 0) {
    message("No enriched pathways for ", celltype_name)
    return(NULL)
  }
  
  top20 <- all_res %>% arrange(padj) %>% head(20)
  top_pathways <- top20$pathway
  
  nes_list <- lapply(fgsea_res_list, function(res) {
    nes <- setNames(res$NES, res$pathway)
    nes[top_pathways]
  })
  nes_mat <- do.call(cbind, nes_list)
  rownames(nes_mat) <- top_pathways
  colnames(nes_mat) <- names(fgsea_res_list)
  nes_mat[is.na(nes_mat)] <- 0
  
  # Save PDF
  pdf_file <- file.path(outdir, paste0(celltype_name, "_top20_pathways.pdf"))
  pdf(pdf_file, width = 10, height = 8)
  pheatmap(nes_mat, scale = "row",
           show_rownames = TRUE, show_colnames = TRUE,
           main = paste(celltype_name, "- Top 20 Pathways"))
  dev.off()
  
  # Save PNG as fallback
  png_file <- file.path(outdir, paste0(celltype_name, "_top20_pathways.png"))
  png(png_file, width = 1200, height = 1000)
  pheatmap(nes_mat, scale = "row",
           show_rownames = TRUE, show_colnames = TRUE,
           main = paste(celltype_name, "- Top 20 Pathways"))
  dev.off()
  
  message("Saved top20 pathways PDF and PNG for ", celltype_name,
          "\nPDF: ", pdf_file,
          "\nPNG: ", png_file)
}

# ---- Run one cluster at a time ----
res_c0_cd8 <- run_fgsea_subclusters(mfs2_c0_cd8, resolution = 0.6)
plot_top20_pathways(res_c0_cd8, "c0_CD8")







Mfs2 cluster subseted seurat done 

###repeat same analysis but with mfs2_clean this time and by firstly do fgsea first on each cell type(Mfs2 cluster subseted seurat done already, use seurat mfs2_c0_cd8, mfs2_c1_cd4, mfs2_c2_nk, mfs2_c4_nk, mfs2_c4_cd8, mfs2_c5_macro, mfs2_c5_mono, mfs2_c6_cycling, identified the top 5 pathways(possibly plot top 20 pathways first by fgsea) and then do plots on respective cell types with top 5 pathways + cell type annotation within one heatmap for each cell type, note all the markers here are for mfs2_clean seurat project(first half of code no need to change these markers)


mfs2_c0_cd8, 

mfs2_c1_cd4, 
mfs2_c2_nk, 
mfs2_c4_nk, 
mfs2_c4_cd8, 
mfs2_c5_macro, 
mfs2_c5_mono, 
mfs2_c6_cycling


library(Seurat)
library(fgsea)
library(dplyr)
library(pheatmap)

# ---- Function: select pathway genes ----
select_pathway_genes <- function(fgsea_res_list, top_n = 20, select_n = 10, genes_per_pathway = 3) {
  all_res <- do.call(rbind, fgsea_res_list)
  if (is.null(all_res) || nrow(all_res) == 0) return(NULL)
  
  top_res <- all_res %>% arrange(padj) %>% head(top_n)
  selected_pathways <- head(top_res$pathway, select_n)
  
  pathway_genes <- list()
  for (pw in selected_pathways) {
    leading_genes <- top_res$leadingEdge[top_res$pathway == pw][[1]]
    pathway_genes[[pw]] <- head(leading_genes, genes_per_pathway)
  }
  return(pathway_genes)
}

# ---- Function: plot heatmap ----
plot_celltype_heatmap <- function(seurat_obj, marker_genes, pathway_genes, filename, title) {
  genes_to_plot <- unique(c(marker_genes, unlist(pathway_genes)))
  genes_present <- intersect(genes_to_plot, rownames(seurat_obj))
  
  avg_expr <- AverageExpression(seurat_obj, features = genes_present,
                                assays = "RNA", slot = "data")$RNA
  avg_expr <- avg_expr[rowSums(avg_expr, na.rm = TRUE) > 0, , drop = FALSE]
  
  # Limit to top 40 variable genes for clarity
  gene_var <- apply(avg_expr, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(40, length(gene_var))]
  avg_expr <- avg_expr[top_genes, , drop = FALSE]
  
  # Row annotation: pathway vs marker
  pathway_annotation <- data.frame(Category = rep("Marker", nrow(avg_expr)),
                                   row.names = rownames(avg_expr))
  for (pw in names(pathway_genes)) {
    pathway_annotation$Category[rownames(avg_expr) %in% pathway_genes[[pw]]] <- pw
  }
  
  annotation_col <- data.frame(Cluster = colnames(avg_expr),
                               row.names = colnames(avg_expr))
  
  # Save high‑resolution PNG
  png(filename, width = 2000, height = 2200, res = 200)
  pheatmap(avg_expr,
           scale = "row",
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_row = pathway_annotation,
           annotation_col = annotation_col,
           main = title,
           fontsize_row = 10,
           fontsize_col = 12,
           angle_col = 45,
           color = colorRampPalette(c("navy","white","firebrick"))(50))
  dev.off()
}

# ---- Example: CD8 cluster ----
fgsea_cd8 <- res_c0_cd8$fgsea_res
pathway_genes_cd8 <- select_pathway_genes(fgsea_cd8, top_n = 20, select_n = 10, genes_per_pathway = 3)

plot_celltype_heatmap(res_c0_cd8$seurat_obj,
                      C0_cd8_markers,
                      pathway_genes_cd8,
                      "CD8_heatmap.png",
                      "Top Pathways + CD8 Markers")

# ---- Example: CD8 cluster ----
fgsea_cd8_c4<- res_c4_cd8$fgsea_res
pathway_genes_cd8_c4 <- select_pathway_genes(fgsea_cd8_c4, top_n = 20, select_n = 10, genes_per_pathway = 3)

plot_celltype_heatmap(res_c4_cd8$seurat_obj,
                      c4_cd8_markers,
                      pathway_genes_cd8_c4,
                      "c4_CD8_heatmap.png",
                      "Top Pathways + CD8 Markers")


# CD4 cluster
fgsea_cd4 <- res_c1_cd4$fgsea_res
pathway_genes_cd4 <- select_pathway_genes(fgsea_cd4, top_n = 20, select_n = 10, genes_per_pathway = 3)

plot_celltype_heatmap(res_c1_cd4$seurat_obj,
                      C1_cd4_markers,
                      pathway_genes_cd4,
                      "CD4_heatmap.png",
                      "Top Pathways + CD4 Markers")

# NK cluster
fgsea_nk <- res_c2_nk$fgsea_res
pathway_genes_nk <- select_pathway_genes(fgsea_nk)

plot_celltype_heatmap(res_c2_nk$seurat_obj,
                      C2_nk_markers,
                      pathway_genes_nk,
                      "NK_heatmap.png",
                      "Top Pathways + NK Markers")


# NK cluster
fgsea_nk_c4 <- res_c4_nk$fgsea_res
pathway_genes_nk_c4 <- select_pathway_genes(fgsea_nk_c4)

plot_celltype_heatmap(res_c4_nk$seurat_obj,
                      c4_nk_markers,
                      pathway_genes_nk_c4,
                      "NK_heatmap.png",
                      "Top Pathways + NK Markers")


fgsea_macro <- res_c5_macro$fgsea_res
pathway_genes_macro <- select_pathway_genes(fgsea_macro)

plot_celltype_heatmap(res_c5_macro$seurat_obj,
                      C5_macrophage_core,
                      pathway_genes_macro,
                      "Macro_heatmap.png",
                      "Top Pathways + Macrophage Markers")

fgsea_mono <- res_c5_mono$fgsea_res
pathway_genes_mono <- select_pathway_genes(fgsea_mono)

plot_celltype_heatmap(res_c5_mono$seurat_obj,
                      C5_classical_mono,
                      pathway_genes_mono,
                      "Mono_heatmap.png",
                      "Top Pathways + Monocyte Markers")

fgsea_cycling <- res_c6_cycling$fgsea_res
pathway_genes_cycling <- select_pathway_genes(fgsea_cycling)

plot_celltype_heatmap(res_c6_cycling$seurat_obj,
                      C6_cycling_cd8_markers,
                      pathway_genes_cycling,
                      "CyclingCD8_heatmap.png",
                      "Top Pathways + Cycling CD8 Markers")






##InferCNV


install.packages("pak")
library(pak)
library(mvtnorm)

library(rjags)

pak::pak("broadinstitute/infercnv")
install.packages("mvtnorm", type = "binary")
pkgbuild::check_build_tools(debug = TRUE)

system("python --version")
system("python -c \"import argparse, json, simplejson; print('OK')\"")
options(python_cmd = "D:/download/miniconda/python.exe")
Sys.setenv(PYTHON = "D:/download/miniconda/python.exe")
system("D:/download/miniconda/python.exe --version")
system("D:/download/miniconda/python.exe -c \"import argparse, json, simplejson; print('OK')\"")
library(reticulate)
use_python("D:/download/miniconda/python.exe", required = TRUE)
findpython::find_python_cmd(required_modules = c("argparse","json","simplejson"))

Sys.setenv(PYTHON = "D:/download/miniconda/python.exe")
Sys.setenv(PYTHON3 = "D:/download/miniconda/python.exe")
Sys.setenv(PYTHON_CMD = "D:/download/miniconda/python.exe")

library(infercnv)



library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
library(Seurat)
counts_mfs1 <- GetAssayData(mfs1_clean, assay = "RNA", layer = "counts")

# Combine all gene symbols from both Seurat objects
all_genes <- unique(c(rownames(mfs1_clean), rownames(mfs2_clean)))

# Query org.Hs.eg.db for chromosome positions
gene_pos_all <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = all_genes,
                                      columns = c("SYMBOL","CHR","CHRLOC","CHRLOCEND"),
                                      keytype = "SYMBOL")

# Remove rows with missing chromosome info
gene_pos_all <- gene_pos_all[!is.na(gene_pos_all$CHR), ]

# Convert negative coordinates to positive
gene_pos_all$CHRLOC <- abs(gene_pos_all$CHRLOC)
gene_pos_all$CHRLOCEND <- abs(gene_pos_all$CHRLOCEND)

# Remove duplicates
gene_pos_all <- gene_pos_all[!duplicated(gene_pos_all$SYMBOL), ]

# Keep only the columns inferCNV needs
gene_pos_final <- gene_pos_all[, c("SYMBOL","CHR","CHRLOC","CHRLOCEND")]

# Save to file
write.table(gene_pos_final, "gene_pos_combined.txt", sep="\t", quote=FALSE, row.names=FALSE)


colnames(gene_pos_combined)
gene_pos_combined <- read.table("gene_pos_combined.txt", header=TRUE)

colnames(gene_pos_combined) <- c("hgnc_symbol","chromosome_name","start_position","end_position")

# Now check overlap
sum(rownames(counts_mfs1) %in% gene_pos_combined$hgnc_symbol)
sum(rownames(counts_mfs2) %in% gene_pos_combined$hgnc_symbol)

head(gene_pos_combined$hgnc_symbol, 20)
head(rownames(counts_mfs1), 20)

head(rownames(counts_mfs2), 20)
library(org.Hs.eg.db)

ids2 <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = rownames(counts_mfs2),
                              columns = "SYMBOL",
                              keytype = "ENSEMBL")

# Replace rownames with gene symbols
rownames(counts_mfs2) <- ids2$SYMBOL[match(rownames(counts_mfs2), ids2$ENSEMBL)]

gene_pos_combined <- read.table("gene_pos_combined.txt", header=TRUE)
sum(rownames(counts_mfs1) %in% gene_pos_combined$hgnc_symbol)
sum(rownames(counts_mfs2) %in% gene_pos_combined$hgnc_symbol)

annots1 <- data.frame(
  Cell = colnames(counts_mfs1),
  Cluster = mfs1_clean$SingleR.labels,
  stringsAsFactors = FALSE
)

# Strip whitespace just in case
annots1$Cell <- trimws(annots1$Cell)

write.table(annots1,
            file = "mfs1_annotations.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

read.delim("mfs1_annotations.txt", nrows = 5)
identical(read.delim("mfs1_annotations.txt")$Cell, colnames(counts_mfs1))

length(colnames(counts_mfs1))
nrow(annots1)
test_annots <- read.delim("mfs1_annotations.txt", header=TRUE, sep="\t")
head(test_annots)
identical(test_annots$Cell, colnames(counts_mfs1))

write.table(annots1,
            file = "mfs1_annotations.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "UTF-8")

# Force clean ASCII, Unix line endings, no row names
write.table(annots1,
            file = "mfs1_annotations.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "ASCII",
            eol = "\n")

readLines("mfs1_annotations.txt", 3)
all(colnames(counts_mfs1) %in% annots1$Cell)


test_annots <- read.delim("mfs1_annotations.txt", header=TRUE, sep="\t")
identical(test_annots$Cell, colnames(counts_mfs1))

anyDuplicated(annots1$Cell)
class(counts_mfs1)

unique(annots1$Cluster)
setdiff(c("CD4+ T-cells","CD8+ T-cells","B-cells","NK cells","Monocytes","Macrophages","DC"),
        unique(annots1$Cluster))

table(annots1$Cluster)



gene_pos_combined <- read.table("gene_pos_combined.txt", header=TRUE)

# Rename columns correctly
colnames(gene_pos_combined) <- c("hgnc_symbol","chromosome_name","start_position","end_position")

# Force numeric conversion
gene_pos_combined$start_position <- as.numeric(gene_pos_combined$start_position)
gene_pos_combined$end_position   <- as.numeric(gene_pos_combined$end_position)

# Write back clean file
write.table(gene_pos_combined,
            file = "gene_pos_combined.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# Remove rows with NA in start or end position
gene_pos_combined_clean <- gene_pos_combined[!is.na(gene_pos_combined$start_position) &
                                               !is.na(gene_pos_combined$end_position), ]

# Double-check structure
str(gene_pos_combined_clean)

# Write back clean file
write.table(gene_pos_combined_clean,
            file = "gene_pos_combined.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
sum(is.na(gene_pos_combined_clean$start_position))
sum(is.na(gene_pos_combined_clean$end_position))

str(gene_pos_combined)

# Clean chromosome names
gene_pos_combined_clean$chromosome_name <- gsub("^chr", "", gene_pos_combined_clean$chromosome_name)

# Keep only standard chromosomes
valid_chroms <- c(as.character(1:22), "X", "Y", "MT")
gene_pos_combined_clean <- gene_pos_combined_clean[gene_pos_combined_clean$chromosome_name %in% valid_chroms, ]

# Write back clean file
write.table(gene_pos_combined_clean,
            file = "gene_pos_combined.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

unique(gene_pos_combined_clean$chromosome_name)

# Convert chromosome_name to factor with correct ordering
valid_chroms <- c(as.character(1:22), "X", "Y", "MT")
gene_pos_combined_clean$chromosome_name <- factor(
  gene_pos_combined_clean$chromosome_name,
  levels = valid_chroms,
  ordered = TRUE
)

# Write back clean file
write.table(gene_pos_combined_clean,
            file = "gene_pos_combined.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


gene_pos_combined_clean$chromosome_name <- as.character(gene_pos_combined_clean$chromosome_name)

# Map chromosomes to numeric values, keep X/Y/MT as strings
gene_pos_combined_clean$chromosome_name[gene_pos_combined_clean$chromosome_name == "X"] <- "23"
gene_pos_combined_clean$chromosome_name[gene_pos_combined_clean$chromosome_name == "Y"] <- "24"
gene_pos_combined_clean$chromosome_name[gene_pos_combined_clean$chromosome_name == "MT"] <- "25"

# Convert to numeric
gene_pos_combined_clean$chromosome_name <- as.numeric(gene_pos_combined_clean$chromosome_name)

# Write back clean file
write.table(gene_pos_combined_clean,
            file = "gene_pos_combined.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

str(gene_pos_combined_clean)
unique(gene_pos_combined_clean$chromosome_name)

is.numeric(gene_pos_combined_clean$chromosome_name)
any(is.na(gene_pos_combined_clean$chromosome_name))

length(rownames(counts_mfs1))
sum(rownames(counts_mfs1) %in% gene_pos_combined_clean$hgnc_symbol)
setdiff(head(rownames(counts_mfs1), 20), gene_pos_combined_clean$hgnc_symbol)

library(biomaRt)

library(biomaRt)

library(AnnotationHub)


# Option 2: Ensembl EnsDb (release 86)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
gene_ranges <- genes(EnsDb.Hsapiens.v86)
head(gene_ranges)

library(EnsDb.Hsapiens.v86)

gene_ranges <- genes(EnsDb.Hsapiens.v86)
gene_pos <- as.data.frame(gene_ranges)[, c("symbol","seqnames","start","end")]
colnames(gene_pos) <- c("hgnc_symbol","chromosome_name","start_position","end_position")

gene_pos$chromosome_name <- gsub("^chr", "", gene_pos$chromosome_name)
gene_pos$chromosome_name[gene_pos$chromosome_name == "X"] <- "23"
gene_pos$chromosome_name[gene_pos$chromosome_name == "Y"] <- "24"
gene_pos$chromosome_name[gene_pos$chromosome_name == "MT"] <- "25"
gene_pos$chromosome_name <- as.numeric(gene_pos$chromosome_name)

gene_pos <- gene_pos[gene_pos$hgnc_symbol %in% rownames(counts_mfs1), ]
# Remove duplicated gene symbols
gene_pos <- gene_pos[!duplicated(gene_pos$hgnc_symbol), ]

gene_pos$start_position <- as.numeric(gene_pos$start_position)
gene_pos$end_position   <- as.numeric(gene_pos$end_position)
gene_pos <- gene_pos[!is.na(gene_pos$chromosome_name) &
                       !is.na(gene_pos$start_position) &
                       !is.na(gene_pos$end_position), ]
gene_pos <- gene_pos[gene_pos$hgnc_symbol %in% rownames(counts_mfs1), ]
gene_pos <- gene_pos[!duplicated(gene_pos$hgnc_symbol), ]
write.table(gene_pos, "gene_pos_mfs1.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
sum(rownames(counts_mfs1) %in% gene_pos$hgnc_symbol)
anyDuplicated(gene_pos$hgnc_symbol)
str(gene_pos)
matched_genes <- intersect(rownames(counts_mfs1), gene_pos$hgnc_symbol)
counts_mfs1 <- counts_mfs1[matched_genes, ]
gene_pos <- gene_pos[gene_pos$hgnc_symbol %in% matched_genes, ]


# Sort gene order file by chromosome and start
gene_pos <- gene_pos[order(gene_pos$chromosome_name, gene_pos$start_position), ]

# Reorder counts matrix to match gene order
counts_mfs1 <- counts_mfs1[gene_pos$hgnc_symbol, ]

# Save gene order file again
write.table(gene_pos, "gene_pos_mfs1.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

anyDuplicated(gene_pos$hgnc_symbol)


library(org.Hs.eg.db)
mapped <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = rownames(counts_mfs1),
                                columns = "ENSEMBL",
                                keytype = "SYMBOL")
head(rownames(counts_mfs1))

str(gene_pos)

sum(rownames(counts_mfs1) %in% gene_pos$hgnc_symbol)

sum(rownames(counts_mfs1) %in% gene_pos$hgnc_symbol)

annots1 <- data.frame(
  colnames(counts_mfs1),
  mfs1_clean$SingleR.labels
)

counts_mfs1 <- as.matrix(counts_mfs1)
storage.mode(counts_mfs1) <- "numeric"
str(counts_mfs1[1:5, 1:5])

write.table(annots1,
            "mfs1_annotations.txt",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
readLines("mfs1_annotations.txt", n=5)


# Ensure perfect alignment
matched_genes <- intersect(rownames(counts_mfs1), gene_pos$hgnc_symbol)
counts_mfs1 <- counts_mfs1[matched_genes, ]
gene_pos <- gene_pos[match(matched_genes, gene_pos$hgnc_symbol), ]

# Check alignment
stopifnot(all(rownames(counts_mfs1) == gene_pos$hgnc_symbol))

# Save gene order file
write.table(gene_pos, "gene_pos_mfs1.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
readLines("mfs1_annotations.txt", n=5)
unique(annots1$Cluster)

cluster_vector <- mfs1_clean$SingleR.labels

annots1 <- data.frame(
  Cell = colnames(counts_mfs1),
  Cluster = cluster_vector,   # replace with your cluster assignments
  stringsAsFactors = FALSE
)

# Save with tab delimiter, no header
write.table(annots1,
            "mfs1_annotations.txt",
            sep="\t", quote=FALSE,
            row.names=FALSE,
            col.names=FALSE,
            fileEncoding="UTF-8")
annots_check <- read.table("mfs1_annotations.txt", sep="\t", header=FALSE)
str(annots_check)
head(annots_check)
unique(annots_check$V2)
setdiff(c("CD4+ T-cells","CD8+ T-cells","B-cells","NK cells","Monocytes","Macrophages","DC"),
        unique(annots_check$V2))
gene_pos$chromosome_name <- paste0("chr", gene_pos$chromosome_name)

write.table(gene_pos, "gene_pos_mfs1.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# Convert numeric chromosome names to "chrX" format
gene_pos$chromosome_name <- paste0("chr", gene_pos$chromosome_name)

# Save again
write.table(gene_pos, "gene_pos_mfs1.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


head(readLines("mfs1_annotations.txt", n=5))
unique(annots1$Cluster)
str(gene_pos)
# Keep only standard human chromosomes
# Define standard human chromosomes
standard_chrs <- paste0("chr", c(1:22, "X", "Y"))

# Keep only those
gene_pos <- gene_pos[gene_pos$chromosome_name %in% standard_chrs, ]

# Save again
write.table(gene_pos, "gene_pos_mfs1.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

unique(gene_pos$chromosome_name)
# Should show chr1–chr22, chrX, chrY only

# Reorder columns to match InferCNV expectation
gene_pos_clean <- gene_pos[, c("hgnc_symbol", "chromosome_name", "start_position", "end_position")]

# Save with tab delimiter, no header
write.table(gene_pos_clean, "gene_pos_mfs1.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Check for duplicates
sum(duplicated(gene_pos$hgnc_symbol))

# Check for NA values
colSums(is.na(gene_pos))

# Check chromosome names
unique(gene_pos$chromosome_name)

library(infercnv)
counts_mfs1 <- counts_mfs1[rowSums(counts_mfs1) > 10, ]

# Retry InferCNV
infercnv_obj1 <- CreateInfercnvObject(
  raw_counts_matrix = counts_mfs1,
  annotations_file = "mfs1_annotations.txt",
  gene_order_file = "gene_pos_mfs1.txt",
  ref_group_names = c("CD4+ T-cells","CD8+ T-cells","B-cells","NK cells","Monocytes","Macrophages","DC")
)

infercnv_obj1 <- infercnv::run(
  infercnv_obj1,
  cutoff = 0.1,
  out_dir = "infercnv_output_mfs1",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE
)
# If it's an RDS file
tumor_subclusters <- readRDS("infercnv_output_mfs1/15_tumor_subclusters.leiden.infercnv_obj")
class(tumor_subclusters)
str(tumor_subclusters)

# See what slots are available
slotNames(tumor_subclusters)

# Look at the structure
str(tumor_subclusters, max.level = 1)

# Access the tumor_subclusters slot
clusters <- tumor_subclusters@tumor_subclusters

# Example: get all Fibroblast subclusters
fibro_clusters <- clusters$Fibroblasts
fibro_cells <- unlist(fibro_clusters)
fibro_ids <- names(fibro_cells)   # cell barcodes

all_clusters <- tumor_subclusters@tumor_subclusters
tumor_cells <- unlist(all_clusters)
tumor_ids <- names(tumor_cells)

flatten_clusters <- function(obj) {
  res <- list()
  for (group in names(obj)) {
    sublist <- obj[[group]]
    if (is.list(sublist)) {
      for (sub in names(sublist)) {
        vec <- sublist[[sub]]
        # Case 1: vec is a named integer vector
        if (is.integer(vec) && !is.null(names(vec))) {
          ids <- names(vec)
          res[[length(res) + 1]] <- data.frame(
            cell_id   = ids,
            group     = group,
            subcluster= sub,
            stringsAsFactors = FALSE
          )
        }
        # Case 2: vec is itself a list of named integer vectors
        if (is.list(vec)) {
          for (sub2 in names(vec)) {
            vec2 <- vec[[sub2]]
            if (is.integer(vec2) && !is.null(names(vec2))) {
              ids <- names(vec2)
              res[[length(res) + 1]] <- data.frame(
                cell_id   = ids,
                group     = group,
                subcluster= sub2,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }
  if (length(res) > 0) {
    do.call(rbind, res)
  } else {
    NULL
  }
}

# Apply to the slot
tumor_df <- flatten_clusters(tumor_subclusters@tumor_subclusters)
head(tumor_df)

# CNV scores per cell
cnv_matrix <- tumor_subclusters@expr.data
cnv_scores <- colMeans(cnv_matrix, na.rm = TRUE)

# Align to Seurat cells
mfs1_clean$CNV_score <- cnv_scores[colnames(mfs1_clean)]

# Check cluster metadata
head(mfs1_clean$seurat_clusters)
boxplot(CNV_score ~ seurat_clusters, data = mfs1_clean@meta.data,
        xlab = "Cluster", ylab = "CNV score",
        main = "CNV score distribution across Seurat clusters")

library(ggplot2)
library(ggpubr)

df <- mfs1_clean@meta.data
df$seurat_clusters <- factor(df$seurat_clusters)

# Define some comparisons you want to highlight
comparisons <- list(c("1","2"), c("7","6"), c("5","2"), c("0","5"),c("13","5"))

ggplot(df, aes(x = seurat_clusters, y = CNV_score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "CNV score",
       title = "CNV score distribution across Seurat clusters") +
  # Global Kruskal-Wallis test annotation
  stat_compare_means(method = "kruskal.test", label.y = max(df$CNV_score)*1.05) +
  # Pairwise Wilcoxon tests with significance stars
  stat_compare_means(method = "wilcox.test", comparisons = comparisons,
                     label = "p.signif", p.adjust.method = "bonferroni")
# Zoom out y-axis range
coord_cartesian(ylim = c(0.99, 1.01)) +
  # Add a vivid color palette
  scale_fill_brewer(palette = "Set2")


library(Seurat)
library(ggplot2)

# Representative DE genes per fibroblast-like clusters
genes_to_plot <- c("DCN",     # C0
                   "IFIT3",   # C1
                   "PDGFRA",  # C2
                   "COL5A3",  # C3
                   "TIMP1",   # C4
                   "ISG15",   # C6
                   "IFIT1",   # C7
                   "RUNX2",   # C8
                   "MKI67",   # C9
                   "FGF7")    # C11

# Set cluster identities
Idents(mfs1_clean) <- "seurat_clusters"

# Define fibroblast vs reference clusters
fibro_clusters <- c("0","1","2","3","4","6","7","8","9","11")
ref_clusters   <- c("5","13")
clusters_to_show <- c(fibro_clusters, ref_clusters)

# Loop through genes and plot violin plots
for (gene in genes_to_plot) {
  if (gene %in% rownames(mfs1_clean)) {
    p <- VlnPlot(mfs1_clean,
                 features = gene,
                 group.by = "seurat_clusters",
                 pt.size = 0) +
      ggtitle(paste("Expression of", gene, "Fibroblast vs Reference")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_x_discrete(limits = clusters_to_show)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in Seurat object"))
  }
}

# Look for exact matches
"DCN" %in% rownames(mfs1_clean)
"RUNX2" %in% rownames(mfs1_clean)

# Search for partial matches (sometimes Ensembl IDs or aliases are stored)
grep("DCN", rownames(mfs1_clean), value = TRUE)
grep("RUNX2", rownames(mfs1_clean), value = TRUE)

summary(mfs1_clean[["RNA"]]@data["DCN", ])
summary(mfs1_clean[["RNA"]]@data["RUNX2", ])

rna_data <- GetAssayData(mfs1_clean, assay = "RNA", layer = "data")

genes_to_plot <- c("COL1A1","COL1A2","COL3A1","FN1","ACTA2","PDGFRA","MKI67","ISG15","IFIT3")

for (gene in genes_to_plot) {
  if (gene %in% rownames(rna_data) && sum(rna_data[gene, ]) > 0) {
    p <- VlnPlot(mfs1_clean, features = gene, group.by = "seurat_clusters", pt.size = 0) +
      ggtitle(paste("Expression of", gene, "across clusters")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  } else {
    message(paste("Gene", gene, "not expressed in dataset"))
  }
}



library(Seurat)
library(ggplot2)

# Genes that are expressed in your dataset
genes_to_plot <- c("ACTA2","COL3A1","PDGFRA","IFIT3","FN1","COL1A1","COL1A2","SPARC","VCAN","FBN1",
                   "CD68","CD14","CD3D","CD8A","NKG7")

Idents(mfs1_clean) <- "seurat_clusters"

for (gene in genes_to_plot) {
  if (gene %in% rownames(mfs1_clean)) {
    p <- VlnPlot(mfs1_clean, features = gene, group.by = "seurat_clusters", pt.size = 0) +
      ggtitle(paste("Expression of", gene, "across clusters")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset"))
  }
}



# Run inferCNV for mfs2_clean
infercnv_obj2 <- CreateInfercnvObject(
  raw_counts_matrix = counts_mfs2,
  annotations_file = "mfs2_annotations.txt",
  gene_order_file = "gene_pos_combined.txt",
  ref_group_names = c("T cells","B cells","NK cells","Monocytes","Macrophages","Dendritic cells")
)

infercnv_obj2 <- infercnv::run(
  infercnv_obj2,
  cutoff = 0.1,
  out_dir = "infercnv_output_mfs2",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE
)


library(infercnv)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
library(Seurat)
library(EnsDb.Hsapiens.v86)

# Extract raw counts from Seurat object
counts_mfs2 <- GetAssayData(mfs2_clean, assay = "RNA", layer = "counts")


# Map ENSEMBL IDs to gene symbols if needed
ids2 <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = rownames(counts_mfs2),
                              columns = "SYMBOL",
                              keytype = "ENSEMBL")

# Replace rownames with gene symbols
rownames(counts_mfs2) <- ids2$SYMBOL[match(rownames(counts_mfs2), ids2$ENSEMBL)]

# Build gene position file using EnsDb
gene_ranges <- genes(EnsDb.Hsapiens.v86)
gene_pos <- as.data.frame(gene_ranges)[, c("symbol","seqnames","start","end")]
colnames(gene_pos) <- c("hgnc_symbol","chromosome_name","start_position","end_position")

# Clean chromosome names
gene_pos$chromosome_name <- gsub("^chr", "", gene_pos$chromosome_name)
gene_pos$chromosome_name[gene_pos$chromosome_name == "X"] <- "23"
gene_pos$chromosome_name[gene_pos$chromosome_name == "Y"] <- "24"
gene_pos$chromosome_name[gene_pos$chromosome_name == "MT"] <- "25"
gene_pos$chromosome_name <- as.numeric(gene_pos$chromosome_name)

# Keep only genes present in counts
gene_pos <- gene_pos[gene_pos$hgnc_symbol %in% rownames(counts_mfs2), ]
gene_pos <- gene_pos[!duplicated(gene_pos$hgnc_symbol), ]

# Remove NA values
gene_pos <- gene_pos[!is.na(gene_pos$chromosome_name) &
                       !is.na(gene_pos$start_position) &
                       !is.na(gene_pos$end_position), ]

# Sort by chromosome and start
gene_pos <- gene_pos[order(gene_pos$chromosome_name, gene_pos$start_position), ]

# Reorder counts matrix to match gene order
matched_genes <- intersect(rownames(counts_mfs2), gene_pos$hgnc_symbol)
counts_mfs2 <- counts_mfs2[matched_genes, ]
gene_pos <- gene_pos[match(matched_genes, gene_pos$hgnc_symbol), ]

# Save gene order file
write.table(gene_pos, "gene_pos_mfs2.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Build annotation file
annots2 <- data.frame(
  Cell = colnames(counts_mfs2),
  Cluster = mfs2_clean$SingleR.labels,
  stringsAsFactors = FALSE
)

# Save annotation file
write.table(annots2,
            "mfs2_annotations.txt",
            sep="\t", quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)

# Filter lowly expressed genes
library(Matrix)
counts_mfs2 <- counts_mfs2[rowSums(counts_mfs2) > 30, ]
counts_mfs2 <- as(counts_mfs2, "dgCMatrix")   # sparse format

# Create InferCNV object
infercnv_obj2 <- CreateInfercnvObject(
  raw_counts_matrix = counts_mfs2,
  annotations_file = "mfs2_annotations.txt",
  gene_order_file = "gene_pos_mfs2.txt",
  ref_group_names = c("CD4+ T-cells","CD8+ T-cells","B-cells","NK cells","Monocytes","Macrophages","DC")
)

unlink("infercnv_output_mfs2", recursive = TRUE)

# Run InferCNV
infercnv_obj2 <- infercnv::run(
  infercnv_obj2,
  cutoff = 0.1,
  out_dir = "infercnv_output_mfs2",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE
)



# Load InferCNV tumor subclusters object
# Load InferCNV tumor subclusters object for mfs2
tumor_subclusters2 <- readRDS("infercnv_output_mfs2/15_tumor_subclusters.leiden.infercnv_obj")
class(tumor_subclusters2)
str(tumor_subclusters2, max.level = 1)

# See what slots are available
slotNames(tumor_subclusters2)

# Access the tumor_subclusters slot
clusters2 <- tumor_subclusters2@tumor_subclusters

# Example: get all Fibroblast subclusters
fibro_clusters2 <- clusters2$Fibroblasts
fibro_cells2 <- unlist(fibro_clusters2)
fibro_ids2 <- names(fibro_cells2)   # cell barcodes

all_clusters2 <- tumor_subclusters2@tumor_subclusters
tumor_cells2 <- unlist(all_clusters2)
tumor_ids2 <- names(tumor_cells2)

# Flatten clusters into a dataframe
flatten_clusters <- function(obj) {
  res <- list()
  for (group in names(obj)) {
    sublist <- obj[[group]]
    if (is.list(sublist)) {
      for (sub in names(sublist)) {
        vec <- sublist[[sub]]
        if (is.integer(vec) && !is.null(names(vec))) {
          ids <- names(vec)
          res[[length(res) + 1]] <- data.frame(
            cell_id   = ids,
            group     = group,
            subcluster= sub,
            stringsAsFactors = FALSE
          )
        }
        if (is.list(vec)) {
          for (sub2 in names(vec)) {
            vec2 <- vec[[sub2]]
            if (is.integer(vec2) && !is.null(names(vec2))) {
              ids <- names(vec2)
              res[[length(res) + 1]] <- data.frame(
                cell_id   = ids,
                group     = group,
                subcluster= sub2,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }
  if (length(res) > 0) do.call(rbind, res) else NULL
}

# Apply to the slot
tumor_df2 <- flatten_clusters(tumor_subclusters2@tumor_subclusters)
head(tumor_df2)

# CNV scores per cell
cnv_matrix2 <- tumor_subclusters2@expr.data
cnv_scores2 <- colMeans(cnv_matrix2, na.rm = TRUE)

  # Align to Seurat cells
  mfs2_clean$CNV_score <- cnv_scores2[colnames(mfs2_clean)]
  
  # Quick boxplot
  boxplot(CNV_score ~ seurat_clusters, data = mfs2_clean@meta.data,
          xlab = "Cluster", ylab = "CNV score",
          main = "CNV score distribution across Seurat clusters (mfs2_clean)")
  library(ggplot2)
  library(ggpubr)
  
  df2 <- mfs2_clean@meta.data
  df2$seurat_clusters <- factor(df2$seurat_clusters)
  
  # Define comparisons of interest
  comparisons2 <- list(c("5","2"), c("5","10"), c("5","8"))
  
  ggplot(df2, aes(x = seurat_clusters, y = CNV_score, fill = seurat_clusters)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    labs(x = "Cluster", 
         y = "CNV score (zoomed to 0.99–1.01)", 
         title = "CNV score distribution across Seurat clusters (mfs2_clean)") +
    # Global Kruskal-Wallis test annotation
    stat_compare_means(method = "kruskal.test", 
                       label.y = 1.0095) +   # place annotation near top of new range
    # Pairwise Wilcoxon tests with significance stars
    stat_compare_means(method = "wilcox.test", comparisons = comparisons2,
                       label = "p.signif", p.adjust.method = "bonferroni",
                       label.y = c(1.009, 1.0085, 1.008, 1.0075, 1.007)) + 
    # Zoom out y-axis range
    coord_cartesian(ylim = c(0.99, 1.01)) +
    # Add a vivid color palette
    scale_fill_brewer(palette = "Set2")
  
  library(Seurat)
  library(ggplot2)
str(mfs2_clean)
levels(Idents(mfs2_clean))
length(levels(Idents(mfs2_clean)))
# Access UMAP embeddings
head(Embeddings(mfs2_clean, "umap"))

# Get cluster assignment per cell
head(mfs2_clean$seurat_clusters)

DimPlot(mfs2_clean, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of mfs2_clean (12 clusters)")

DimPlot(mfs1_clean, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of mfs1_clean (17 clusters)")

head(mfs2_clean@meta.data)

mfs1_clean$ClusterCellType <- cluster_map1[as.character(mfs1_clean$seurat_clusters)]

DimPlot(mfs1_clean, reduction = "umap", group.by = "ClusterCellType",
        label = TRUE, repel = TRUE) +
  ggtitle("UMAP of mfs1_clean (17 clusters, annotated)")



library(dplyr)
mfs2_clean <- AddMetaData(
  mfs2_clean,
  metadata = cluster_map2_vec[as.character(mfs2_clean$seurat_clusters)],
  col.name = "ClusterCellType"
)


library(biomaRt)

# Suppose your Seurat object has Ensembl IDs as rownames
ensembl_ids <- rownames(mfs2_clean)

# Map Ensembl → SYMBOL
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = ensembl_ids,
                                  columns = "SYMBOL",
                                  keytype = "ENSEMBL")

head(gene_map)

# Remove NA symbols
gene_map <- gene_map[!is.na(gene_map$SYMBOL), ]

# Deduplicate: keep one symbol per Ensembl ID
gene_map_unique <- gene_map[!duplicated(gene_map$ENSEMBL), ]

# Create mapping vector
symbol_map <- setNames(gene_map_unique$SYMBOL, gene_map_unique$ENSEMBL)

# Replace rownames with SYMBOLs where available
rownames(mfs2_clean) <- ifelse(rownames(mfs2_clean) %in% names(symbol_map),
                               symbol_map[rownames(mfs2_clean)],
                               rownames(mfs2_clean))

VlnPlot(mfs2_clean, features = c("PDGFRB", "COL1A2", "PDGFRA", "NOTCH3", "FAP", "EGFR"),
        group.by = "seurat_clusters",
        pt.size = 0.1)


library(Seurat)
library(ggplot2)
library(ggpubr)

library(org.Hs.eg.db)
# Define fibroblast clusters and immune reference clusters
fibro_clusters <- c("7","9","10")
immune_clusters <- c("0","1","2","4","5","6")

# Representative fibroblast markers and oncogenes
fibroblast_genes <- c("PDGFRB", "COL1A2", "PDGFRA")
oncogenes <- c("NOTCH3", "FAP", "EGFR")

# Combine into one list
genes_to_plot <- c(fibroblast_genes, oncogenes)

# Make sure genes exist in dataset
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(mfs2_clean)]

# Violin plots
VlnPlot(mfs2_clean, features = genes_to_plot,
        group.by = "seurat_clusters",
        pt.size = 0.1,
        cols = c("grey70","grey80","grey90","grey60","grey50","grey40",
                 "firebrick","darkorange","steelblue","forestgreen","purple","gold")) +
  ggtitle("Fibroblast markers and oncogenes across clusters (mfs2_clean)")




library(dplyr)

# For mfs1_clean
df1 <- data.frame(
  Cluster   = mfs1_clean$seurat_clusters,
  CellType  = mfs1_clean$SingleR.labels
)

summary1 <- df1 %>%
  group_by(Cluster, CellType) %>%
  summarise(Cells = n(), .groups = "drop") %>%
  mutate(Percent = round(Cells / ncol(mfs1_clean) * 100, 2))

summary1

# For mfs2_clean
df2 <- data.frame(
  Cluster   = mfs2_clean$seurat_clusters,
  CellType  = mfs2_clean$SingleR.labels
)

summary2 <- df2 %>%
  group_by(Cluster, CellType) %>%
  summarise(Cells = n(), .groups = "drop") %>%
  mutate(Percent = round(Cells / ncol(mfs2_clean) * 100, 2))

summary2


library(dplyr)
library(Seurat)


library(ggplot2)


summary1_df <- summary1 %>% ungroup() %>% as.data.frame()
str(summary1_df)
colnames(summary1_df)
library(ggplot2)

# Add a label column for plotting
summary1_df$Label <- paste0(summary1_df$SingleR.labels,
                            "\n", summary1_df$Cells,
                            " (", summary1_df$Percent, "%)")

library(ggplot2)

ggplot(summary1_df, aes(x = reorder(SingleR.labels, -Cells), y = Cells, fill = SingleR.labels)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  ggtitle("mfs1_clean Cell Type Composition") +
  xlab("Cell Type") +
  ylab("Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_text(aes(label = paste0(Cells, " (", Percent, "%)")),
            vjust = -0.3, size = 3)


summary2_df <- summary2 %>% ungroup() %>% as.data.frame()
str(summary2_df)
colnames(summary2_df)

library(dplyr)
library(ggplot2)

# Collapse across clusters: sum counts and recompute percentages
summary2_total <- summary2_df %>%
  group_by(CellType) %>%
  summarise(Cells = sum(Cells), .groups = "drop") %>%
  mutate(Percent = round(Cells / sum(Cells) * 100, 2))

# Plot histogram
ggplot(summary2_total, aes(x = reorder(CellType, -Cells), y = Cells, fill = CellType)) +
  geom_col(color = "black") +
  theme_minimal() +
  ggtitle("mfs2_clean Cell Type Composition (Total)") +
  xlab("Cell Type") +
  ylab("Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_text(aes(label = paste0(Cells, " (", Percent, "%)")),
            vjust = -0.3, size = 3)



library(Seurat)
library(ggplot2)
library(dplyr)

# 1. Plot colored by SingleR cell type
p <- DimPlot(mfs2_clean, group.by = "SingleR.labels", label = FALSE) +
  theme(legend.position = "right")

# 2. Extract UMAP coordinates and cluster IDs
umap_data <- Embeddings(mfs2_clean, "umap")
cluster_ids <- mfs2_clean$seurat_clusters

# 3. Compute cluster centers
centers <- aggregate(umap_data, by = list(cluster = cluster_ids), FUN = mean)
colnames(centers)[2:3] <- c("UMAP_1", "UMAP_2")

# 4. Overlay cluster numbers at centers
p <- p + geom_text(data = centers,
                   aes(x = UMAP_1, y = UMAP_2, label = cluster),
                   size = 6, fontface = "bold", color = "black")

# 5. Display
p

