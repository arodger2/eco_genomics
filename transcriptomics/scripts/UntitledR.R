### Code for analyzing RNAseq data using DESeq2

# load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
setwd("~/Projects/eco_genomics/transcriptomics/")

# Import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header=TRUE, row.names = 1)
# DESeq2 doesn't like decimals but Salmon creates decimals
countsTableRound <- round(countsTable)

# Experimental details
conds <-read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                   header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

dds <- DESeqDataSetFromMatrix (countData= countsTableRound, colData=conds,
                              design= ~DevTemp+FinalTemp)
