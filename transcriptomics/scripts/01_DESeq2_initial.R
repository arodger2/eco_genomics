### Code for analyzing RNAseq data using DESeq2

# load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
options(bitmapType = "cairo")

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

##########################
#
# Explore counts matrix
#
###########################

# Let's see how many reads we have from each sample

colSums(countsTableRound)
mean (colSums(countsTableRound))
# 20 million reads before filtering is the goal

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names=0.5, las =2, ylim=c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col="blue4", lwd=2)


# The average number of case counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # 3244.739
median(rowSums(countsTableRound)) #64


apply(countsTableRound, 2, mean) # across columns
apply(countsTableRound, 1, mean) # across rows (not useful), gives a sense of variation in sequencing effort across samples


#####################
#
# DESeq2
#
#####################

dds <-DESeqDataSetFromMatrix(countData=countsTableRound, colData=conds,
                             design= ~DevTemp + FinalTemp)

dim(dds)

# Filtering based on average number of transcripts 
# If there aren't more than 10 reads for at least 15 out of 21 samples, we are ignoring
dds<- dds[rowSums(counts(dds) >= 10)>= 15,]
nrow(dds) # 35,527 transcripts (way more transcripts than genes because of splicing)

# Run the DESeq model to test for global differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
# [1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28" 
# [4] "FinalTemp_BASE_vs_A28"

# Visualize out global gene expression patterns using PCA!
# First, transform the data for plotting using variance stabilization (function in DESeq2-normalizes)

vsd <- vst(dds, blind = F)

pcaData <- plotPCA(vsd, intgroup=c("DevTemp","FinalTemp"), returnData=T)
percentVar <- round(100*attr(pcaData, "percentVar"))

final_temp_colors <- c("BASE"="grey", "A28"="hotpink", "A33"="red")
shapes_choose <- c("D18"= 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color=FinalTemp, shape=DevTemp))+
  geom_point(size=5)+
  scale_shape_manual(values=shapes_choose)+
  scale_color_manual(values=final_temp_colors)+
  labs(x= paste0("PC1: ", percentVar[1], "%"), 
       y= paste0("PC2: ", percentVar[2], "%")) +
  theme_bw(base_size = 16)

p
  