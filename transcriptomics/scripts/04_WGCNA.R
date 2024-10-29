### 04_WGCNA.R
### Script for analyzing and visualizing gene correlation networks

library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringsAsFactors=FALSE);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)

options(bitmapType = "cairo")

setwd("~/Projects/eco_genomics/transcriptomics/")

# Step 1: Import counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header=TRUE, row.names = 1)
countsTableRound <- round(countsTable)
conds <-read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                   header=TRUE, stringsAsFactors = TRUE, row.names=1)
traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header=T, row.names=1)

# Filter the matrix to just base data (because those are the data for which we have traits measured)

filtered_count_matrix_BASEonly <-countsTable[ ,conds$FinalTemp== "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp== "BASE", ]

rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)


# Step 2: Detecting outliers
# detect outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes)
# 
#FALSE  TRUE 
#37235 82203 

table(gsg$goodSamples)

#TRUE 
#7 

# Filter out bad genes
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes==TRUE, ]

dim(data_WGCNA)

# use clustering with a tree dendrogram to identify outlier samples
htree <- hclust(dist(t(data_WGCNA)), method="average")
plot(htree)

#Leaving in the outlier for now, because all samples were "good"

# PCA- outlier detection method
pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
# transform into a data frame
pca_data <- as.data.frame(pca_data)

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits=2)

ggplot(pca_data, aes(PC1,PC2))+
  geom_point()+
  geom_text(label=rownames(pca_data))+
  labs (x=paste0('PC1: ', pca.var.percent[1], "%"),
        y=paste0('PC2: ', pca.var.percent[2], "%"))

# Step 3: Normalization
colData <- row.names(filtered_sample_metadata_BASEonly)

# Run DESeq without any model defined
dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA, colData = filtered_sample_metadata_BASEonly,
                                    design=~1)# no specified groups

dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >=6,]

nrow(dds_WGCNA_75) #filtered down to 29559 transcripts

dds_norm <- vst(dds_WGCNA_75) # perform variance stabilization

# get and save normalized counts to use below
norm.counts <-assay(dds_norm) %>% 
  t()

# Step 4: Network Construction!

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from=12, to=50, by=2))

# Call the network topology analysis function (takes a couple min to run)
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose=5) # only genes positively correlates, not anticorrelated

sft.data <- sft$fitIndices

# plot to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label= Power))+
  geom_point()+
  geom_text(nudge_y=0.1)+
  geom_hline(yintercept = 0.8, color='red')+
  labs(x="Power",y="Scale free topology model fit, signed, R^2")+
  theme_classic()

a1
# Scale free topology is a type of network, power is the strength of the correlation
# We want to know how well certain powers fit a scale free topology

a2<- ggplot(sft.data, aes(Power, mean.k., label= Power))+
  geom_point()+
  geom_text(nudge_y=0.1)+
  geom_hline(yintercept = 0.8, color='red')+
  labs(x="Power",y="Mean connectivity")+
  theme_classic()

grid.arrange(a1, a2, nrow=2)

# Choose a power of 26
# power is strength/weakness of correlation between genes
# want to pick a power that has a more "biological topology" while not getting too low

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor

# this sets the temp cor function to use WGCNA's correlation function

norm.counts[] <- sapply(norm.counts, as.numeric)

# The command below creates a network and identifies modules based on the parameters
# that we chose
bwnet26 <- blockwiseModules(norm.counts,
                            maxBlockSize = 30000,
                            TOMType = "signed",
                            power = soft_power,
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            randomSeed = 1234,
                            verbose=3) # TOMtype signed means only focus on positive correlations

cor <- temp_cor # This resets the cor function to basse R's cor (instead of using WGCNA's)

saveRDS(bwnet26, file = "outputs/bwnet26.rds")
# To load the bwnet file in later, use:
# bwnet26 <- readRDS("outputs/bwnet26.rds") #


# STEP 5: Explore Module Eigengenes
module_eigengenes <- bwnet26$MEs
head(module_eigengenes)
dim(module_eigengenes)

# get the number of genes for each module using table function (each module is named for a color)
table(bwnet$colors)

# Plot the dendrogram and the module colors
plotDendroAndColors(bwnet26$dendrograms[[1]], cbind(bwnet26$unmergedColors, bwnet26$colors),
                    c("unmerged", "merged"),
                    dendroLabels= FALSE,
                    addGuide= TRUE,
                    hang= 0.03,
                    guidehang= 0.05)


# STEP 6: Correlation of modules with traits
# Define the numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Test for a correlation between module eigen genes and trait data
module.trait.corr <- cor(module_eigengenes, traitData, use = "p")

# Calculate p values for each correlation
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Visualise module trait association as a heatmap!
heatmap.data <- merge(module_eigengenes, traitData, by = "row.names")
head(heatmap.data)

# address error of row.names not being numeric
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var="Row.names")

names(heatmap.data)

# Make pretty heatmap of correlations
CorLevelPlot(heatmap.data,
             x= names(heatmap.data)[42:44], # these values may need to change based on 
             y = names(heatmap.data) [1:41],   # number of eigengenes
             col = c("blue2", "skyblue", "white", "pink", "red"))




            




