
library(tidyverse)
library(vcfR)
library(SNPfiltR)

# New library for PCA
library(LEA)

setwd("~/Projects/eco_genomics/population_genomics/")
vcf <- read.vcfR("outputs/vcf.filtered.vcf.gz")

# Thinning SNPs for linkage disequilibrium- thinning spatially
# Necessary step for PCA and Admixture (assumption of independence)

vcf.thin <-distance_thin(vcf, min.distance = 500)

# Most SNPs are stacked (near eachother)-- 
# Typical of GBS because we sequence fragments

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)

# Again, make sure dimensions of meta matches our new vcf dimensions
meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[, -1]) , ]


write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

# We have to uncompress this file for LEA
# this will make too big of a file for github- hide outside of repo
# system means its as if you're on command line in BASH
system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

# PCA with thinned data
# Call LEA then PCA to get correct pca function

geno <- vcf2geno(input.file = "/gpfs1/home/a/r/arodger/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.vcf.geno")
CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.vcf.geno", scale=TRUE)

# Using base R
# Porjections is first 2 pc axes
plot(CentPCA$projections,
     col=as.factor(meta2$region))
legend("bottomright", legend=as.factor(unique(meta2$region)),
       fill=as.factor(unique(meta2$region)))