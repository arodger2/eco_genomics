

library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)


vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz",
                    type="vcf")
vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")
#vcfR <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf.filtered.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]), ]

# componentwise test for selection separately on each axis
# rare freq alleles are represented as minor allele frequency
# min.maf in PCA avoids really rare allels
# PCA general assumes no linkage disequelibrium, but PCAdapt can break this assumption
pcadapt.pca <- pcadapt(vcf,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))
summary(pcadapt.pca)
plot(pcadapt.pca, options="stat.distribution",
     pop=meta2$region,
    i=1, j=2,
    K=2)

# Get manhattan plot into similar format to before (qqman)
View(head(vcfR@fix))

# making a data frame
# again separating by chromosome and assigning chromosome number

vcfR.fix <- as.data.frame(vcfR@fix[,1:2])
chr.main <- unique(vcfR.fix$CHROM)[1:8]
chrnum <- as.data.frame(cbind(chr.main, seq (1,8,1)))

# Getting p values
Pval <- pcadapt.pca$pvalues
pcadapt.MHplot <- cbind(vcfR.fix, Pval)


pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, join_by(chr.main==CHROM))
pcadapt.MHplot <- pcadapt.MHplot %>%
  mutate(SNP=paste0(chr.main,"_",POS))

pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1 = as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2 = as.numeric(pcadapt.MHplot[,5])

pcadapt.MHplot <- pcadapt.MHplot %>% 
  drop_na(pPC1)

manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log 10 p-value",
          genomewideline = F,
          main="PCAdapt genome scan for selection (PC1)")

View(pcadapt.MHplot %>% 
       filter(pPC1<quantile(pcadapt.MHplot$pPC1, 0.001)) %>%
       select(chr.main, POS, pPC1))

