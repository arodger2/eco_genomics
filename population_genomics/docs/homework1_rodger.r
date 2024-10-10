# Loading in library
################# 
library(vcfR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")


# Loading files
###############
vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")


# Starting filtering
#####
DP <- extract.gt(vcf, element = "DP", as.numeric = T)
DP[DP==0] <- NA
quantile(DP, na.rm = T)

library(SNPfiltR)

# filter low depths
vcf.filt <- hard_filter(vcf, depth= 3)

# filter high depths (2x the mean)
vcf.filt <- max_depth(vcf.filt, maxdepth = 60)

# load in metadata (individual IDs, sites, populations, etc.)
meta <- read.csv("metadata/meta4vcf.csv", header=T) 
#meta begins as 629 observations of 8 variables

# to take only the columns "id" and "pop"...
# meta2 <- meta[,c(1,4)]
# names(meta2) <- c("id","pop")
# meta2$id =as.factor(meta2$id)
# meta2$pop =as.factor(meta2$pop)

# Setting cutoff for individual-wise missing-ness
sixtyfive <- missing_by_sample(vcf.filt,
                            popmap = meta2,
                            cutoff=0.65)
eightyfive <-missing_by_sample(vcf.filt,
                                popmap = meta2,
                                cutoff=0.85)
# Filtering out any alleles which aren't biallelic
sixtyfive.filt <- filter_biallelic(sixtyfive)
eightyfive.filt <- filter_biallelic(eightyfive)

# sets the minimum number of times to see SNP to 1
sixtyfive.filt <- min_mac(sixtyfive.filt, min.mac = 1)
eightyfive.filt <- min_mac(eightyfive.filt, min.mac = 1)

# Snp-wise missing alleles
sixtyfive.filt.SNPMiss <- missing_by_snp(sixtyfive.filt, cutoff= 0.5)
eightyfive.filt.SNPMiss <- missing_by_snp(eightyfive.filt, cutoff= 0.5)

# Write out the vcfs
write.vcf(sixtyfive.filt.SNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/sixtyfive.filt.SNPMiss.gz")
write.vcf(eightyfive.filt.SNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/eightyfive.filt.SNPMiss.gz")

# Script 2
######
library(vcfR)
library(tidyverse)
library(qqman)

X11.options(type="cairo")

##### begining with code for 65% filter
sixtyfive.vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/sixtyfive.filt.SNPMiss.gz")
meta65 <-meta[meta$id %in% colnames(sixtyfive.vcf@gt[,-1]),]
vcf65.div <- genetic_diff(sixtyfive.vcf,
                        pops=as.factor(meta65$region),
                        method="nei")

# Assigning SNPs to chromosomes
chr.main <- unique(vcf65.div$CHROM) [1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

vcf65.div.MHplot <- left_join(chrnum, vcf65.div, join_by(chr.main==CHROM))
vcf65.div.MHplot <- vcf65.div.MHplot %>%
  filter(Gst>0) %>% 
  mutate(SNP=paste0(chr.main, "_", POS))

# Make sure numbers are viewed as numbers
vcf65.div.MHplot$V2 = as.numeric(vcf65.div.MHplot$V2)
vcf65.div.MHplot$POS = as.numeric(vcf65.div.MHplot$POS)

# Manhattan plot with qqman
manhattan(vcf65.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          main="Diversity among regions (65% missingness filter)",
          suggestiveline = quantile(vcf65.div.MHplot$Gst, 0.50))
# These individuals generally share a large amount of SNPs

##### Repeating for 85% missingness
eightyfive.vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/eightyfive.filt.SNPMiss.gz")

meta85 <-meta[meta$id %in% colnames(eightyfive.vcf@gt[,-1]),]

# calculate diversity stats using the genetic_diff fxn in vcfR
vcf85.div <- genetic_diff(eightyfive.vcf,
                          pops=as.factor(meta85$region),
                          method="nei")

# Assigning to chromosomes
chr.main <- unique(vcf85.div$CHROM) [1:8]
chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))
vcf85.div.MHplot <- left_join(chrnum, vcf85.div, join_by(chr.main==CHROM))

vcf85.div.MHplot <- vcf85.div.MHplot %>%
  filter(Gst>0) %>% 
  mutate(SNP=paste0(chr.main, "_", POS))

# Make sure numbers are viewed as numbers
vcf85.div.MHplot$V2 = as.numeric(vcf85.div.MHplot$V2)
vcf85.div.MHplot$POS = as.numeric(vcf85.div.MHplot$POS)

# Manhattan plot with qqman
manhattan(vcf85.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          main="Diversity among regions (85% missingness filter)",
          suggestiveline = quantile(vcf85.div.MHplot$Gst, 0.50))

# Looking at Hs values
# Hs values are stored in columns 4:9
# Density plot
# Starting with 65
vcf65.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x=value, fill=name))+
  geom_histogram(position="identity",alpha=1, bins=50)+
  labs(title="Genome-wide expected heterozygosity (Hs)--65% missingness", 
       fill="Regions",
       x="Gene diversity within regions",
       y="Counts of SNPs")

# Now 85
vcf85.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x=value, fill=name))+
  geom_histogram(position="identity",alpha=1, bins=50)+
  labs(title="Genome-wide expected heterozygosity (Hs)--85% missingness", 
       fill="Regions",
       x="Gene diversity within regions",
       y="Counts of SNPs")

# Making a summary table showing averages of the data we just plotted
diversitysummary65 <- vcf65.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

lowerlim65=diversitysummary65$avgHs-diversitysummary65$StdDev_Hs
upperlim65=diversitysummary65$avgHs+diversitysummary65$StdDev_Hs

ggplot(diversitysummary65, aes(x=name, y=avgHs, color=name))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lowerlim65, ymax=upperlim65))+
  labs(title="Genome-wide expected heterozygosity (Hs)--65% missingness", 
       fill="Regions",
       x="Region",
       y="Average Expected Heterozygosity")

diversitysummary85 <- vcf85.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

lowerlim85=diversitysummary85$avgHs-diversitysummary85$StdDev_Hs
upperlim85=diversitysummary85$avgHs+diversitysummary85$StdDev_Hs

ggplot(diversitysummary85, aes(x=name, y=avgHs, color=name))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lowerlim85, ymax=upperlim85))+
  labs(title="Genome-wide expected heterozygosity (Hs)--85% missingness", 
       fill="Regions",
       x="Region",
       y="Average Expected Heterozygosity")
vcf85.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n())
  

##### Script 3
library(LEA)
setwd("~/Projects/eco_genomics/population_genomics/")
vcf65 <- read.vcfR("outputs/sixtyfive.filt.SNPMiss.gz")
vcf85 <- read.vcfR("outputs/eightyfive.filt.SNPMiss.gz")


# Thinning SNPs which are closer than 500 bp apart because of PCA's assumption of independence
vcf65.thin <-distance_thin(vcf65, min.distance = 500)
vcf85.thin <- distance_thin(vcf85, min.distance = 500)

# Estimating individuals and SNP loci after all filtering
individuals65 <- ncol(vcf65.thin@gt [, -1])
Snps65 <- nrow (vcf65.thin@gt)

individuals85 <- ncol(vcf85.thin@gt [, -1])
Snps85 <- nrow (vcf85.thin@gt)

Estimates <- matrix(c(individuals65, Snps65, individuals85, Snps85), nrow=2, ncol=2) 
colnames(Estimates) <- c("65% missing data","85% missing data")
rownames(Estimates) <- c("Individuals","SNPs")
Estimatestab <- as.table(Estimates)
Estimatestab

# Make sure dimensions of meta matches our new vcf dimensions
meta65 <- meta[meta$id %in% colnames(vcf65.thin@gt[, -1]) , ]
meta85 <- meta[meta$id %in% colnames(vcf85.thin@gt[, -1]) , ]

write.vcf(vcf65.thin, "outputs/vcf65_final.filtered.thinned.vcf.gz")
write.vcf(vcf85.thin, "outputs/vcf85_final.filtered.thinned.vcf.gz")

system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf65_final.filtered.thinned.vcf.gz > ~/vcf65_final.filtered.thinned.vcf")

# PCA with thinned data
geno65 <- vcf2geno(input.file = "/gpfs1/home/a/r/arodger/vcf65_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf65_final.filtered.thinned.vcf.geno")
CentPCA65 <- LEA::pca("outputs/vcf65_final.filtered.thinned.vcf.geno", scale=TRUE)
# Calculating percentage of variation explained by each PC
CentPCA65$eigenvalues[1]/sum(CentPCA65$eigenvalues)
CentPCA65$eigenvalues[2]/sum(CentPCA65$eigenvalues)

ggplot(as.data.frame(CentPCA65$projections),
       aes(x=V1, y=V2, color=meta65$region, shape=meta65$continent)) +
  geom_point() +
  labs(title = "Centaurea genetic PCA (65% missingness)", x="PC1 (2.4%)",y="PC2 (1.1%)", color="Region", shape="Continent")
ggsave("figures/CentPCA65_PC1vPC2.pdf")

# Repeat with 85% missingness
system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf85_final.filtered.thinned.vcf.gz > ~/vcf85_final.filtered.thinned.vcf")
geno85 <- vcf2geno(input.file = "/gpfs1/home/a/r/arodger/vcf85_final.filtered.thinned.vcf",
                   output.file = "outputs/vcf85_final.filtered.thinned.vcf.geno")
CentPCA85 <- LEA::pca("outputs/vcf85_final.filtered.thinned.vcf.geno", scale=TRUE)
CentPCA85$eigenvalues[1]/sum(CentPCA85$eigenvalues)
CentPCA85$eigenvalues[2]/sum(CentPCA85$eigenvalues)

ggplot(as.data.frame(CentPCA85$projections),
       aes(x=V1, y=V2, color=meta85$region, shape=meta85$continent)) +
  geom_point() +
  labs(title = "Centaurea genetic PCA (85% missingness)", x="PC1 (2.2%)",y="PC2(1.0%)", color="Region", shape="Continent")
ggsave("figures/CentPCA85_PC1vPC2.pdf")

### Selection
library(pcadapt)
vcf65 <- read.pcadapt("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf65.gz",
                    type="vcf")
pcadapt.pca65 <- pcadapt(vcf65,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))

vcfR65 <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf65.gz")
vcfR.fix65 <- as.data.frame(vcfR65@fix[,1:2])
chr.main65 <- unique(vcfR.fix65$CHROM)[1:8]
chrnum65 <- as.data.frame(cbind(chr.main65, seq (1,8,1)))

# Getting p values
Pval65 <- pcadapt.pca65$pvalues
pcadapt.MHplot65 <- cbind(vcfR.fix65, Pval65)
pcadapt.MHplot65 <- left_join(chrnum65, pcadapt.MHplot65, join_by(chr.main65==CHROM))
pcadapt.MHplot65 <- pcadapt.MHplot65 %>%
  mutate(SNP=paste0(chr.main65,"_",POS))

pcadapt.MHplot65$V2 = as.numeric(pcadapt.MHplot65$V2)
pcadapt.MHplot65$POS = as.numeric(pcadapt.MHplot65$POS)

pcadapt.MHplot65$pPC1 = as.numeric(pcadapt.MHplot65[,4])
pcadapt.MHplot65$pPC2 = as.numeric(pcadapt.MHplot65[,5])

pcadapt.MHplot65 <- pcadapt.MHplot65 %>% 
  drop_na(pPC1)

manhattan(pcadapt.MHplot65,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log 10 p-value",
          genomewideline = F,
          main="PCAdapt genome scan for selection (PC1)- 65%")


# Repeating for 85%
vcf85 <- read.pcadapt("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf85.gz",
                      type="vcf")
pcadapt.pca85 <- pcadapt(vcf85,
                         K=2,
                         method="componentwise",
                         min.maf=0.01,
                         LD.clumping = list(size=500, thr=0.2))

vcfR85 <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf85.gz")
vcfR.fix85 <- as.data.frame(vcfR85@fix[,1:2])
chr.main85 <- unique(vcfR.fix85$CHROM)[1:8]
chrnum85 <- as.data.frame(cbind(chr.main85, seq (1,8,1)))

# Getting p values
Pval85 <- pcadapt.pca85$pvalues
pcadapt.MHplot85 <- cbind(vcfR.fix85, Pval85)
pcadapt.MHplot85 <- left_join(chrnum85, pcadapt.MHplot85, join_by(chr.main85==CHROM))
pcadapt.MHplot85 <- pcadapt.MHplot85 %>%
  mutate(SNP=paste0(chr.main85,"_",POS))

pcadapt.MHplot85$V2 = as.numeric(pcadapt.MHplot85$V2)
pcadapt.MHplot85$POS = as.numeric(pcadapt.MHplot85$POS)

pcadapt.MHplot85$pPC1 = as.numeric(pcadapt.MHplot85[,4])
pcadapt.MHplot85$pPC2 = as.numeric(pcadapt.MHplot85[,5])

pcadapt.MHplot85 <- pcadapt.MHplot85 %>% 
  drop_na(pPC1)

manhattan(pcadapt.MHplot85,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log 10 p-value",
          genomewideline = F,
          main="PCAdapt genome scan for selection (PC1)- 85%")