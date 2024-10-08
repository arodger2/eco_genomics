library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

# Fixed section is all SNPs (alt compared to reference genome)
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

# Annotation file
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

#
DP <- extract.gt(vcf, element = "DP", as.numeric = T)

# Starting filtering
DP[DP==0] <- NA
quantile(DP, na.rm = T)

# Begin filtering by depth
library(SNPfiltR)

# filter low depths
vcf.filt <- hard_filter(vcf, depth= 3)

# filter high depths (2x the mean)
vcf.filt <- max_depth(vcf.filt, maxdepth = 60)

meta <- read.csv("metadata/meta4vcf.csv", header=T)

# subsetting meta
meta2 <- meta[,c(1,4)]

names(meta2) <- c("id","pop")
meta2$id =as.factor(meta2$id)
meta2$pop =as.factor(meta2$pop)

# Looking at individual missingness, and setting cutoff
sixtyfive <- missing_by_sample(vcf.filt,
                            popmap = meta2,
                            cutoff=0.65)
eightyfive <-missing_by_sample(vcf.filt,
                                popmap = meta2,
                                cutoff=0.85)

sixtyfive.filt <- filter_biallelic(sixtyfive)
eightyfive.filt <- filter_biallelic(eightyfive)

# sets the minimum number of times to see SNP
sixtyfive.filt <- min_mac(sixtyfive.filt, min.mac = 1)
eightyfive.filt <- min_mac(eightyfive.filt, min.mac = 1)

# Snp-wise missing alleles
sixtyfive.filt.SNPMiss <- missing_by_snp(sixtyfive.filt, cutoff= 0.5)
eightyfive.filt.SNPMiss <- missing_by_snp(eightyfive.filt, cutoff= 0.5)

write.vcf(sixtyfive.filt.SNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/sixtyfive.filt.SNPMiss.gz")

write.vcf(eightyfive.filt.SNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/eightyfive.filt.SNPMiss.gz")

#### Script 2

library(vcfR)
library(tidyverse)
library(qqman)

X11.options(type="cairo")

##### begining with code for 65% filter
sixtyfive.vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/sixtyfive.filt.SNPMiss.gz")

meta265 <-meta2[meta2$id %in% colnames(sixtyfive.vcf@gt[,-1]),]

# calculate diversity stats using the genetic_diff fxn in vcfR
vcf65.div <- genetic_diff(sixtyfive.vcf,
                        pops=as.factor(meta265$region),
                        method="nei")
# Could also be interesting to look at other groupings than region?

# Taking just first 8 entries (8 main chromsomes)
chr.main <- unique(vcf65.div$CHROM) [1:8]

# Tell it which chromosome numbers
chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

# merge diversity stats with chrnum- assign to chromosome numbers

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
# Interpreting manhattan plot- these individuals generally share a large amount of SNPs
# Write out file

write.csv(vcf65.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion65.csv",
          quote=F,
          row.names=F)

##### Repeating for 85% missingness
eightyfive.vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/eightyfive.filt.SNPMiss.gz")

meta285 <-meta[meta$id %in% colnames(eightyfive.vcf@gt[,-1]),]

# calculate diversity stats using the genetic_diff fxn in vcfR
vcf85.div <- genetic_diff(eightyfive.vcf,
                          pops=as.factor(meta285$region),
                          method="nei")
# Could also be interesting to look at other groupings than region?

# Taking just first 8 entries (8 main chromsomes)
chr.main <- unique(vcf85.div$CHROM) [1:8]

# Tell it which chromosome numbers
chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

# merge diversity stats with chrnum- assign to chromosome numbers

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
# Interpreting manhattan plot- these individuals generally share a large amount of SNPs
# Write out file

write.csv(vcf85.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion85.csv",
          quote=F,
          row.names=F)

# Somewhat different, but may not be the most valuable graphs^

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

# To me these plots look pretty similar in terms of trends
# value means generic value, coloring by column name
# alpha is transparency

# Saving the graph
# ggsave("Histogram_GenomeDiversity_byRegion.pdf",
       # path="~/Projects/eco_genomics/population_genomics/figures/")

# Making a summary table showing averages of the data we just plotted
vcf65.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

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

vcf65.thin <-distance_thin(vcf65, min.distance = 500)
vcf85.thin <- distance_thin(vcf85, min.distance = 500)

# Again, make sure dimensions of meta matches our new vcf dimensions
meta265 <- meta[meta$id %in% colnames(vcf65@gt[, -1]) , ]
meta285 <- meta[meta$id %in% colnames(vcf85@gt[, -1]) , ]

write.vcf(vcf65.thin, "outputs/vcf65_final.filtered.thinned.vcf.gz")
write.vcf(vcf85.thin, "outputs/vcf85_final.filtered.thinned.vcf.gz")

# We have to uncompress this file for LEA
# this will make too big of a file for github- hide outside of repo
# system means its as if you're on command line in BASH
system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf65_final.filtered.thinned.vcf.gz > ~/vcf65_final.filtered.thinned.vcf")

# PCA with thinned data
# Call LEA then PCA to get correct pca function

geno65 <- vcf2geno(input.file = "/gpfs1/home/a/r/arodger/vcf65_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf65_final.filtered.thinned.vcf.geno")
CentPCA65 <- LEA::pca("outputs/vcf65_final.filtered.thinned.vcf.geno", scale=TRUE)

# To reload a previous PCA
CentPCA65 <- load.pcaProject("vcf65_final.filtered.thinned.vcf.pcaProject")

ggplot(as.data.frame(CentPCA65$projections),
       aes(x=V1, y=V2, color=meta265$region, shape=meta265$continent)) +
  geom_point() +
  labs(title = "Centaurea genetic PCA (65% missingness)", x="PC1",y="PC2", color="Region", shape="Continent")

ggsave("figures/CentPCA65_PC1vPC2.pdf")
# Repeat with 85% missingness
system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf85_final.filtered.thinned.vcf.gz > ~/vcf85_final.filtered.thinned.vcf")

geno85 <- vcf2geno(input.file = "/gpfs1/home/a/r/arodger/vcf85_final.filtered.thinned.vcf",
                   output.file = "outputs/vcf85_final.filtered.thinned.vcf.geno")
CentPCA85 <- LEA::pca("outputs/vcf85_final.filtered.thinned.vcf.geno", scale=TRUE)

# To reload a previous PCA
CentPCA85 <- load.pcaProject("vcf85_final.filtered.thinned.vcf.pcaProject")

ggplot(as.data.frame(CentPCA85$projections),
       aes(x=V1, y=V2, color=meta285$region, shape=meta285$continent)) +
  geom_point() +
  labs(title = "Centaurea genetic PCA (85% missingness)", x="PC1",y="PC2", color="Region", shape="Continent")

ggsave("figures/CentPCA85_PC1vPC2.pdf")

### Selection

library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)
