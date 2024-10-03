library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
list.files("variants/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

# Fixed section is all SNPs (alt compared to reference genome)

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

# Annotation file

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

# Using chromosome 1 as an example

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1)

# Saving out files to the repo
pdf(file="~/Projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()

# 9/17/24
#
DP <- extract.gt(vcf, element = "DP", as.numeric = T)
dim(DP)
DP[1:5,1:10]


# Starting filtering
DP[DP==0] <- NA
quantile(DP, na.rm = T)

# Visual tool of depth and missingness

heatmap.bp(DP, rlabels=F, clabels=F)

# Begin filtering by depth

library(SNPfiltR)

# filter low depths
vcf.filt <- hard_filter(vcf, depth= 3) # could explore other depths, 5, 10

# filter high depths (2xthe mean)
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

DP2sixtyfive <-extract.gt(sixtyfive.filt.SNPMiss,
                 element="DP",
                 as.numeric = T)

DP2eightyfive <-extract.gt(eightyfive.filt.SNPMiss,
                          element="DP",
                          as.numeric = T)

heatmap.bp(DP2sixtyfive, rlabels = F, clabels = F)
heatmap.bp(DP2eightyfive, rlabels = F, clabels =F)

write.vcf(sixtyfive.filt.SNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/sixtyfive.filt.SNPMiss.gz")

write.vcf(eightyfive.filt.SNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/eightyfive.filt.SNPMiss.gz")


str(sixtyfive.filt.SNPMiss)
str(eightyfive.filt.SNPMiss)



#### Script 2

library(vcfR)
library(tidyverse)
library(qqman)

X11.options(type="cairo")

# begining with code for 65% filter
sixtyfive.vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/sixtyfive.filt.SNPMiss.gz")

meta2 <-meta[meta$id %in% colnames(sixtyfive.vcf@gt[,-1]),]

# calculate diversity stats using the genetic_diff fxn in vcfR
vcf65.div <- genetic_diff(sixtyfive.vcf,
                        pops=as.factor(meta2$region),
                        method="nei")
# Could also be interesting to look at other groupings than region?

# Taking just first 8 entries (8 main chromsomes)
chr.main <- unique(vcf65.div$CHROM) [1:8]

# Tell it which chromosome numbers

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))
chrnum

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
          suggestiveline = quantile(vcf65.div.MHplot$Gst, 0.50))
# Interpreting manhattan plot- these individuals generally share a large amount of SNPs
# Write out file

write.csv(vcf65.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion65.csv",
          quote=F,
          row.names=F)

# Looking at Hs values
names(vcf.div.MHplot)
# Hs values are stored in columns 4:9

# Density plot

vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x=value, fill=name))+
  geom_histogram(position="identity",alpha=0.5, bins=50)+
  labs(title="Genome-wide expected heterozygosity (Hs)", 
       fill="Regions",
       x="Gene diversity within regions",
       y="Counts of SNPs")
# value means generic value, coloring by column name
# alpha is transparency

# Saving the graph
ggsave("Histogram_GenomeDiversity_byRegion.pdf",
       path="~/Projects/eco_genomics/population_genomics/figures/")

# Making a summary table showing averages of the data we just plotted
vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n())
