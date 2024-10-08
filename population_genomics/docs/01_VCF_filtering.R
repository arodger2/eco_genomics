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
vcf.filt.indMiss <-missing_by_sample(vcf.filt,
                                     popmap = meta2,
                                     cutoff=0.75)


vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)

# sets the minimum number of times to see SNP
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

# Snp-wise missing alleles
vcf.filt.ind.SNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)

DP2 <-extract.gt(vcf.filt.ind.SNPMiss,
                 element="DP",
                 as.numeric = T)
heatmap.bp(DP2, rlabels = F, clabels = F)

write.vcf(vcf.filt.ind.SNPMiss, 
          "~/Projects/eco_genomics/population_genomics/outputs/vcf.filtered.vcf.gz")


