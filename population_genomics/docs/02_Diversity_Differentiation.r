# Estimating diversity and differentiation

library(vcfR)
library(tidyverse)
library(qqman)

X11.options(type="cairo")

# Read in vcf file 
vcf <-read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf.filtered.vcf.gz")

# Read in metadata
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
head(meta)

# vcf has 593 samples
dim(meta)

# making them compatable

meta2 <-meta[meta$id %in% colnames(vcf@gt[,-1]),]
dim(meta2)

# calculate diversity stats using the genetic_diff fxn in vcfR
vcf.div <- genetic_diff(vcf,
                        pops=as.factor(meta2$region),
                        method="nei")

str(vcf.div)

# Taking just first 8 entries (8 main chromsomes)
chr.main <- unique(vcf.div$CHROM) [1:8]

# Tell it which chromosome numbers

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))
chrnum

# merge diversity stats with chrnum- assign to chromosome numbers

vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))
vcf.div.MHplot <- vcf.div.MHplot %>%
  filter(Gst>0) %>% 
  mutate(SNP=paste0(chr.main, "_", POS))

# Make sure numbers are viewed as numbers
vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)
  
# Manhattan plot with qqman
manhattan(vcf.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.50))
# Interpreting manhattan plot- these individuals generally share a large amount of SNPs
# Write out file

write.csv(vcf.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.csv",
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
  filter(value!=0) %>% 
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

# In introduced areas, there's more opportunities for mixing?



 

