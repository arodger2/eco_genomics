######
#
# HW #2 Ecological Genomics
#
######

# Load libraries
library(DESeq2)
library(eulerr)
library(ggplot2)
options(bitmapType = "cairo")

setwd("~/Projects/eco_genomics/transcriptomics/")

# Import counts matrix and round values
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header=TRUE, row.names = 1)
countsTableRound <- round(countsTable)

# Import experimental details
conds <-read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                   header=TRUE, stringsAsFactors = TRUE, row.names=1)
# DESeq 2
dds <-DESeqDataSetFromMatrix(countData=countsTableRound, colData=conds,
                             design= ~DevTemp + FinalTemp)

# Filtering based on average number of transcripts 
# If there aren't more than 10 reads for at least 15 out of 21 samples, we are ignoring
dds<- dds[rowSums(counts(dds) >= 10)>= 15,]

dds <- DESeq(dds)

# Part 1: 18_BASE versus 18_A28 and 18_BASE versus 18_A33
res_A_BASE_vs_A28 <- results(dds, name="FinalTemp_BASE_vs_A28", alpha=0.5)

# Start by making groups
dds$group <-factor(paste0(dds$DevTemp, dds$FinalTemp))

# sort dds by group
design(dds) <- ~group
dds <- DESeq(dds)
resultsNames(dds)

res_D18_BASE_D18_A28 <- results(dds, contrast = c("group","D18BASE","D18A28"), alpha=0.05)

res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[!is.na(res_D18_BASE_D18_A28$padj),]

res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[order(res_D18_BASE_D18_A28$padj),]

summary(res_D18_BASE_D18_A28)
#11 genes upregulated, 30 downregulated

# 18_BASE versus 18_A33
res_D18_BASE_D18_A33 <- results(dds, contrast = c("group","D18BASE","D18A33"), alpha=0.05)

res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[!is.na(res_D18_BASE_D18_A33$padj),]

res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[order(res_D18_BASE_D18_A33$padj),]

summary(res_D18_BASE_D18_A33)# 92 upregulated, 240 downregulated

#
degs_D18_BASE_D18_A28 <- row.names(res_D18_BASE_D18_A28[res_D18_BASE_D18_A28$padj < 0.05,])
degs_D18_BASE_D18_A33 <- row.names(res_D18_BASE_D18_A33[res_D18_BASE_D18_A33$padj < 0.05,])

length(degs_D18_BASE_D18_A28) #41
length (degs_D18_BASE_D18_A33) #332

length(intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33))#34

#unique DEGs at A28
41-34 #7

#unique DEGs at A33
332-34 #298

# Euler plot
EulerD18 <-euler(c("A28"=7, "A33"=298, "A28&A33"=34))
plot(EulerD18, lty=1:3, quantities=TRUE, fill=c("steelblue","gold","seagreen"), 
     alpha=0.4)

# combine results matrices
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D18BASE","D18A28"), alpha=0.05))
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group","D18BASE","D18A33"), alpha=0.05))
res_df18 <- merge(res_D18_BASEvsA28, res_D18_BASEvsA33, by="row.names", suffixes=c(".28",".33"))
rownames(res_df18) <- res_df18$Row.names


# color code genes significant in both contrasts
library(dplyr)
library(tidyverse)
res_df18 <- res_df18 %>% 
  mutate(fill=case_when(
    padj.28<0.05 & stat.28<0 ~ "turquoise2",
    padj.28<0.05 & stat.28>0 ~ "magenta2",
    padj.33<0.05 & stat.33<0 ~ "blue2",
    padj.33<0.05 & stat.33>0 ~ "red"
  ))



#Count the number of points per 
color_counts18 <- res_df18 %>% 
  group_by(fill) %>% 
  summarise(count=n())

label_positions18 <-data.frame(
 fill = c("blue2","magenta2","red","turquoise2"),
 x_pos=c(1,5,0,-7.5),
 y_pos=c(-5,0,9,3))

label_data18 <- merge(color_counts18, label_positions18, by="fill")

x_pos=c(1,5,0,-7.5)
y_pos=c(-5,0,9,3)

label_data18 <- merge(color_counts18, label_positions18, by="fill")

#Plot
plot18 <- ggplot(res_df18, aes(x=log2FoldChange.28, y=log2FoldChange.33, color=fill))+
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data18, aes(x=x_pos, y=y_pos, label=count, color=fill), size=5)+
  geom_abline(intercept=0, slope=1, linetype= "dashed", color="black")+
  geom_abline(intercept=0, slope=-1, linetype= "dashed", color= "grey")+
  xlim(-10,10)+ ylim(-10,10)+
  labs(x="log2Fold 28 vs Base",
       y="log2Fold 33 vs Base",
       title="Individuals developed at 18 response to 28 vs 33") +
  theme_minimal()
plot18

# Part 2: 22_BASE versus 22_A28 and 22_BASE versus 22_A33
res_D22_BASE_D22_A28 <- results(dds, contrast = c("group","D22BASE","D22A28"), alpha=0.05)

res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[!is.na(res_D22_BASE_D22_A28$padj),]

res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[order(res_D22_BASE_D22_A28$padj),]

summary(res_D22_BASE_D22_A28)#274 upregulated, 15 downregulated

res_D22_BASE_D22_A33 <- results(dds, contrast = c("group","D22BASE","D22A33"), alpha=0.05)

res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[!is.na(res_D22_BASE_D22_A33$padj),]

res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[order(res_D22_BASE_D22_A33$padj),]

summary(res_D22_BASE_D22_A33) #1176 upregulated, 388 downregulated


degs_D22_BASE_D22_A28 <- row.names(res_D22_BASE_D22_A28[res_D22_BASE_D22_A28$padj < 0.05,])
degs_D22_BASE_D22_A33 <- row.names(res_D22_BASE_D22_A33[res_D22_BASE_D22_A33$padj < 0.05,])

length(degs_D22_BASE_D22_A28)#289
length(degs_D22_BASE_D22_A33)#1564
length(intersect(degs_D22_BASE_D22_A28, degs_D22_BASE_D22_A33)) #144

#unique to 28
289-144 #145
#unique to 33
1564-144 #1420

# Euler plot
EulerD22 <-euler(c("A28"=145, "A33"=1420, "A28&A33"=144))
plot(EulerD22, lty=1:3, quantities=TRUE, fill=c("steelblue","gold","seagreen"), 
     alpha=0.4)

# combine results matrices
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D22BASE","D22A28"), alpha=0.05))
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group","D22BASE","D22A33"), alpha=0.05))
res_df22 <- merge(res_D22_BASEvsA28, res_D22_BASEvsA33, by="row.names", suffixes=c(".28",".33"))
rownames(res_df22) <- res_df22$Row.names

# color code genes significant in both contrasts
library(dplyr)
library(tidyverse)
res_df22 <- res_df22 %>% 
  mutate(fill=case_when(
    padj.28<0.05 & stat.28<0 ~ "turquoise2",
    padj.28<0.05 & stat.28>0 ~ "magenta2",
    padj.33<0.05 & stat.33<0 ~ "blue2",
    padj.33<0.05 & stat.33>0 ~ "red"
  ))

#Count the number of points per 
color_counts22 <- res_df22 %>% 
  group_by(fill) %>% 
  summarise(count=n())

label_positions22 <-data.frame(
  fill = c("blue2","magenta2","red","turquoise2"),
  x_pos=c(1,5,0,-7.5),
  y_pos=c(-5,0,9,3))

label_data22 <- merge(color_counts22, label_positions22, by="fill")

#Plot
plot22 <- ggplot(res_df22, aes(x=log2FoldChange.28, y=log2FoldChange.33, color=fill))+
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data22, aes(x=x_pos, y=y_pos, label=count, color=fill), size=5)+
  geom_abline(intercept=0, slope=1, linetype= "dashed", color="black")+
  geom_abline(intercept=0, slope=-1, linetype= "dashed", color= "grey")+
  xlim(-10,10)+ ylim(-10,10)+
  labs(x="log2Fold 28 vs Base",
       y="log2Fold 33 vs Base",
       title="Individuals developed at 22 response to 28 vs 33") +
  theme_minimal()
plot22


