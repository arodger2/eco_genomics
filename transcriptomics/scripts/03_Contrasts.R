
# Loading new library
library(eulerr)

# Start by making groups
dds$group <-factor(paste0(dds$DevTemp, dds$FinalTemp))

# sort dds by group
design(dds) <- ~group
dds <- DESeq(dds)
resultsNames(dds)

# [1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28"
# [4] "group_D22A28_vs_D18A28"  "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

# 1. Setting up first contrsast, creating an object with our results for first baseline

res_D18_BASE_D22_BASE <- results(dds, contrast = c("group","D18BASE","D22BASE"), alpha=0.05)

res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]

res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]

summary(res_D18_BASE_D22_BASE) 

# 556 genes upregulated
# 1379 genes downregulated

# in our deseq object, each row is a gene
# make a list of which genes in our comparisons of interest are differentially expressed (DEGs)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])


# MA of this contrast 
plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

# 2. Comparing at A28
res_D18_A28_D22_A28 <- results(dds, contrast = c("group","D18A28","D22A28"), alpha=0.05)

res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]

res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]

summary(res_D18_A28_D22_A28) 

# 146 genes upregulated
# 150 genes downregulated

# in our deseq object, each row is a gene
# make a list of which genes in our comparisons of interest are differentially expressed (DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < 0.05,])


# MA of this contrast 
plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))


# 3. Comparing at A33
res_D18_A33_D22_A33 <- results(dds, contrast = c("group","D18A33","D22A33"), alpha=0.05)

res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]

res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),]

summary(res_D18_A33_D22_A33) 

# 47 genes upregulated
# 31 genes downregulated

# Make a list
degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj < 0.05,])


# MA of this contrast 
plotMA(res_D18_A33_D22_A33, ylim=c(-4,4))

# Now, build eulerr plot using lists we just made (all DEGs, not separated by up vs down regulation)

length(degs_D18_BASE_D22_BASE)# 1935
length(degs_D18_A28_D22_A28) # 296
length(degs_D18_A33_D22_A33) # 78

# Look at overlap between which genes are differentially expressed in multiple contrasts
length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)) #107
length(intersect(degs_D18_A28_D22_A28, degs_D18_A33_D22_A33)) #29
length(intersect(degs_D18_A33_D22_A33, degs_D18_BASE_D22_BASE)) #44

length(intersect(degs_D18_BASE_D22_BASE, intersect(degs_D18_A28_D22_A28, degs_D18_A33_D22_A33)))
#23

# Calculate the number of unique genes in each portion of the euler plot

1935-107-44+23 #1807 (unique) differentially expressed between 18 and 22 at base conditions

296-107-29+23 #183 genes diff expressed unique at A28

78-44-29+23 #28 genes differentially expressed uniquely at A33

# Degree of overlap between
107-23 #84 DEGs overlapping in base and A28
44-23 #21 genes unique to base and A33
29-23 #6 unique to A28 and A33

# Ready to plot!
myEuler <-euler(c("BASE"=1807, "A28"=183, "A33"=28, "BASE&A28"=84, 
                  "BASE&A33"=21, "A28&A33"=6, "BASE&A28&A33"=23))

plot(myEuler, lty=1:3, quantities=TRUE)

##################
#
# Make a scatter plot of responses to A28 (Then A33) when copepods develop at 18 vs 22
#
##################


# contrast D18_A28 vs BASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D18BASE","D18A28"), alpha=0.05))

# Contrast D22_A28 vs BASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group","D22BASE","D22A28"), alpha=0.05))

# Merge these data frames
res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by="row.names", suffixes=c(".18",".22"))

rownames(res_df28) <- res_df28$Row.names

# keep all rows get rid of first column
res_df28 <- res_df28[,-1]

library(dplyr)
library(tidyr)

# Define our color mapping logic with the mutate function (showing up and down regulated genes in each conditrin)
res_df28 <- res_df28 %>% 
  mutate(fill=case_when(
    padj.18<0.05 & stat.18<0 ~ "turquoise2",
    padj.18<0.05 & stat.18>0 ~ "magenta1",
    padj.22<0.05 & stat.22<0 ~ "blue2",
    padj.22<0.05 & stat.22>0 ~ "red"
  ))

#Count the number of points per 
color_counts1 <- res_df28 %>% 
  group_by(fill) %>% 
  summarise(count=n())

label_positions1 <-data.frame(
  fill = c("blue2","magenta2","red","turquoise2"),
  x_pos=c(1,5,0,-7.5),
  y_pos=c(-5,0,9,3))

label_data1 <- merge(color_counts1, label_positions1, by="fill")


#Plot
plot28 <- ggplot(res_df28, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill))+
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data1, aes(x=x_pos, y=y_pos, label=count, color=fill), size=5)+
  geom_abline(intercept=0, slope=1, linetype= "dashed", color="black")+
  geom_abline(intercept=0, slope=-1, linetype= "dashed", color= "grey")+
  xlim(-10,10)+ ylim(-10,10)+
  labs(x="log2FoldChange 28 vs. BASE at 18",
       y="log2FoldChange 28 vs. BASE at 22",
       title="How does response to 28C vary by dev temp?") +
  theme_minimal()
plot28
##################

#
# Make a scatter plot of responses to A33 when copepods develop at 18 vs 22
#
##################


# contrast D18_A33 vs BASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group","D18BASE","D18A33"), alpha=0.05))

# Contrast D22_A33 vs BASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group","D22BASE","D22A33"), alpha=0.05))

# Merge these data frames
res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by="row.names", suffixes=c(".18",".22"))

rownames(res_df33) <- res_df33$Row.names

# keep all rows get rid of first column
res_df33 <- res_df33[,-1]

library(dplyr)
library(tidyr)

# Define our color mapping logic with the mutate function (showing up and down regulated genes in each conditrin)
res_df33 <- res_df33 %>% 
  mutate(fill=case_when(
    padj.18<0.05 & stat.18<0 ~ "turquoise2",
    padj.18<0.05 & stat.18>0 ~ "magenta1",
    padj.22<0.05 & stat.22<0 ~ "blue2",
    padj.22<0.05 & stat.22>0 ~ "red"
  ))
#Count the number of points per 
color_counts2 <- res_df33 %>% 
  group_by(fill) %>% 
  summarise(count=n())

label_positions2 <-data.frame(
  fill = c("blue2","magenta2","red","turquoise2"),
  x_pos=c(1,5,0,-7.5),
  y_pos=c(-5,0,9,3))

label_data2 <- merge(color_counts2, label_positions2, by="fill")




#Plot
plot33 <- ggplot(res_df33, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill))+
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data2, aes(x=x_pos, y=y_pos, label=count, color=fill, size=5))+
  geom_abline(intercept=0, slope=1, linetype= "dashed", color="black")+
  geom_abline(intercept=0, slope=-1, linetype= "dashed", color= "grey")+
  xlim(-10,10)+ ylim(-10,10)+
  labs(x="log2FoldChange 33 vs. BASE at 18",
       y="log2FoldChange 33 vs. BASE at 22",
       title="How does response to 33C vary by dev temp?") +
  theme_minimal()

library(gridExtra)

# Put the two plots together in a two panel plot
combined_plot <- grid.arrange(plot28, plot33, ncol=2)

ggsave("~/Projects/eco_genomics/transcriptomics/figures/combined_scatter_plot.png", combined_plot, width=12, height=6)