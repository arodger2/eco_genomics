
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

