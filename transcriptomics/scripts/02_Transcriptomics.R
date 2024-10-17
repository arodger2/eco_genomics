
# Load in day 1 transcriptomic script and run DESeq object

library(pheatmap)
resultsNames(dds) #dds file contains all comparisons of gene expression
options(bitmapType = "cairo")

# [1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28" 
# [4] "FinalTemp_BASE_vs_A28"

# Saving the comparisons of interest into a "results file"
# First, comparing developmental temperature 22 vs 18

res_D22vsD18 <- results(dds, name="DevTemp_D22_vs_D18", alpha=0.5)

# Order by significance of difference in expression between developmental temperatures
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]

# Look at counts of a specific top gene we're interested in to validate that the model is working

d <- plotCounts(dds, gene="TRINITY_DN140854_c0_g5_i2", int=(c("DevTemp", "FinalTemp")),
                returnData = T)
d

p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp))+
  theme_minimal()+theme(text=element_text(size=20), panel.grid.major=(element_line(colour="gray")))+ 
  geom_point(position = position_jitter(w=0.2, h=0), size=3)

#  Making an MA plot (m is another term for logfold change, a means average)
plotMA(res_D22vsD18, ylim=c(-4,4))

# Volcano plot
# convert our deseq results object into a data frame

res_df <- as.data.frame(res_D22vsD18)
# add a column to this data frame to label whether a gene is significantly differentially expressed
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange)> 1, "Significant", "Not Significant")


# Plot 
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color= Significant)) +
  geom_point(alpha= 0.8)+
  scale_color_manual(values=c("slateblue", "tomato"))+
  labs (x="Log2 Fold Change", y="log 10 adjusted p value", title="Volcano")+
  theme_linedraw()+
  theme(legend.position = "top")+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="orange")+
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="orange")

# Heatmap to look at expression levels
vsd <- vst(dds, blind=FALSE)

# Focusing on top genes in results data frame (most significant)
topgenes <- head(rownames(res_D22vsD18), 20)

# Make a matrix of gene expression data across all samples in data ste
mat <- assay(vsd)[topgenes, ]
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T, cluster_rows=T)




