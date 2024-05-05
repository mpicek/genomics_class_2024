library(DESeq2)
library(pheatmap)
alpha = 0.05

# Load count matrix
x = read.table("ex_3/sample.counts", row.names=1, header=T, sep=",") # count matrix
s = read.table("ex_3/sample.info", header=T, row.names=1, colClasses=c("character", "factor"))

# Create DESeq2 object
dds = DESeqDataSetFromMatrix(countData = x, colData = s, design = ~ condition)

# Run a differential expression analysis (Tumour vs. Normal) using a log-fold change threshold of 1
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
dds <- DESeq(dds)
results <- results(dds, contrast=c("condition","Tumour","Normal"), alpha=alpha, lfcThreshold=1)

head(results[order(results$padj),], n=10)

# Generate an MA-plot
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results

plotMA(results, alpha=alpha)
resultsLFC <- lfcShrink(dds, coef="condition_Tumour_vs_Normal", type="apeglm")
plotMA(resultsLFC, alpha=alpha)

# SAVING THE PLOTS FOR THE REPORT:
# png("ex_3/MA_plot_noise.png", width = 800, height = 600)
# plotMA(results, alpha=alpha)
# dev.off()
# png("ex_3/MA_plot.png", width = 800, height = 600)
# plotMA(resultsLFC, alpha=alpha)
# dev.off()

# Plot the normalized counts for the GJB2 gene
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts
plotCounts(dds, gene="GJB2", intgroup="condition")

# SAVING THE PLOTS FOR THE REPORT:
# png("ex_3/GJB2_normalized.png", width = 600, height = 600)
# plotCounts(dds, gene="GJB2", intgroup="condition")
# dev.off()

# Generate a PCA plot of the samples using the transformed count data
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extracting-transformed-values
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples


# Visualize the differential gene expression results as a heatmap
# Take the top 20 genes according to the adjusted p-value
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix


# Export the significant results (padj < 0.01) to a CSV file
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exporting-results-to-csv-files

