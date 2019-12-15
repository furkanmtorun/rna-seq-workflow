#
# +-----------------------------------------------------------------------------------------+
# +-----------------------------------------------------------------------------------------+
# RNA-Seq Workflow by @furkanmtorun 
    # E-mail: furkanmtorun@gmail.com  
    # GitHub: https://github.com/furkanmtorun 
    # Google Scholar: https://scholar.google.com/citations?user=d5ZyOZ4AAAAJ 
    # Personal Website: https://furkanmtorun.github.io/    
# +-----------------------------------------------------------------------------------------+
# +-----------------------------------------------------------------------------------------+
#


# +--------------------------------------------------+
# Install and import required libraries & packages
# +--------------------------------------------------+
# install.packages("dplyr")
# install.packages("calibrate")
# BiocManager::install("DESeq2")
library(dplyr)
library(ggplot2)
library(DESeq2)
library(calibrate)
library(pheatmap)
library(gplots)
library(RColorBrewer)


# +--------------------------------------------------+
# Argument parsing
# +--------------------------------------------------+

options(echo=FALSE)
args<-commandArgs(TRUE)

# Arguments order must be: counts_csv_file_path meta_data_csv_file_path sample_name_for_control sample_name_for_sample
if (length(args) != 4) {
    stop("At least four arguments must be supplied!")
}


# +--------------------------------------------------+
# Set and define the files & DESeq2
# +--------------------------------------------------+
counts_csv <- args[1]
meta_data_csv <- args[2]

# Import countdata
countData = read.csv(counts_csv, row.names=1) %>%
    dplyr::select(-length) %>% 
    as.matrix()

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]

# Import meta_data
colData = read.csv(meta_data_csv, row.names=1)

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds = DESeq(dds)

# Get results for the sample versus control, sort the file by p-value. 
res = results(dds, contrast=c("condition", args[4], args[3]))
res = res[order(res$pvalue),]
write.csv(as.data.frame(res), file="files/results/DESeq2.csv")

#PCA Plot
png(filename = "files/results/PCA_plot.png")
plotPCA(rlog(dds))
dev.off()

# Sparsity Plot
png(filename = "files/results/Sparsity_plot.png")
plotSparsity(estimateSizeFactors(dds))
dev.off()

# MA Plot
png(filename = "files/results/MA_plot.png")
plotMA(results(dds))
dev.off()

# Gene based Plot Counts
png(filename = "files/results/Particular_gene_plot.png")
plotCounts(dds, gene="ENSG00000103196", intgroup="condition")
dev.off()
