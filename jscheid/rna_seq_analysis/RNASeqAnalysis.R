if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library("RColorBrewer")
library("DESeq2")
library("sva")
library("rafalib")
library("ggplot2")
library("stringr")


setwd("~/Documents/Uni/Master/Thesis/RNA")

rlog = read.csv('RNAseq_results_June/differential_gene_expression/gene_counts_tables/rlog_transformed_gene_counts.tsv', sep = '\t')
vst = read.csv('RNAseq_results_June/differential_gene_expression/gene_counts_tables/vst_transformed_gene_counts.tsv', sep = '\t')
#de_genes = read.csv('RNAseq_results_June/differential_gene_expression/DE_genes_tables/DE_contrast_condition_tissues_recurrence_vs_primary.tsv', sep = '\t')
final = read.csv('RNAseq_results_June/differential_gene_expression/final_gene_table/final_DE_gene_list.tsv', sep='\t')
colnames(vst)[3:46] = ifelse(grepl("Prim", colnames(vst)[3:46] ), paste0("GBM", gsub(".*_GBM(.+)_Prim.", "\\1", colnames(vst)[3:46] ), "P"),
                                                                  paste0("GBM", gsub(".*_GBM(.+)_Rec.", "\\1", colnames(vst)[3:46] ), "R"))

colnames(rlog)[3:46] = ifelse(grepl("Prim", colnames(rlog)[3:46] ), paste0("GBM", gsub(".*_GBM(.+)_Prim.", "\\1", colnames(rlog)[3:46] ), "P"),
                                                                    paste0("GBM", gsub(".*_GBM(.+)_Rec.", "\\1", colnames(rlog)[3:46] ), "R"))
#DE genes with padj < 0.05 and |log2FC| > 1
de_genes = which(final['padj_condition_tissues_recurrence_vs_primary'] < 0.05 & abs(final['log2FoldChange_condition_tissues_recurrence_vs_primary']) > 1)

# Scatterplot of GBM2 primary vs recurrent
plot(unlist(rlog['GBM2P']), unlist(rlog['GBM2R']), pch=19, cex = 0.2)
points(unlist(rlog[de_genes,]['GBM2P']), unlist(rlog[de_genes,]['GBM2R']), pch=19, cex = 0.5, col='red')

#heatmap
heatmap(as.matrix(rlog[de_genes,-c(1:2)]), labRow = "", cexCol = 0.6)
heatmap(as.matrix(vst[de_genes,-c(1:2)]), labRow = "", cexCol = 0.6) # expot with size 12x12

# plot displaying de gene expression increase or decrease between prim and rec
colnames(rlog)
for (i in 2:23) {
  print(colnames(rlog)[grep(colnames(rlog)[], colnames(rlog))])
}

###########################################################################
#           From PCA exploration derived batch corrected data             #
###########################################################################
# Read in raw count file
rawData <- read.table('RNAseq_results_June/differential_gene_expression/gene_counts_tables/raw_gene_counts.tsv',
                      sep = '\t', header=T)
geneSub <- rawData[,1:2]
rawData <- rawData[,-c(1:2)]
rownames(rawData) <- geneSub$Ensembl_ID
for (i in 1:ncol(rawData)){
  colnames(rawData)[i] <- gsub('[.]','',paste(tail(unlist(strsplit(colnames(rawData)[i],
                                                                   split = "_")),2), collapse = "_"))
} 
# Read metadata file
target <- read.csv2('PCA/PCA_metadata.csv')
target$PatientID <- gsub("_.*","", target$SampleID)
# Sort target and raw data such that they can be used in downstream analysis
target <- target[order(target$SampleID),]
rawData <- rawData[,order(colnames(rawData))]
# Create a DESeq2 model based on raw counts
dds <- DESeqDataSetFromMatrix(countData = rawData, colData = target, design = ~ PatientID + phenotype)
## Filter lowely expressed genes (genes with a rowsum > 1 are kept)
dds <- dds[ rowSums(counts(dds)) > 1, ]
## The data is normalized based on the vst transformation
vsd <- vst(dds, blind = FALSE)

# The counts are saved to a new vector
counts_deseq <- assay(vsd)
# A H1 is formulated, meaning that there might be variation based on age and gender
mod <- model.matrix(~ phenotype, target)
# The H0 is created
mod0 <- model.matrix(~ 1, target)
# svaseq makes sure that the noise that is depicted in the H0 is removed (depleted) in H1
svseq <- svaseq(counts_deseq, mod, mod0)
# The folowing function enables the generation of a new dataframe based on the noise that was depicted with svaseq
cleaningY = function(y, mod, svaobj) {
  X = cbind(mod, svaobj$sv)
  Hat = solve(t(X)%*%X)%*%t(X)
  beta = (Hat%*%t(y))
  P = ncol(mod)
  cleany = y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
# A new dataframe (without batch effect) is generated
cleandat = cleaningY(counts_deseq, mod, svseq)

###########################################################################
#           DE analysis of BC data                                        #
###########################################################################
# Check which test should be applied to the BC data:
# - Check if the data is parametric: The mean is a proper 
# representation of the data distribution
hist(cleandat, breaks = 20)
abline(v=mean(cleandat))
# Compare it to the not BC data
hist(as.matrix(vst[,3:ncol(vst)]), breaks = 20)
abline(v=mean(as.matrix(vst[,3:ncol(vst)])))
# Conclusion: Comparing both plots the distribution and mean are similar located
#             therefore t-test is a proper test for the underlying distribution

# Following segment inspired by: http://genomicsclass.github.io/book/pages/adjusting_with_linear_models.html
res <- genefilter::rowttests(cleandat,factor(as.numeric(grepl("Rec", colnames(cleandat)))))
res$padj <- res$p.value * nrow(cleandat)
res$log2fc <- apply(cleandat, MARGIN = 1, function(row){
  prim_col <- which(grepl("Prim", colnames(cleandat)))
  rec_col <- which(grepl("Rec", colnames(cleandat)))
  return(log2(mean(row[rec_col]) / mean(row[prim_col])))
})

# Vulcano plot
ggplot(data=res, aes(x = res$dm, y = -log10(res$padj))) +
  geom_point(alpha=0.5, col = ifelse(res$padj > 0.05, "black", brewer.pal('RdBu', n = 4)[1])) +
  geom_hline(yintercept = -log10(0.05), alpha=0.3) +
  theme_classic() +
  xlab("Difference of group means") +
  ylab(expression("-log"[10]*"(p"[adj]*")"))
  
bc_de_genes <- geneSub[match(rownames(res)[which(res$padj < 0.05)], geneSub$Ensembl_ID), ]

match(final$gene_name[de_genes], bc_de_genes)


