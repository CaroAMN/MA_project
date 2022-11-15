################################################################################
###                                Purpose                                   ###
################################################################################
## DESCRIPTIOM: This script performs differential expression analysis on Salmon
##              quantification files, based on the DESeq2 package.
##
##
## INPUT: Salmon quantification files
## OUTPUT: for each pairwise comparison: - count table with DE expression
##                                         Annotation --> SIG, DOWN, UP
##                                       - Volcano plot
##                                       - MA plot
################################################################################
###                             References                                   ###
################################################################################
## Differential expression analysis and testing is based on the DEseq2 package
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html



################################################################################
###                             Manual Steps                                 ###
################################################################################
# clear R environment
rm(list = ls())
# set directory of the input files for the analysis
input_dir <- "/Users/cschwitalla/Documents/transcriptomics_results/rnaseq_results_GRCh38/"
input_dir_xoutlier <- "/Users/cschwitalla/Documents/transcriptomics_results/rnaseq_results_GRCh38_without outlier/"
# set output directory for results table and plots
output_dir <- "/Users/cschwitalla/Documents/RNAseq_analysis/"
# set directory of the metadata file
metadata_file <- "/Users/cschwitalla/Documents/RNAseq_analysis/Metadata_GB.tsv"
# set directory of the script that contains RNAseq analysis functions
source("/Users/cschwitalla/git/students/cschwitalla/RNAseq_Analysis/functions_RNA.R")
# set lfc and padju values for significance thresholds
lfc <- 3
padju <- 0.01
################################################################################
###                           Load libraries                                 ###
################################################################################
required_libs <- c(
  "DESeq2", "tximport", "ashr", "ggVennDiagram", "ggplot2",
  "dplyr", "readxl", "tidyr", "pheatmap", "RColorBrewer",
  "EnsDb.Hsapiens.v86", "ggrepel", "utils", "EnhancedVolcano", "limma"
)


suppressMessages(invisible(lapply(required_libs, library, character.only = TRUE)))


################################################################################
###                            Load Data                                     ###
################################################################################
# Load meta data --> Metadata_GB.tsv in workdir
metadata <- read.table(file = metadata_file, sep = "\t", header = TRUE)
metadata2 <- metadata[-grep(("QATLV139AX"),metadata$QBiC.Code),]
# get filenames of inputDir
file_names <- list.files(path = input_dir)
# get file names without outlier
file_names_xoutl <- list.files(path = input_dir_xoutlier)


# files without ben + outlier sample
#filnames_excl <- grep(("NEC|INF|T1"), file_names, value = TRUE)
#filnames_excl <- filnames_excl[c(1:7,9:45)]
# TO DO: create/ import tx2gene dataframe
#tx2gene <- create_tx2gene()
tx2gene <- read.csv(paste0(input_dir, "../salmon_tx2gene.tsv"),sep = "\t", header = FALSE)
# import Salmon quant files with tximport and tx2gene
txi <- tximport(paste0(input_dir_xoutlier, file_names_xoutl), type = "salmon", tx2gene = tx2gene)
#txi2 <- tximport(paste0(input_dir, filnames_excl), type = "salmon", tx2gene = tx2gene)

# make Deseq2 object of raw read counts and specify metadata and project design
dds <- DESeqDataSetFromTximport(txi,
  colData = metadata2,
  design = ~ Tumor_region
)
# dds2 <- DESeqDataSetFromTximport(txi2,
#                                 colData = metadata2,
#                                 design = ~Tumor_region
# )
# pre filtering, keep only rows with at least 10 reads total
# WHY:reduces computation time and memory usage
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
#dds2 <- dds2[keep, ]
# model read counts, estimate LFCs and test for DE expressed genes
# with default parameters
dds_default <- DESeq(dds)
#dds_minor <- DESeq(dds2)

dds_default <- estimateSizeFactors(dds_default)
counts_data <- counts(dds_default, normalized = TRUE)
################################################################################
###                            DE_Analysis                                   ###
################################################################################
# list of genes that are diff expressed
#TODO: safe all DEgenes per comparison 
comparisons <- c()
results_list <- list()
DE_list <- list()
DE_genes <- c()
# get each pairwise comparison for condition that should be tested
for (com in apply(combn(levels(dds_default$Tumor_region), 2), 2, paste, collapse = " vs ")) {
  # get condition vector for DE analysis function
  comparison <- unlist(strsplit(com, " vs "))
  comparisons <- append(comparisons, com)
  # make results data frame
  res_df <- DE_analysis(dds_default, c("Tumor_region", comparison), lfc, padju)
  #store results df
  results_list <- append(results_list, list(res_df))
  # get DE genes and append them into DE_genes list
  DE_genes <- append(DE_genes, res_df$Names[res_df$significant == "SIG"])
  DE_list <- append(DE_list, list(res_df@listData[c("diffexp", "Names", "log2FoldChange", "padj")]))
  # write results as tsv into the workdir to stor results
  write.table(res_df, file = paste(output_dir, com, ".tsv"), sep = "\t", quote = FALSE)
  # make MA plot and save it as png in output dir
   pdf(paste0(output_dir, "MAplot_", com, ".pdf"))
   DESeq2::plotMA(res_df, xlim = c(1, 1e5), ylim = c(-4, 4), main = com)
   dev.off()
  # make Volcano plot and save it as pdf in output dir
  ggsave(make_volcano(res_df, "padj", com, lfc, padju, TRUE),file = paste0(output_dir, "Volcano_", com, ".pdf"),width = 15, height = 13) 
  #dev.off()
}
names(DE_list) <- comparisons
names(results_list) <- comparisons
################################################################################
###                              Volcano                                     ###
################################################################################
write.table(unique(DE_genes), "/Users/cschwitalla/Documents/RNAseq_analysis/DE_gene_list.tsv", quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)

#
# exome_genes <- c("TXNIP", "DNAH8", "OR2A1", "PAX3", "ANKRD30A")
# highlight_genes <- c("EGFR", "LNP1", "OR8U1", "CHIT1", "CSF-1R", "VEGF",
#                      "IL-8", "HIF-1", "SOX2", "OCT4", "CD133", "TXNIP", "DNAH8", "OR2A1", "PAX3", "ANKRD30A")
# 
# test <- data.frame(DE_genes)
genes_nec_t1 <- c("TRNAL19", "LOC101928747", "TRNAL-CAA", "NFE4", "DEFA1", "SNORA17",
                   "CXCR2", "RNU11", "VTRNA2−1", "BCL2A1", "MME", "IL18R1",
                  "MTRNR2L9", "C19orf59", "VTRNA2-1")
library(gridExtra)
library(grid)

v1 <- make_volcano(results_list$`BEN vs INF`,"padj", "BEN vs INF", lfc, padju, NULL)
#v2 <- make_volcano(results_list$`INF vs T1`,"padj", "INF vs T1", lfc, padju, NULL)
#v3 <-make_volcano(results_list$`NEC vs T1`,"padj", "NEC vs T1", lfc, padju, genes_nec_t1)
recurent_genes <- c("MYO1G", "LIF", "CA9", "STC2", "TGFBI", "ESM1", "GPRC5A", "CHRDL2", "VEGFA")
necben <- data.frame(DE_list$`BEN vs NEC`)
necben <- subset(necben, Names %in% recurent_genes)
t1ben <- data.frame(results_list$`BEN vs T1`)
t1ben <- subset(t1ben, Names %in% recurent_genes)
infben <- data.frame(DE_list$`BEN vs INF`)
infben <- subset(infben, Names %in% recurent_genes)
length(unique(DE_genes))

vb1 <- make_volcano(results_list$`BEN vs INF`,"padj", "BEN vs INF", lfc, padju, NULL)
vb2 <- make_volcano(results_list$`BEN vs T1`,"padj", "BEN vs T1", lfc, padju, NULL)
vb3 <- make_volcano(results_list$`BEN vs NEC`,"padj", "BEN vs NEC", lfc, padju, NULL)

T1_vs_NEC <- DE_analysis(dds_default, c("Tumor_region", c("T1", "NEC")), lfc, padju)


v1 <- make_volcano(results_list$`BEN vs INF`,"padj", "BEN vs INF", lfc, padju, NULL)
v2 <- make_volcano(results_list$`INF vs T1`,"padj", "INF vs T1", lfc, padju, NULL)
v3 <- make_volcano(T1_vs_NEC,"padj", "T1 vs NEC", lfc, padju, genes_nec_t1)

make_volcano(T1_vs_NEC,"padj", "T1 vs NEC", lfc, padju, NULL)

make_volcano(results_list$`NEC vs T1`,"padj", "NEC vs T1", lfc, padju, NULL)
  
grid.arrange(v1, v2, v3,
             ncol = 3)

grid.arrange(vb1, vb2, vb3,
             ncol = 3)

#make_volcano(results_list$`INF vs NEC`, "padj", com, lfc, padju, NULL)
nt <- data.frame(T1_vs_NEC)
#######################################solution for outliers####################
# https://www.biostars.org/p/9481664/#9481786
res.new <- results_list$`INF vs T1`
cutoff <- 10

up <- rownames(subset(results_list$`INF vs T1`, log2FoldChange >= cutoff))
down <- rownames(subset(results_list$`INF vs T1`, log2FoldChange <= cutoff*-1))

# impute log2 FCs in the object to be max 2.5 or min -2.5  
res.new$log2FoldChange <- ifelse(res.new$log2FoldChange >= cutoff, cutoff,
                                 ifelse(res.new$log2FoldChange <= cutoff * -1, cutoff * -1,
                                        res.new$log2FoldChange))
max(res.new$log2FoldChange, na.rm = TRUE)
min(res.new$log2FoldChange, na.rm = TRUE)

# custom shapes for the points
customshape <- rep(19, nrow(res.new))
names(customshape) <- rep('group1', nrow(res.new))
customshape[which(rownames(res.new) %in% up)] <- -9658
names(customshape)[which(rownames(res.new) %in% up)] <- 'group2'
customshape[which(rownames(res.new) %in% down)] <- -9668
names(customshape)[which(rownames(res.new) %in% down)] <- 'group3'

customsize <- rep(2.0, nrow(res.new))
customsize [which(rownames(res.new) %in% up)] <- 8
customsize [which(rownames(res.new) %in% down)] <- 8


make_volcano(res.new, "padj", com, lfc, padju, NULL, customshape, customsize)


test <- data.frame(results_list$`NEC vs T1`@listData)
###############################################################################
test <- data.frame(DE_list$`BEN vs T1`)
test2 <- test[test$Names == "TXNIP",]

test <- assay(vsd)
test2 <- test[rownames(test) %in% c("TXNIP", "DNAH8"),]
################################################################################
###                              Heatmap                                     ###
################################################################################
set.seed(8)
# varianze stabilizing transformation of the data to plot as heatmap
vsd <- vst(dds_default, blind = FALSE)

# create list of annotation colors for heatmap plotting
annotation_color <- create_annotation_color(
  metadata2$Patient_ID,
  metadata2$Tumor_region,
  metadata2$Sex,
  metadata2$MGMT,
  #gene_cluster$gene_cluster
)

# plot heatmap with bacth corection for Patientens
cluster_col_data <- make_heatmap(de_genes$V1, vsd, vsd$Patient_ID, NULL, NULL)
cluster_row_data <- make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, k = NULL)

#clusters <- cutree(tree = as.dendrogram(cluster_row_data$tree_row), k = 7)
#clusters <- data.frame(clusters )

#clusters 
#get gene cliszers
gene_cluster <- sort(cutree(cluster_row_data$tree_row, k = 6))
#plot(cluster_row_data$tree_row)
gene_cluster <- data.frame(gene_cluster)
gene_cluster$gene_cluster <- as.character(gene_cluster$gene_cluster)
gene_cluster <- split(gene_cluster, gene_cluster$gene_cluster)

genes_c1 <- rownames(gene_cluster$`1`)
write.csv(genes_c1, paste0(output_dir,"_c1_genes.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE )
genes_c2 <- rownames(gene_cluster$`2`)
write.csv(genes_c2, paste0(output_dir,"_c2_genes.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_c3 <- rownames(gene_cluster$`3`)
write.csv(genes_c3, paste0(output_dir,"_c3_genes.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_c4 <- rownames(gene_cluster$`4`)
write.csv(genes_c4, paste0(output_dir,"_c4_genes.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_c5 <- rownames(gene_cluster$`5`)
write.csv(genes_c5, paste0(output_dir,"_c5_genes.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_c6 <- rownames(gene_cluster$`6`)
write.csv(genes_c6, paste0(output_dir,"_c6_genes.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_c7 <- rownames(gene_cluster$`7`)
write.csv(genes_c7, paste0(output_dir,"_c7_genes.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
#mak row annotation 

#make heatmap with row annotation 

make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, gene_cluster)
make_heatmap(DE_genes, vsd, NULL, annotation_color, k = NULL)

#get tree in more detail 
plot(cluster_col_data$tree_col)
plot(cluster_col_data$tree_row)
dend <- as.dendrogram(cluster_col_data$tree_col)
dend_row <- as.dendrogram(cluster_col_data$tree_row)
labels(dend) <- vsd$Patient_ID
plot(dend_row)
# plot without batch correction
make_heatmap(DE_genes, vsd, batch = NULL, annotation_color, k = NULL)

################
annotation_df <- data.frame(
  "Tumor_region" = vsd@colData@listData$Tumor_region,
  "Patient_ID" = vsd@colData@listData$Patient_ID,
  "Sex" = vsd@colData@listData$Sex,
  "MGMT_methylation" = vsd@colData@listData$MGMT
)

#################get cluster data 

cluster_row_data <- make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, gene_cluster)
plot(cluster_row_data$tree_row)
dev.off()
clusters_rows <- sort(cutree(cluster_row_data$tree_row, k = 6))
plot(cluster_row_data$tree_row)
clusters_rows <- data.frame(clusters_rows)
clusters_rows$genes <- row.names(clusters_rows)
clusters_rows <- split(clusters_rows, clusters_rows$clusters_rows)

cvsd$Patient_ID
cluster_col_data <- make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, k = NULL)
cluster_col_data$tree_col$labels <- paste0(vsd$Tumor_region, "_",vsd$Patient_ID)
clusters_cols <- sort(cutree(cluster_col_data$tree_row, k = 5))
plot(cluster_col_data$tree_col)
library("dendextend")
color_col_dend <- color_branches(cluster_col_data$tree_col, k = 5)
clusters_cols <- data.frame(clusters_cols)
clusters_cols$genes <- row.names(clusters_cols)
clusters_cols <- split(clusters_cols, clusters_cols$clusters_cols)
plot(color_col_dend)
make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, k = NULL, NULL, color_col_dend)


go_list_cols <- list()
for (cluster in clusters_cols){
  go_enrichment <- enrich_go(cluster$genes)
  go_list_cols <- append(go_list_cols, go_enrichment)
}
names(go_list_cols) <- c("1", "2", "3", "4", "5")
dotplot(go_list_cols$`1`)
dotplot(go_list_cols$`2`)
dotplot(go_list_cols$`3`)
dotplot(go_list_cols$`4`)
dotplot(go_list_cols$`5`)

##-----------------------------go annotation
# GO overrepresentation for each 6 clusters
go_list <- list()
for (cluster in clusters_rows){
  go_enrichment <- enrich_go(cluster$genes)
  go_list <- append(go_list, go_enrichment)
}

names(go_list) <- c("1", "2", "3", "4", "5", "6")
dotplot(go_list$`1`)
dotplot(go_list$`2`)
dotplot(go_list$`3`)
dotplot(go_list$`4`)
dotplot(go_list$`5`)
dotplot(go_list$`6`)


#make_heatmap(DE_genes, vsd, vsd$Patient_ID, NULL, k = NULL, test_anno)
#make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, k = 4)
################################################################################
###                               PCA                                        ###
################################################################################
# plot pca with batch correction
plot_pca(dds_default = dds_default, batch = vsd$Patient_ID)
# plot oca without batch correction
plot_pca(dds_default = dds_default, batch = NULL)



# Inspection--------------------------------------
#par(mar=c(8,5,2,2))
df <- log10(assays(dds_default)[["cooks"]])
names_ <- data.frame(paste0(metadata$Patient_ID, "-", metadata$Tumor_region))
colnames(df) <- c(names_$paste0.metadata.Patient_ID.......metadata.Tumor_region.)
df <- as.data.frame(df)
boxplot(df, range=0, las=2, cex.axis =0.7, ylab = "Log10 of Cook's distances")







res <- results(dds_default)
plotDispEsts(dds_default)
summary(res)
# TODO: make a function out of it  ---------------------------------------------
sampleDists <- dist(t(assay(vsd)))


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Tumor_region, vsd$Patient_ID, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap::pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)


################################################################################
###                             Combine exome                                ###
################################################################################


#metascape 





# test for enrichment analyiss and stuff----------------------------------------
# ensembla and org db sind nicht vergleichbar und deshlab lieber nicht mischen
# muss hier ne gute lösung finden

################################################################################
###                             GSEA                 #                        ##
################################################################################
# library("clusterProfiler")
# library("enrichplot")
# library("org.Hs.eg.db")
# 
# keytypes(org.Hs.eg.db)
# test2 <- data.frame(DE_list$`BEN vs NEC`)[ !is.na(data.frame(DE_list$`BEN vs NEC`$diffexp)),]
# test2 <- test2[c("Names", "log2FoldChange")]
# test2_2 <- unique(merge(test2,test[c("SYMBOL", "ENTREZID")], by.x = "Names", by.y = "SYMBOL", all.x = TRUE))
# test2_list <- test2$log2FoldChange
# test2_2list <- test2_2$log2FoldChange
# names(test2_2list) <- test2_2$ENTREZID
# test2_2sorted <- sort(test2_2list, decreasing = TRUE)
# gse <- gseGO(geneList = test_sorted,
#              ont = "ALL",
#              keyType = "SYMBOL",
#              minGSSize = 3, 
#              maxGSSize = 800, 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = org.Hs.eg.db, 
#              )
# #-------------------------------------------------------------------------------
# 
enrich_go <- function(gene_list) {
  go_enrichment <- enrichGO(gene = gene_list,
                            OrgDb  = org.Hs.eg.db,
                            ont     = "BP",
                            keyType = "SYMBOL",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.05)
  return(go_enrichment)
}
cluster2 <-  enrichGO(gene = clusters_rows$`2`$genes,
                 OrgDb  = org.Hs.eg.db,
                 ont     = "BP",
                 keyType = "SYMBOL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

goplot(cluster2)
dotplot(cluster2)
# #-------------------------------------------------------------------------------
# ggo2 <- groupGO(gene     = DE_genes,
#                OrgDb    = org.Hs.eg.db,
#                ont      = "MF",
#                keyType = "SYMBOL",
#                level    = 5)
# 
# head(ggo)
# 
# unlist(strsplit(ggo@result["GO:0016477", "geneID"], "/"))
# 
# 
# test_anno <- unique(data.frame("Genes" = DE_genes))
# test_anno$ann <- NA
# #cell migartion
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0016477", "geneID"], "/")),"cell migration", test_anno$ann)
# #	apoptotic process
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0006915", "geneID"], "/")),"apoptotic process", test_anno$ann)
# #	blood vessel morphogenesis
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0048514", "geneID"], "/")),"blood vessel morphogenesiss", test_anno$ann)
# #	regulation of cell differentiation
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0045595", "geneID"], "/")),"regulation of cell differentiation", test_anno$ann)
# 
# # angiogenesis
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0048514", "geneID"], "/")),"angiogenesis", test_anno$ann)
# #	T cell activation
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0042110", "geneID"], "/")),"T cell activation", test_anno$ann)
# 
# 
# rownames(test_anno) <- test_anno$Genes
# test_anno$Genes <- NULL
# 
# #inflammatory response
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0006954", "geneID"], "/")),"inflammatory response",test_anno$ann)
# # cellular response to hypoxia
# test_anno$ann <- ifelse(test_anno$Genes %in% unlist(strsplit(ggo@result["GO:0071456", "geneID"], "/")),"cellular response to hypoxia",test_anno$ann)
# test_anno$ann <- ifelse(is.na(test_anno$ann), "NO",test_anno$ann )
# 
# 
# #---------------------------------------------KEGG------------------------------
# 
# human <- search_kegg_organism("Homo sapiens", by = "scientific_name")
# 
# 
# kk <- enrichKEGG(gene = test1, organism = "hsa", pvalueCutoff = 0.05)
# 
# 
# test <- AnnotationDbi::select(org.Hs.eg.db, keys=DE_genes, columns=c("SYMBOL","ENTREZID","ENSEMBL"), keytype="SYMBOL") 
# test1 <- test$ENTREZID
# names(test1) <- test$SYMBOL
# 
# test1 <- test1[!is.na(test1)]
# 
# 
# hsa04668 <- pathview(gene.data  = test1,
#          pathway.id = "	hsa04668",
#          species    = "hsa",
#          limit      = list(gene=max(abs(geneList)), cpd=1))
# 
# 
# 
# 
# library(ReactomePA)
# 
# test2_2sorted <- test2_2sorted[!is.na(test2_2sorted)]
# x <- enrichPathway(gene=test1, pvalueCutoff = 0.05, readable=TRUE)
# 
# xx <- gsePathway(test2_2sorted, 
#                  pvalueCutoff = 0.2,
#                  pAdjustMethod = "BH", 
#                  verbose = FALSE)
#   
#   
# viewPathway("Signaling by PDGF", 
#               readable = TRUE, 
#               foldChange = test2_2sorted)  
# 
# #.------------------------------------------------------------------------------
# 
# library(DOSE)
# x <- enrichDO(gene          = test2_2sorted,
#               ont           = "DO",
#               pvalueCutoff  = 0.05,
#               pAdjustMethod = "BH",
#               universe      = names(geneList),
#               minGSSize     = 5,
#               maxGSSize     = 500,
#               qvalueCutoff  = 0.05,
#               readable      = FALSE)
# 
# ################################################################################
# ###                               Test Consensus clustering                  ###
# ################################################################################
# 
# library("M3C")
# 
# # remove batch 
# #filter significant genes 
# #remove outliers
# 
# 
# 
