################################################################################
###                                Purpose                                   ###
################################################################################
## Title: Functions for RNAseq Analysis 
## Author: Carolin Schwitalla
## 
## Description: This script contains all functions that are used for the RNAseq
##              analysis for my master thesis 
################################################################################
###                                                                          ###
################################################################################

# create tx2gene df that maps transcripts to genenames for tximport
# according to DESeq2 documentation
create_tx2gene <- function(){
  edb <- EnsDb.Hsapiens.v86
  k <- keys(edb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(edb, k, "GENENAME", "TXNAME")
  return(tx2gene)
}


## DESCRIPTION: function that performs DE analysis
## PARAMETERS:
#             - deseq2_object: DESeq(dds) 
#             - contrast: vector of strings, containing condition column 
#                         and 2 particular levels that should be compared
#             - lfc: numeric, threshold for log2FoldChange
#             - padj: numeric, thershold for p adjusted value
##
DE_analysis <- function(deseq2_object, contrast, lfc, padj){
  # get shrunken results of pairwise comparison
  # shrinkage reduces noisy lfcs for low count genes
  # contrast is used because of more than 2 conditions --> 
  # need to specify which conditions will be compared
  # ashr type of shrinkage because by using contrast it's the only option 
  DE_df <- lfcShrink(deseq2_object, contrast = contrast, type = "ashr")
  # include columns for UP , DOWN , SIG , expressed genes
  DE_df$significant <- "NO"
  DE_df$significant[abs(DE_df$log2FoldChange) > lfc & DE_df$padj < padj] <- "SIG"
  DE_df$diffexp[DE_df$log2FoldChange > 0 & DE_df$significant == "SIG"] <- "UP"
  DE_df$diffexp[DE_df$log2FoldChange < 0 & DE_df$significant == "SIG"] <- "DOWN"
  DE_df$Names <- DE_df@rownames
  return(DE_df)
  
}

## DESCRIPTION: function that makes Volcano plot
## PARAMETERS:
#             - res_df: DE dataframe from pairwise comparison 
#             - y_val: string --> either padju or pval
#             - title: string --> title of volcano plot 
#             - lfc: numeric, log2FoldChange threshold
#             - padj: numeric, thershold for p adjusted value
##
makeVolcano <- function(res_df, y_val, title, lfc, padju){
  EnhancedVolcano(res_df, lab = rownames(res_df),
                  x = "log2FoldChange",
                  y = y_val,
                  FCcutoff = lfc, 
                  pCutoff = padju,
                  ylim = c(0, max(log10(res_df$padj))),
                  xlim = c(min(res_df$log2FoldChange), 
                           max(res_df$log2FoldChange)),
                  ylab =  bquote(paste(-Log[10],.(y_val))),
                  selectLab = NA, # if TRUE: gene names will be visible as labels
                  labSize = 3, # gene label size
                  drawConnectors = TRUE, # gene and corresponding labels are 
                                         # connected with arrows
                  title = title,
                  legendPosition = "bottom",
                  legendLabels = c("NS", "Log2 FC",
                                   y_val,
                                   paste(y_val,"& Log2 FC")))
}


## DESCRIPTION: function that plots heatmap for a geneselection
## PARAMETERS:
#             - gene_selection: vecotr of genenames that will be plotted
#             - vsd: transformed count data object from DESeq2
##
makeHeatmap <- function(gene_selection ,vsd ){
  # get counts df for gene selection
  vsd_selection <- as.data.frame(assay(vsd)) %>% filter(row.names(assay(vsd)) %in% gene_selection)
  # scaling on rows, scale as default scales on columns, 
  #for rows we need transpose , and then transpose again to get original matrix 
  df_scale <- t(scale(t(vsd_selection)))
  # annotation for heatmap
  anno_df <- data.frame("Tumor_region" = vsd@colData@listData$Tumor_region,
                        "Patient" = vsd@colData@listData$Patient_ID,
                        "Sex" = vsd@colData@listData$Sex,
                        "MGMT_methylation" = vsd@colData@listData$MGMT)
  # annotation colors
  patientcolor = c("#543005", "#8c510a", "#A6691C", "#bf812d", "#CFA255", "#dfc27d",
                   "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#5BB2A8", "#35978f",
                   "#1B7F77", "#01665e", "#003c30")
  names(patientcolor2) = unique(c(vsd@colData@listData$Patient_ID))
  regioncolor = c("#9e0142", "#fdae61", "#74add1", "#4B6C22") 
  names(regioncolor) = c("NEC", "T1", "INF", "BEN")
  sexcolor = c("#565E71","#D5BD9E")
  names(sexcolor) = c("M", "F")
  methcolor = c("#74add1", "#4B6C22","#333A49")
  names(methcolor) = c("Unmethylated", "Methylated", "methylated")
  anno_colors = list(Tumor_region = regioncolor,
                           Sex = sexcolor, MGMT_methylation = methcolor, 
                           Patient =  patientcolor2)
  #draw heatmap
  pheatmap(df_scale, show_rownames = FALSE, annotation_col = anno_df, annotation_colors = anno_colors)
}

