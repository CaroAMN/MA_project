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
create_tx2gene <- function(){
  edb <- EnsDb.Hsapiens.v86
  k <- keys(edb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(edb, k, "GENENAME", "TXNAME")
  name2ID <- AnnotationDbi::select(edb,k, "GENEID", "TXNAME")
  mapping_df <- merge(tx2gene, name2ID, by = c("TXNAME","TXID"))
  return(mapping_df)
}


# diff expression analysis 
DE_analysis <- function(deseq2_object, contrast, lfc, padj){
  # get shrunken results of pairwise comparison
  DE_df <- lfcShrink(deseq2_object, contrast = contrast, type = "ashr")
  # include columns for UP , DOWN , SIG , expressed genes
  DE_df$significant <- "NO"
  DE_df$significant[abs(DE_df$log2FoldChange) > lfc & DE_df$padj < padj] <- "SIG"
  DE_df$diffexp[DE_df$log2FoldChange > 0 & DE_df$significant == "SIG"] <- "UP"
  DE_df$diffexp[DE_df$log2FoldChange < 0 & DE_df$significant == "SIG"] <- "DOWN"
  DE_df$Names <- DE_df@rownames
  return(DE_df)
  
}

# make volcano
makeVolcano <- function(res_df, y_val, title){
  EnhancedVolcano(res_df, lab = rownames(res_df),
                  x = "log2FoldChange",
                  y = y_val,
                  FCcutoff = 3,
                  pCutoff = 0.01,
                  ylim = c(0, max(log10(res_df$padj))),
                  xlim = c(min(res_df$log2FoldChange), 
                           max(res_df$log2FoldChange)),
                  ylab =  bquote(paste(-Log[10],.(y_val))),
                  selectLab = NA,
                  labSize = 3,
                  drawConnectors = TRUE,
                  title = title,
                  legendPosition = "bottom",
                  legendLabels = c("NS", "Log2 FC",
                                   y_val,
                                   paste(y_val,"& Log2 FC")))
  
}



