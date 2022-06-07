################################################################################
###                                Purpose                                   ###
################################################################################
## DESCRIPTIOM: This script performs differential expression analysis on Salmon
##              quantification files, based on the DESeq2 package.
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## 
## 
## INPUT: Salmon quantification files 
## OUTPUT: for each pairwise comparison: - count table with DE expression
##                                         Annotation --> SIG, DOWN, UP
##                                       - Volcano plot
##                                       - MA plot
##
################################################################################
###                             Manual Steps                                 ###
################################################################################
#TODO workdir and write table directory are the same --> write tables into output dir
# setwd will be unused when file names contain komplete path to files for imort
# clear R enviroment
rm(list = ls())
# set working dir
setwd("/Users/cschwitalla/Documents/transcriptomics_results/Pipeline_running/transcriptomics_results/Quant_files/")
# set directory of the input files for the analysis 
inputDir <- "/Users/cschwitalla/Documents/transcriptomics_results/Pipeline_running/transcriptomics_results/Quant_files/"
# set output directory for results table and plots
outputDir <- "/Users/cschwitalla/Documents/RNAseq_analysis/"
# set directory of the metadata file
metadataFile <- "./Metadata_GB.tsv"
# set directory of the script that contains RNAseq analysis functions
source("/Users/cschwitalla/git/students/cschwitalla/RNAseq_Analysis/functions_RNA.R")
# set lfc and padju values for significance thresholds
lfc <- 3
padju <- 0.01
################################################################################
###                           Load libraries                                 ###
################################################################################
required_Libs <- c("DESeq2", "tximport", "ashr", "ggVennDiagram", "ggplot2",
                   "dplyr", "readxl", "tidyr","pheatmap", "RColorBrewer",
                   "EnsDb.Hsapiens.v86","ggrepel", "utils", "EnhancedVolcano")

# invisible(lapply(necessaryLibs, BiocManager::install, update = F, ask = F))
suppressMessages(invisible(lapply(required_Libs, library, character.only = T)))


################################################################################
###                            Load Data                                     ###
################################################################################
# Load meta data --> Metadata_GB.tsv in workdir
metadata <- read.table(file = metadataFile ,sep = "\t", header = TRUE)
# get filenames of inputDir
file_names  <- list.files(path = inputDir)
# TO DO: create/ import tx2gene dataframe
tx2gene <- create_tx2gene()
# import Salmon quant files with tximport and tx2gene
txi <- tximport(file_names, type = "salmon", tx2gene = tx2gene)
# make Deseq2 object 
dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata, 
                                design = ~ Tumor_region)
# TO DO: why? --> pre filtering, keep only rows with at least 10 reads total
# WHY: 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# TO DO: what doese this step 
dds <- DESeq(dds)
################################################################################
###                            DE_Analysis                                   ###
################################################################################
# list of genes that are diff expressed
DE_genes <- c()
# get each pairwise comparison for condition that should be tested
for (com in apply(combn(levels(dds$Tumor_region), 2), 2, paste, collapse = " vs ")) {
  # get condition vector for DE analysis function
  comparison <- unlist(strsplit(com, " vs "))
  # make results data frame
  res_df <- DE_analysis(dds, c("Tumor_region", comparison), lfc, padju)
  # get DE genes and append them into DE_genes list 
  DE_genes <- append(DE_genes, res_df$Names[res_df$significant =="SIG"])
  # write results as tsv into the workdir to stor results
  write.table(res_df, file = paste(com,".tsv" ), sep = "\t", quote = FALSE)
  # make MA plot and save it as png in output dir 
  png(paste0(outputDir,"MAplot_", com, ".png")) 
  print({plotMA(res_df, xlim = c(1,1e5), ylim = c(-4,4), main = com)})
  dev.off()
  # make Volcano plot and save it as png in output dir 
  png(paste0(outputDir, "Volcano_", com, ".png"))
  print({makeVolcano(res_df, "padj", com, lfc, padju)})
  dev.off()
}



################################################################################
###                              Heatmap                                     ###
################################################################################
# vsd transformation of the data --> why ??
vsd <- vst(dds, blind = FALSE)
# plot heatmap
makeHeatmap(DE_genes,vsd)




