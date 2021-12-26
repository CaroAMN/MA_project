# -----------------------------------
# ----------- PURPOSE ---------------
# -----------------------------------
# Obtain the pre-processed TPM values, processed are derived as follows in Nelde et al. 2019:
# Expression levels were quantified as Ei,j = log2(TPMi,j/10+1), where TPMi,j refers to transcript-per-million
# for gene i in sample j, as calculated by RSEM (Li and Dewey, 2011). TPM values were divided by 10 since we 
# estimate the complexity of single cell libraries in the order of 100,000 transcripts and would like to avoid
# counting each transcript ∼10 times, as would be the case with TPM, which may inflate the difference between
# the expression level of a gene in cells in which the gene is detected and those in which it is not detected
# -----------------------------------
# Additionally tailor the data for the Input of Scaden and Cibersort according to the 4 subtypes Nelde et al proposed
# -----------------------------------

ss2log2TPM <- read.csv("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv", sep = "\t")
celltypes <- read.csv("neftel_all_celltypes.txt", sep = '\t')
# R exchanges some dashes with "neg" or dots
celltypes$X <- gsub("neg","-", celltypes$X)
colnames(ss2log2TPM) <- gsub("neg","-", colnames(ss2log2TPM))
colnames(ss2log2TPM) <- gsub("\\.","-", colnames(ss2log2TPM))
rownames(ss2log2TPM) <- ss2log2TPM$GENE
ss2log2TPM$GENE <- NULL
# Only obtain the 4 glioblastoma subtypes
subtypes <- celltypes[which(is.na(match(celltypes$Celltype, c("AClike","MESlike","NPClike","OPClike"))) == F),]
write.table(subtypes,"neftel_celltypes.txt", quote = F, sep = '\t', row.names = F)
ss2log2TPM <- ss2log2TPM[,match(subtypes$X, colnames(ss2log2TPM))]
# To obtain the TPM compute (2^(TPM) - 1) * 10
ss2TPM <- (2^(ss2log2TPM) - 1) * 10

# -----------------------------------
# ------------ SCADEN ---------------
# -----------------------------------
# Rearrange dfs for Scaden input => Genes as columns, samples/cells as rows
t_ss2TPM <- data.frame(t(ss2TPM))
# Rearrange matrix according to celltype annotation
t_ss2TPM <- t_ss2TPM[match(subtypes$X, rownames(t_ss2TPM)),]
# Adjust according to example scaden input
rownames(t_ss2TPM) <- seq(length(rownames(t_ss2TPM)))
X <- seq(length(rownames(t_ss2TPM))) - 1
t_ss2TPM <- cbind(X,t_ss2TPM)
write.table(t_ss2TPM,"neftel_counts.txt", sep = "\t", quote = F)
# Adjust bulkRNA to Scaden input
bulkRNA <- read.table("../../RNA/gene_counts_TPM.txt", sep = '\t', header=T)
write.table(bulkRNA[-1], "../scaden/neftel_model/bulkRNA_TPM.txt", quote= F, row.names = F ,sep = '\t')

# -----------------------------------
# ---------- CIBERSORTx -------------
# -----------------------------------
# Step 1: Prepare single cell reference sample file
# All single cells must be pre-labeled with the cell’s phenotype or cluster identifier
# colnames need to be included into first row, because colnames need to be unique
Samples <- colnames(ss2TPM)
scRef <- rbind(Samples, ss2TPM)
Gene <- c("Gene",rownames(ss2TPM))
scRef <- cbind(Gene, scRef)
# Exchange samples with celltype annotation of Nelde et al. 2019
scRef[1,2:length(colnames(scRef))] <- subtypes$Celltype 
write.table(scRef, "../cibersort/single_cell_reference_cibersort.txt", quote = F, sep = '\t', row.names = F, col.names = F)
# Prepare bulkRNA data for cibersort input -> only remove ensemble column
bulkRNA <- read.table("../../RNA/gene_counts_TPM.txt", sep = '\t', header=T)
write.table(bulkRNA[-1],"../cibersort/bulkRNA_cibersort.txt", quote = F, sep = '\t', row.names = F)
