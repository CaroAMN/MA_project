#---------------------------------------------------------------------
#                            Purpose                                
#---------------------------------------------------------------------
# This script aims to calculate the raw counts of a expression matrix
# into transcripts-per-million (TPM) values

#---------------------------------------------------------------------
#                            Preperation                             
#---------------------------------------------------------------------
# We need to know the gene lengths to calculate the TPM values
# Gene lengths can be obtained by a reference genome, e.g. using:
# awk -F "\t" '$3=="transcript"{print $9, $5-$4} hg38.gtf > transcripts.txt'
# Then genes and gene lengths need to be extracted
#   -> Did this in bash
# ==> there are some ambigouties with additional columns which will 
#     then be solved
#------------------------------------------------------------------
# Dependencies
#install.packages("/Users/jonas/Downloads/hugo_1.0.0.tar.gz", repos=NULL, type="source")

setwd("~/Documents/Uni/Master/Thesis/RNA")
library("hugo")
library("data.table")

# Read gene lengths file
geneLen = read.csv("hg38_len.txt", header = FALSE, sep=';')

# Declare colnames and correct geneID col
geneLen$V1 <- NULL
colnames(geneLen) <- c("geneID", "geneLength")
geneLen$geneID = gsub(' ','', geneLen$geneID)

# Read count matrix
counts = read.csv("RNAseq_results_June/differential_gene_expression/gene_counts_tables/raw_gene_counts.tsv", sep = '\t')
# Rearrange indices
rownames(counts) <- counts$Ensembl_ID
counts$Ensembl_ID <- NULL

# Load gene annotation db HUGO
db <- hugo.load()
# Get HUGO annotations of ensembl ids
counts$hugoGenes <- hugo.convert(db, rownames(counts), from = "ensembl", to="symbol",ignore.case = TRUE, multiple= "collapse")
geneLen$hugoGenes <- hugo.convert(db, geneLen$geneID, from = "ensembl", to="symbol",ignore.case = TRUE, multiple= "collapse")

#--------------------------
#       Quick Stats
#-------------------------
# Mapping Hugo annotations to raw count matrix and genome and then matching their genes obtaines more genes
# sum(is.na(match(rownames(counts),geneLen$geneID))) = 4502
# sum(is.na(match(counts$hugoGenes,geneLen$hugoGenes))) = 138
# length(unique(counts$hugoGenes)) = 34247
# length(unique(geneLen$hugoGenes)) = 37302
# sum(duplicated(counts$geneName)) = 1700
#--------------------------

# Extract only the hugo genes
hugoCounts <- counts[which(is.na(match(counts$hugoGenes, geneLen$hugoGenes)) == F), ]
hugoCounts$geneName <- hugoCounts$hugoGenes
hugoCounts$hugoGenes <- NULL
# Get mean transcript lengths of each gene
meanTranscriptLength <- aggregate(geneLen$geneLength, list(geneLen$hugoGenes), mean)
colnames(meanTranscriptLength) <- c("geneName","geneLength")
# Add gene length to counts
hugoCounts$geneLength <- meanTranscriptLength$geneLength[match(hugoCounts$geneName, meanTranscriptLength$geneName)]
# Remove NA's in gene columns => only work with hugo genes
hugoCounts <- na.omit(hugoCounts)

#-------------Calculate TPM-------------------
# Adapted from: https://gist.github.com/slowkow/c6ab0348747f86e2748b
# Step 1: Normalize raw counts for gene length (in kB)-> Compute reads per kilobase (RPK)

rpk <- hugoCounts[1]
for (i in 2:45) {
  rpk[colnames(hugoCounts)[i]] <- hugoCounts[i] / (hugoCounts$geneLength * 0.001) 
}

# Step 2: Normalize for sequencing depth
# -> Sum up all rpk values of a sample and use it as scaling factor for the sample

tpm <- hugoCounts[1]
for (i in 2:45) {
  tpm[colnames(hugoCounts)[i]] <- rpk[,i] / sum(rpk[,i]) * 1e6
}
#---------------------------------------------

write.table(tpm, file = "gene_counts_TPM.tsv", sep="\t", quote = F)

#-----
# Code reminder
View(counts[which(is.na(match(hugoCounts$hugoGenes, geneLen$hugoGenes)) == T),])
make.unique(e)
View(counts[duplicated(counts$geneName),])
length(unique(counts$geneName))
View(counts)

tpm <- read.csv("gene_counts_TPM.tsv", sep = "\t")
sum(log(tpm[,2:45])[,2])
