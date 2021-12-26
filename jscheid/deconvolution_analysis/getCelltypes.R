############################
######## Purpose ###########
############################

# This script aims to prepare cell typed scRNA-seq data as a reference data set for using
# CIBERSORTx or Scaden in order to deconvolute the cell types in the bulk RNA-seq gbm data

##########################################################################################
# Data from Neftle et al 2019 (https://www.sciencedirect.com/science/article/pii/S0092867419306877# )
# Cell annotations obtained from the Broad single-cell portal (https://singlecell.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma)
# TPM matrix obtained from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928)
##########################################################################################


smrtsq2 <- read.csv("Neftel_et_al_2019/IDHwtGBM.processed.SS2.logTPM.txt", sep = '\t')
anno <- read.csv("Neftel_et_al_2019/cell_type_annotation.txt", sep = '\t')

# Remove first row and turn chars into numeric values
anno <- anno[-1,]
for (i in 8:15) {
  anno[,i] <- as.numeric(anno[,i])
}

# Scaden needs a cell type annotation file with 2 columns (cellID, Celltype)
# The annotation file contains the degree of a cell being in one of the subtypes
celltypes <- data.frame(row.names = anno$NAME)
celltypes[1] <- colnames(anno[8:15])[max.col(anno[,8:15])]
colnames(celltypes) <-'Celltype'
# NAs are benign entries -> also provided in metadata
for (i in 1:length(celltypes$Celltype)){
  if (is.na(celltypes$Celltype[i])){
    celltypes$Celltype[i] <- anno$CellAssignment[i]
  }
}
# Combine MESlike1 with MESlike2 and NPClike1 with NPClike2
celltypes$Celltype[grep("MESlike", celltypes$Celltype)] <- 'MESlike'
celltypes$Celltype[grep("NPClike", celltypes$Celltype)] <- 'NPClike'

write.table(celltypes, file = "Neftel_et_al_2019/neftle_celltypes.txt", sep = '\t' , quote=FALSE)



