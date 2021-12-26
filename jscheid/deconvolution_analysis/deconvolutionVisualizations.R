# -----------------------------
# --------- PURPOSE -----------
# -----------------------------
# This script creates figures related to the deconvolutet subset compostion
# of the bulkRNA seq glioblastoma data.
library(ggplot2)
library(RColorBrewer)
# Read in Scaden results
scaden_out <- read.table("scaden/scaden_predictions.txt", sep = '\t', header = T, row.names = 1)
rownames(scaden_out) <- gsub('.*_G','G', rownames(scaden_out))
# Order df
barplot(t(scaden_out))
# Melt data to use violin plot function
violin_data_scaden <- data.frame(cbind(rep(rownames(scaden_out), 4), c(scaden_out$MESlike, scaden_out$OPClike, scaden_out$AClike, scaden_out$NPClike), c(rep("MESlike", dim(scaden_out)[1]), rep("OPClike", dim(scaden_out)[1]), rep("AClike", dim(scaden_out)[1]), rep("NPClike", dim(scaden_out)[1]))))
colnames(violin_data_scaden) <- c("sample", "fraction", "annotation")
violin_data_scaden$fraction <- as.numeric(violin_data_scaden$fraction)
violin_data_scaden$sample <- gsub("\\.", "", violin_data_scaden$sample)
violin_data_scaden$sample <- as.factor(violin_data_scaden$sample)
violin_data_scaden$annotation <- as.factor(violin_data_scaden$annotation)

ggplot(violin_data_scaden, aes(x=annotation, y=fraction, fill=sample)) +
  ggplot2::geom_violin()

str(violin_data_scaden)
# ====> DOESNT WORK ATM

# Sort Df according to sample name
scaden_out <- scaden_out[order(rownames(scaden_out)),c(3,1,4,2)]
# Barplot Scaden
sc <- barplot(t(scaden_out), col=brewer.pal(4, "Pastel1"), xaxt='n', ylab='Cell Composition', xlim = c(0, 60),legend=T)
text(bp, par("usr")[3], labels = rownames(scaden_out), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.7)

# Read in cibersortx results 
cibersortx_out <- read.table("cibersort/cibersortx_output.txt", sep = '\t', header = T, row.names = 1)
rownames(cibersortx_out) <- gsub('.*_G','G', rownames(cibersortx_out))
# Sort Df according to sample name
cibersortx_out <- cibersortx_out[order(rownames(cibersortx_out)),]
# Barplot cibersort
cs <- barplot(t(cibersortx_out[1:4]), col=brewer.pal(4, "Pastel1"), xaxt='n', ylab='Cell Composition', xlim = c(0, 60),legend=T)
text(cs, par("usr")[3], labels = rownames(cibersortx_out), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.7)

# Circular Barplot?
# https://www.r-graph-gallery.com/299-circular-stacked-barplot.html

# Categorize samples into 4 subtypes

