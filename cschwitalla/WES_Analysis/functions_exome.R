################################################################################
###                                Purpose                                   ###
################################################################################
## Title: Functions for WES Analysis
## Author: Carolin Schwitalla
##
## Description: This script contains all functions that are used for the WES
##              analysis for my master thesis
################################################################################
###                                                                          ###
################################################################################

patient_color <- setNames(
  c(
    "#543005", "#8c510a", "#A6691C", "#bf812d", "#CFA255",
    "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1",
    "#5BB2A8", "#35978f", "#1B7F77", "#01665e", "#003c30"
  ),
  c(
    "ZH0984", "ZH0991", "ZH0997", "ZH1002", "ZH1006", "ZH1007",
    "ZH1008", "ZH1015", "ZH1019", "ZH1022", "ZH1028", "ZH1032",
    "ZH1033", "ZH1039", "ZH1041"
  )
)

region_color <- setNames(
  c("#9e0142", "#fdae61", "#74add1", "#4B6C22"),
  c("NEC", "T1", "INF", "BEN")
)

# make list of tumor region and patient colors for annotation in oncoplot
annotation_color <- list(Tumor_region = region_color, Patient_ID = patient_color)
# colors for oncoplot itself -->
vc_cols <- setNames(
  c("#A00202", "#438EA9", "#436211", "#323949", "#EE6D17", "#FFC61A"),
  c(
    "Frame_Shift_Del", "Missense_Mutation", "Nonsense_Mutation",
    "Multi_Hit", "Frame_Shift_Ins", "In_Frame_Del"
  )
)




################################################################################
###                                Functions                                 ###
################################################################################

# Variant allel Frequency
# ----- from Rike
# Calculate Varaint allele frequency
compute_vaf <- function(maf_obj) {
  genes_vaf <- subsetMaf(
    maf = maf_obj, query = "t_depth > 0",
    fields = c("t_depth", "t_alt_count"), mafObj = FALSE
  )
  genes_vaf$VAF <- as.numeric(genes_vaf$t_alt_count) / as.numeric(genes_vaf$t_depth)
  genes_vaf <- genes_vaf[, mean(VAF), Hugo_Symbol]
  colnames(genes_vaf)[2] <- "VAF"
  genes_vaf <- genes_vaf[which(genes_vaf$VAF > 0.02)]
  return(genes_vaf)
}
