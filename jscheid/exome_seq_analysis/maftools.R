# ------------------------
#         PURPOSE
# ------------------------
# This script reads im maf files using maftools and
# plots different figures of exome seq data from different
# perspectives.
# maftool tutorial: 
# http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#13_MultiAssayExperiment
setwd("~/Documents/Uni/Master/Thesis/DNA/maftools")

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")
install.packages('R.utils')
library(maftools)


gbm = read.maf(maf = "all_vcfs_filtered.maf.gz")
prim <- read.maf(maf = "primary_vcfs_filtered.maf.gz")
rec <- read.maf(maf = "recurrent_vcfs_filtered.maf.gz")

# Plot maftools summary
plotmafSummary(maf = gbm, rmOutlier = T, addStat = 'median', dashboard = T, titvRaw = F)

# ----- from Rike
# Calculate Varaint allele frequency
compute_vaf <- function(maf_obj){
  genes_vaf = subsetMaf(maf = maf_obj, query = "t_depth > 0", fields = c("t_depth","t_alt_count"), mafObj = FALSE)
  genes_vaf$VAF <- genes_vaf$t_alt_count/genes_vaf$t_depth
  genes_vaf <- genes_vaf[, mean(VAF), Hugo_Symbol]
  colnames(genes_vaf)[2] = "VAF"
  genes_vaf <- genes_vaf[which(genes_vaf$VAF > 0.05)]
  return(genes_vaf)
}
# -----

# Oncoplot of the top 20 mutated genes
# TODO: Adjust margin of gbm oncoplot
# TODO Color Primary & Recurrent with colorbar: showTumorSampleBarcodes = T 8 etc
oncoplot(maf = gbm, top = 50, legendFontSize = 1.2,fontSize = 0.5, leftBarData = compute_vaf(gbm), logColBar = T,
         draw_titv = T, showTitle = F)
oncoplot(maf = prim, top = 20, legendFontSize = 1.2, leftBarData = compute_vaf(prim), logColBar = T,
         draw_titv = T, showTitle = F)

oncoplot(maf = rec, top = 20, legendFontSize = 1.2, leftBarData = compute_vaf(rec), logColBar = T,
         draw_titv = T, showTitle = F)

# Plot transitions and transversions
prim.titv = titv(maf = prim, plot = FALSE, useSyn = TRUE)
plotTiTv(res = prim.titv)

# Rainfallplot highlighting somatic hypermutation of one sample
sample <- read.maf(maf = "filteredStrelka_QLFGB016AJ_R_vs_QLFGB022AS_N_somatic_snvs.maf.gz")
rainfallPlot(maf = sample, tsb='QLFGB016AJ', detectChangePoints = TRUE, pointSize = 0.35)
# -> only one sample?

# Compate against TCGA DB
gbm.mutload = tcgaCompare(maf = gbm, cohortName = 'GBM_P+R', logscale = TRUE, capture_size = 50)
gbm.mutload = tcgaCompare(maf = prim, cohortName = 'GBM_Primary', logscale = TRUE, capture_size = 50)
gbm.mutload = tcgaCompare(maf = rec, cohortName = 'GBM_Recurrent', logscale = TRUE, capture_size = 50)

# Somatic interactions
somaticInteractions(maf = gbm, top = 20, pvalue = c(0.05, 0.1))
somaticInteractions(maf = prim, top = 20, pvalue = c(0.05, 0.1))
somaticInteractions(maf = rec, top = 20, pvalue = c(0.05, 0.1))

# Identify cancer driver genes
prim.driver = oncodrive(maf = prim, AACol = NULL, minMut = 5, pvalMethod = 'zscore') # ! CARE NOT 100% similar to website
plotOncodrive(res = prim.driver, fdrCutOff = 0.1, useFraction = T, labelSize = 0.6)

# Co-oncoplot with primary and recurrent gbm vcf files
# mafCompare  performs fisher test on all genes between two cohorts
# to detect differentially mutated genes
p_v_r <- mafCompare(m1 = prim, m2 = rec, m1Name = 'Primary', m2Name = 'Recurrent', minMut =5)
top10 <- p_v_r$results$Hugo_Symbol[1:10]
top20 <- p_v_r$results$Hugo_Symbol[1:20]
top30 <- p_v_r$results$Hugo_Symbol[1:30]
coOncoplot(m1 = prim, m2 = rec, m1Name = 'Primary', m2Name = 'Recurrent', genes = top30, removeNonMutated = TRUE, legendFontSize = 1.5)
coBarplot(m1 = prim, m2 = rec, m1Name = "Primary", m2Name = "Recurrent", genes = top20)
forestPlot(mafCompareRes = p_v_r, pVal = 0.05)

# OncogenicPathways function checks for enrichment of known Oncogenic
# Signaling Pathways in TCGA cohorts
OncogenicPathways(maf = prim)
OncogenicPathways(maf = rec)
OncogenicPathways(maf = gbm)

PlotOncogenicPathways(maf = prim, pathways = "RTK-RAS")


