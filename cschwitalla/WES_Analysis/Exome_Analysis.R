################################################################################
###                                Purpose                                   ###
################################################################################
## DESCRIPTIOM: This script performs WES analysis by using mainly the maftools
##              package 
## 
## INPUT: maffiles  
## OUTPUT: Oncoplots
################################################################################
###                               Manual steps                               ###
################################################################################
## Clear the R environment
rm(list = ls())
# set working directory (top folder were the maf fiels are)
setwd("/Users/cschwitalla/Documents/WES_data/")
# list of directories were maffiles for each region are + 
# set regionName as name in same order like dirs to know which maf fiels belogs to which region
regionDir <- setNames(c("./necrotic/maf_converted/", 
               "./t1/maf_converted/", 
               "./peritumoral/maf_converted/",
               "./ben/PBMC_vs_BEN/Annotation/VEP/maf_converted/"),
               c("NEC", "T1", "INF", "BEN") )
# define output directory
outputDir <- "/Users/cschwitalla/Documents/WES_analysis/"
# set path to the metadata file
metadataFile <- "/Users/cschwitalla/Documents/WES_analysis/WES_metadata.tsv"
# set path to functions file required for this script
source("/Users/cschwitalla/git/students/cschwitalla/WES_Analysis/functions_exome.R")
################################################################################
###                             Load libraries                               ###
################################################################################
# Load the libraries
required_Libs <- c("maftools", "ggplot2", "circlize")

# invisible(lapply(necessaryLibs, BiocManager::install, update = F, ask = F))
suppressMessages(invisible(lapply(required_Libs, library, character.only = T)))

################################################################################
###                            Load Data                                     ###
################################################################################
# load metadata file 
metadata <- read.table(file = metadataFile ,sep = "\t", header = TRUE)



# TODO: bild größe anpassen
# TODO: merge all mafs + make oncoplot for all mafs
# list to store the merged maf files for each region
regionMaf_list <- list()
for (dir in regionDir) {
  #get file paths
  files <- list.files(dir, full.names = TRUE, pattern = ".maf")
  # merge maf files + metadata
  maf <- merge_mafs(files, clinicalData = metadata)
  # append maf file into list 
  regionMaf_list  <- append(regionMaf_list ,maf)
  # make png 
  png(paste0(outputDir, names(regionDir)[match(dir,regionDir)], "_oncoplot.png"))
  # make oncoplot for current maf that will exported as png in the output dir
  print(oncoplot(maf = maf,
           draw_titv = FALSE,
           clinicalFeatures = c("Patient_ID", "Tumor_region"),
           leftBarData = compute_vaf(maf), # see function script
           annotationColor = annotation_color, # see function script 
           colors = vc_cols)) # see function script
  dev.off()
  
}
# set names for merged maf files in the list, to see which region is it 
names(regionMaf_list ) <- c("NEC", "T1", "INF", "BEN")
# make venn diagram including all genes from all regions 
# merge all tumor region mafs + metadata
# make oncoplot with all mafs
# -- gene selection befor plotting 


# include lollipopplot maybe 
# circouise plots
circos.initializeWithIdeogram(species = "hg19")

# braiche start und end position + welche mutation 
