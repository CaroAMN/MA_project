################################################################################
###                                Purpose                                   ###
################################################################################
## DESCRIPTIOM: This script performs WES analysis by using mainly the maftools
##              package
##
## INPUT: maf files
## OUTPUT: Oncoplots
################################################################################
###                               Manual steps                               ###
################################################################################
## Clear the R environment
rm(list = ls())
# set working directory (top folder were the maf fiels are)
setwd("/Users/cschwitalla/Documents/WES_data/")
# list of directories were maffiles for each region are +
# set regionName as name in same order like dirs to know which maf file
# belogs to which region
region_dir <- setNames(
  c(
    "./necrotic/maf_converted/",
    "./t1/maf_converted/",
    "./peritumoral/maf_converted/",
    "./ben/PBMC_vs_BEN/Annotation/VEP/maf_converted/"
  ),
  c("NEC", "T1", "INF", "BEN")
)
# define output directory
output_dir <- "/Users/cschwitalla/Documents/WES_analysis/"
# set path to the metadata file
metadata_file <- "/Users/cschwitalla/Documents/WES_analysis/WES_metadata.tsv"
# set path to functions file required for this script
source("/Users/cschwitalla/git/students/cschwitalla/WES_Analysis/functions_exome.R")
################################################################################
###                             Load libraries                               ###
################################################################################
# Load the libraries
required_libs <- c("maftools", "ggplot2", "circlize")

suppressMessages(invisible(lapply(required_Libs, library, character.only = T)))

################################################################################
###                            Load Data                                     ###
################################################################################
# load metadata file
metadata <- read.table(file = metadata_file, sep = "\t", header = TRUE)



# TODO: bild größe anpassen
# TODO: merge all mafs + make oncoplot for all mafs
# list to store the merged maf files for each region
regionmaf_list <- list()
for (dir in region_dir) {
  # get file paths
  files <- list.files(dir, full.names = TRUE, pattern = ".maf")
  # merge maf files + metadata
  maf <- merge_mafs(files, clinicalData = metadata)
  # append maf file into list
  regionmaf_list <- append(regionmaf_list, maf)
  # make png
  png(paste0(output_dir, names(region_dir)[match(dir, region_dir)], "_oncoplot.png"))
  # make oncoplot for current maf that will exported as png in the output dir
  print(oncoplot(
    maf = maf,
    draw_titv = FALSE,
    clinicalFeatures = c("Patient_ID", "Tumor_region"),
    leftBarData = compute_vaf(maf), # see function script
    annotationColor = annotation_color, # see function script
    colors = vc_cols
  )) # see function script
  dev.off()
}
# set names for merged maf files in the list, to see which region is it
names(regionmaf_list) <- c("NEC", "T1", "INF", "BEN")
# make venn diagram including all genes from all regions
# merge all tumor region mafs + metadata
# make oncoplot with all mafs
# -- gene selection befor plotting


# include lollipopplot maybe
# circouise plots
circos.initializeWithIdeogram(species = "hg19")

# braiche start und end position + welche mutation
