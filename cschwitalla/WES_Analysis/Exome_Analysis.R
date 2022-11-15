################################################################################
###                                Purpose                                   ###
################################################################################
## DESCRIPTIOM: This script performs WES analysis by using mainly the maftools
##              package
##
## INPUT: maf files
## OUTPUT: Oncoplots
################################################################################
###                               TODOs                                      ###
################################################################################
# TODO: affected pathways from maftools 
# TODO: appropriate filtering for oncoplot and in general 
# TODO: known mutations in GB 

################################################################################
###                               Manual steps                               ###
################################################################################
## Clear the R environment
rm(list = ls())
# set working directory (top folder were the maf fiels are)
setwd("/Users/cschwitalla/Documents/WES_data/sarek3_annotation/strelka/")

# list of directories were maffiles for each region are +
# set regionName as name in same order like dirs to know which maf file
# belogs to which region

region_dir <- setNames(
  c(
    "./NEC//maf_converted/",
    "./T1/maf_converted/",
    "./INF/maf_converted/",
    "./BEN/maf_converted/"
  ),
  c("NEC", "T1", "INF", "BEN")
)
# define output directory
output_dir <- "/Users/cschwitalla/Documents/WES_analysis/"
# set path to the metadata file
metadata_file <- "/Users/cschwitalla/Documents/WES_analysis/WES_metadata.tsv"
# bed file for capture kit size
bed_file_dir <- "../../hglft_genome_298d9_1b60f0.sort.merged.bed"

# set path to functions file required for this script
source("/Users/cschwitalla/git/students/cschwitalla/WES_Analysis/functions_exome.R")
source("/Users/cschwitalla/git/students/cschwitalla/Ligandomics_Analysis/functions_ligandomics.R")
################################################################################
###                             Load libraries                               ###
################################################################################
# Load the libraries
required_libs <- c("maftools", "ggplot2", "circlize", "UpSetR", "dplyr", "tuple",
                   "ggVennDiagram", "stringr", "plyr")

suppressMessages(invisible(lapply(required_libs, library, character.only = T)))

################################################################################
###                            Load Data                                     ###
################################################################################

#test <- read.maf(paste0("/Users/cschwitalla/Documents/WES_data/sarek3_annotation/strelka/NEC/maf_converted/Strelka_QATLV067AW_vs_QATLV069AE_somatic_indels_VEP.ann_VEP.ann_filtered.vcf.maf"))

# load metadata file
metadata <- read.table(file = metadata_file, sep = "\t", header = TRUE)
#metadata <- metadata %>% arrange(factor(metadata$Tumor_region, levels = c("BEN", "INF", "T1", "NEC")))
#calc capture kit size 
bed_file <- read.csv(file = bed_file_dir, 
                     sep = "\t",
                     header = FALSE)
kit_size = sum(bed_file$V3 - bed_file$V2)

# list to store the merged maf files for each region
regionmaf_list <- list()
all_files <- c()
tumor_files <- c()
most_mutated_genes <- list()
sample_order <- unique(metadata$Tumor_Sample_Barcode)
# looooooop --------------------------------------------------------------------
for (dir in region_dir) {
  # get file paths
  files <- list.files(dir, full.names = TRUE, pattern = ".maf")
  # get all file directories for merging all mafs
  all_files <- append(all_files, files)
  # get only tumor region maffieles for merging 
  if (dir != "./BEN/maf_converted/") {
    tumor_files <- append(tumor_files,files)
  }
  # merge maf files + metadata
  maf <- merge_mafs(files, clinicalData = metadata)
  # get clinical significance 
  clin_sig_region <- maf@data[,c("Hugo_Symbol", "CLIN_SIG")] %>%
    dplyr::filter(!grepl("likely_benign|benign", CLIN_SIG) & CLIN_SIG != "")
  # get novel and unknown variants to dbSNP 
  novel_snvs <- maf@data[,c("Hugo_Symbol", "dbSNP_RS")] %>%
    dplyr::filter(grepl("novel", dbSNP_RS) | dbSNP_RS == "")
  novel_snvs_maf <- subsetMaf(maf, query = "dbSNP_RS == 'novel' | dbSNP_RS == ''")
  novel_snvs_symbols <- novel_snvs_maf@gene.summary[novel_snvs_maf@gene.summary$MutatedSamples >= 2, "Hugo_Symbol"]
  
  # most mutated genes
  if (dir == "./BEN/maf_converted/") {
    most_mutated_genes <- append(most_mutated_genes, list(maf@gene.summary$Hugo_Symbol[maf@gene.summary$MutatedSamples >= 2]))
  }else{
    most_mutated_genes <- append(most_mutated_genes, list(maf@gene.summary$Hugo_Symbol[maf@gene.summary$MutatedSamples >= 3]))
    
  }

  regionmaf_list <- append(regionmaf_list, maf)
  #make annotation list with metadata
  anno_color <- create_annotation_color(patient_column = metadata$Patient_ID,
                                        region_column = metadata$Tumor_region)
  # make pdf
  make_oncoplot(output_dir, region_dir, maf, anno_color, "_all", genes = NULL, sample_order = sample_order)
  # make ocoplot with only clinical sig genes
  make_oncoplot(output_dir, region_dir, maf, anno_color, "_ClinSig_genes", c(unique(clin_sig_region$Hugo_Symbol)),sample_order = sample_order)
  # make oncoplot with novel variants ( not fpind in dbSNP)
  if (length(novel_snvs_symbols$Hugo_Symbol) >= 2) {
    make_oncoplot(output_dir, region_dir, novel_snvs_maf, anno_color, "_novel_snvs", c(novel_snvs_symbols$Hugo_Symbol),sample_order = sample_order)
  }
  # write summary table -> number of genes mutated in number of samples 
  write.table(as.data.frame(table(maf@gene.summary$MutatedSamples)),
              paste0( output_dir, names(region_dir)[match(dir, region_dir)],"_sumary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  #plot maf summary
  plot_summary(output_dir,dir,maf)
  # mutload compare tcga
  tcga_comapre(maf, dir, kit_size, output_dir)
  #onco pathways
  onco_pathways(maf, output_dir, dir)
  #drug interactions
  drug_interactions(output_dir,dir,maf)
}

# set names for merged maf files in the list, to see which region is it
names(regionmaf_list) <- c("NEC", "T1", "INF", "BEN")
names(most_mutated_genes) <- c("NEC", "T1", "INF", "BEN")


# make venn diagram including all genes from all regions
all_genes_set <- setNames(vector("list", length = length(names(regionmaf_list))),
                          c(names(regionmaf_list))
)

for (i in names(regionmaf_list)) {
  all_genes_set[i] <- list(regionmaf_list[[i]]@data$Hugo_Symbol)
  #data <- merge(regionmaf_list[[i]]@data[,c("Hugo_Symbol","ENSP")], regionmaf_list[[i]]@gene)
  write.csv(regionmaf_list[[i]]@data[,c("Hugo_Symbol","ENSP","Tumor_Sample_Barcode")],
            paste0(output_dir, i, "_variants.tsv"),
            row.names = FALSE,
            quote = FALSE )
}



regionmaf_list[[i]]@data[,c("Hugo_Symbol","ENSP")]
test <- regionmaf_list$NEC
#plot custom venn for all genes
all_genes_venn_data <- plot_custom_venn(all_genes_set, "all genes")
#get region exclusive mutations
all_maf_lolipoplot <- subsetMaf()


# merge all tumor region mafs + metadata
all_maf <- merge_mafs(all_files, clinicalData = metadata)
# dataframe for lolipop plots
all_maf_lolipoplot <- dplyr::select(as.data.frame(all_maf@data),"Hugo_Symbol", "Tumor_Sample_Barcode","Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2","Allele", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "PolyPhen")
test <- dplyr::select(all_maf@data,"Tumor_Sample_Barcode", "Hugo_Symbol")
test2 <- merge(test, metadata[c("Tumor_Sample_Barcode", "Patient_ID")], by = "Tumor_Sample_Barcode", all.x = TRUE)
#test3 <- test2 %>% group_by(Patient_ID) %>% summarise(Hugo_Symbol = unique(Hugo_Symbol))
test3 <- test2[,c("Hugo_Symbol", "Patient_ID")]
test3 <- unique(test3)
test4 <- data.frame(table(test3$Hugo_Symbol))
table(test4$Freq)
total <- 1101 + 361 + 379 + 100 +  30 +  29 +   5 +   5   + 2   + 2   + 2 +   1   + 1   + 1   + 1   + 1 

# make oncoplot with all mafs
metadata <- metadata %>% arrange(factor(metadata$Tumor_region, levels = c("BEN", "INF", "T1", "NEC")))
sample_order <- unique(metadata$Tumor_Sample_Barcode)
# overall top 20 most mutated genes
vaf_all <- make_oncoplot(output_dir, "all", all_maf, anno_color, "all", genes = NULL, sample_order = sample_order, NULL, 6)
# most_mut <- all_maf@gene.summary 
# most_mut <- most_mut[order(most_mut$MutatedSamples, decreasing = TRUE),]
# most
most_mut_genes <- test4[order(test4$Freq,decreasing = TRUE),]
most_mut_genes <- most_mut_genes[which(most_mut_genes$Freq >= 3),]

most_mut_genes <- rbind(most_mut_genes, data.frame(list(Var1 = "Dummy1", Freq = 0)))
most_mut_genes <- rbind(most_mut_genes, data.frame(list(Var1 = "Dummy2", Freq = 15)))
# oncoplot for all most mutated genes for ech region
make_oncoplot(output_dir, "all", all_maf, anno_color, "most_mutated_all", genes = most_mut_genes$Var1, sample_order = sample_order)


mutload <- merge(all_maf@variants.per.sample, metadata, by= "Tumor_Sample_Barcode", all.x = TRUE)

#oncoplot with only exclusive genes
make_oncoplot(output_dir, "exclusive_genes", all_maf, anno_color, "exclusive genes",
              genes = c("TXNIP","DNAH8", "ANKRD30A", "PAX3", "OR2A1"),
              sample_order = sample_order)

# make mutational burden plot of all regions maff 
tumor_only_maff <- merge_mafs(tumor_files, clinicalData = metadata)
tcga_comapre(c(regionmaf_list$NEC, regionmaf_list$T1, regionmaf_list$INF, regionmaf_list$BEN),NULL, kit_size, output_dir, c("NEC", "T1", "INF", "BEN"))
tcga_comapre(tumor_only_maff, NULL, kit_size, output_dir, "tumor_regions")
t <- tcgaCompare(tumor_only_maff, capture_size = kit_size/1000000)
b <- tcgaCompare(regionmaf_list$BEN, capture_size = kit_size/1000000)
n <- tcgaCompare(regionmaf_list$NEC, capture_size = kit_size/1000000)
t1 <- tcgaCompare(regionmaf_list$T1, capture_size = kit_size/1000000)
inf <- tcgaCompare(regionmaf_list$INF, capture_size = kit_size/1000000)

# lollipop plots
Lollipop_plot("TXNIP", all_maf_lolipoplot)
Lollipop_plot("DNAH8", all_maf_lolipoplot)
Lollipop_plot("IER5L", all_maf_lolipoplot)
Lollipop_plot("NF1", all_maf_lolipoplot)

nf1_maf <- subsetMaf(all_maf, genes = "NF1")

pik3_maf <- subsetMaf(regionmaf_list$BEN, genes = "PIK3CA")
vaf_or <- or_maf@data$t_alt_count/or_maf@data$t_ref_count

DNAH8_maf <- subsetMaf(all_maf, genes = "DNAH8")
txnip_maf<- subsetMaf(all_maf, genes = "TXNIP")
################################################################################
gene_mafsubset = all_maf_lolipoplot[which(all_maf_lolipoplot$Hugo_Symbol == "DNAH8"),]
gene_protein_pos = gene_mafsubset$Protein_position
total_len = c()
pos = c()
for(i in gene_protein_pos){
  str = strsplit(i, "/", fixed = TRUE)
  str2 = strsplit(str[[1]][1], "-", fixed = TRUE)
  pos = append(pos,str2[[1]][1])
  total_len = append(total_len, str[[1]][2])
}
gene_mafsubset$AA_mut = pos
gene_mafsubset$Total_protein = total_len
mut_occurance = ddply(gene_mafsubset,.(AA_mut,Amino_acids,Variant_Classification, Total_protein), nrow)


gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
gff = readRDS(file = gff)
data.table::setDT(x = gff)
prot = gff[HGNC %in% "DNAH8"]
txs = unique(prot$refseq.ID)
if(length(txs) > 1){
  prot = prot[which(prot$aa.length == max(prot$aa.length)),]
}

prot = prot[,domain_lenght := End - Start][order(domain_lenght, decreasing = TRUE)][,domain_lenght := NULL]
domains = unique(prot[,Label])
domain_cols = brewer.pal(n = 8, name = "Pastel1")

if(length(domains) > length(domain_cols)){
  domain_cols = sample(colours(), size = length(domains), replace = FALSE)
}
domain_cols = domain_cols[1:length(domains)]
names(domain_cols) = domains
prot[, domainCol := domain_cols[prot[, Label]]]


################################################################################
###                     Venndigrams                                          ###
################################################################################

# venn for most mutated genes
venn_most_mut <- plot_custom_venn(most_mutated_genes, "most mutated genes per region")
# venn dataframe
venn_data_mostmut <- setNames(venn_most_mut@region$item,venn_most_mut@region$name)

nec_maf <- subsetMaf(regionmaf_list$NEC, genes = most_mut_genes$Var1)
t1_maf <- subsetMaf(regionmaf_list$T1, genes = most_mut_genes$Var1)
inf_maf <- subsetMaf(regionmaf_list$INF, genes = most_mut_genes$Var1)
ben_maf <- subsetMaf(regionmaf_list$BEN, genes = most_mut_genes$Var1)
# vennof all variants-------most mut
venn_most_mut_all <- plot_custom_venn( list(NEC = nec_maf@data$Hugo_Symbol,
                                            T1 = t1_maf@data$Hugo_Symbol,
                                            INF = inf_maf@data$Hugo_Symbol,
                                            BEN =ben_maf@data$Hugo_Symbol),
                                       "")

venn_data_allmut <- setNames(venn_all_variants@region$item,venn_all_variants@region$name)

(c(venn_data_mostmut$BEN %in% venn_data_allmut$BEN))

venn_data_mostmut$BEN


somaticInteractions(maf = all_maf, top = 25, pvalue = c(0.05, 0.1))

vaf <- compute_vaf(all_maf)
vaf2 <- vaf[vaf$Hugo_Symbol %in% c("TP53", "EGFR", "OR8U1", "CHIT1", "LNP1", "PTEN", "NF1"),]


genes_vaf <- subsetMaf(
  maf = all_maf, query = "t_depth > 0",
  fields = c("t_depth", "t_alt_count"), mafObj = FALSE
)
genes_vaf$VAF <- as.numeric(genes_vaf$t_alt_count) / as.numeric(genes_vaf$t_depth)
genes_vaf <- genes_vaf[, mean(VAF), Hugo_Symbol]
colnames(genes_vaf)[2] <- "VAF"
genes_vaf <- genes_vaf[which(genes_vaf$VAF > 0.02)]



################################################################################
#num of diff genes that are mutated
test <- unique(regionmaf_list$BEN@data$Hugo_Symbol)
# num of variants 
test2 <- sum(regionmaf_list$NEC@variant.type.summary$DEL)

sum(regionmaf_list$INF@variant.type.summary$SNP)
write.table(as.data.frame(table(regionmaf_list$NEC@gene.summary$AlteredSamples)),
            paste0( output_dir, dir ,"_sumary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
################################################################################


test_tcgacomp <- tcgaCompare(regionmaf_list$NEC, capture_size = kit_size/1000000)

regionmaf_list$NEC
regionmaf_list$T1
regionmaf_list$INF
regionmaf_list$BEN
set_exome <- list(
  NEC =regionmaf_list$NEC@data$Hugo_Symbol,
  T1= regionmaf_list$T1@data$Hugo_Symbol,
)




################################################################################
# #dgi = drugInteractions(maf = maf, fontSize = 0.75)
# novel_snvs[,1]
# dim(novel_snvs)
# if (dim(novel_snvs)[1] == 0) {
# print("hello")
# }
# #nix bei raus gekommen :(
# test <- clinicalEnrichment(all_maf, clinicalFeature = "Tumor_region")
# plotEnrichmentResults(enrich_res = test, pVal = 0.3)
# 
# 
# novel_snvs <- NEC_maf@data[,c("Hugo_Symbol", "dbSNP_RS")] %>%
#   filter(grepl("novel", dbSNP_RS) | dbSNP_RS == "")
# 
# make_oncoplot(output_dir, region_dir, maf, anno_color, "_novel_snvs", c(unique(novel_snvs$Hugo_Symbol)),sample_order = NULL)
# 
# test <- maf@gene.summary[maf@gene.summary$MutatedSamples >= 2,c("Hugo_Symbol", "MutatedSamples")] %>%
#   filter(Hugo_Symbol %in% novel_snvs$Hugo_Symbol)
# # -- gene selection befor plotting----------------------------------------------
# NEC_maf <- regionmaf_list$NEC
# NEC_maf@data[,c("Hugo_Symbol", "CLIN_SIG")]
# 
# oncoplot(
#   maf = nec_maf,
#   genes = unique(genes$Hugo_Symbol),
#   draw_titv = FALSE,
#   clinicalFeatures = c("Patient_ID", "Tumor_region"),
#   leftBarData = compute_vaf(nec_maf), # see function script
#   annotationColor = anno_color, # see function script
#   colors = vc_cols,
#   sortByAnnotation = TRUE,
#   gene_mar = 6
# )# see function script
# 
# genes <- NEC_maf@data[,c("Hugo_Symbol", "CLIN_SIG", "PUBMED", "dbSNP_RS", "Tumor_Sample_Barcode", "PolyPhen")] %>%
#   filter(CLIN_SIG != "")
# 
# 
# 
# getGeneSummary(NEC_maf)
# 
# nec_files <- list.files(path = "NEC/maf_converted/", full.names = TRUE, pattern = ".maf")
# nec_maf <- maf <- merge_mafs(nec_files, clinicalData = metadata)
# nec_anno <- nec_maf@data 
# nec_gene <- nec_maf@gene.summary
# #subset maf 
# library("vcfR")
# vcf <- read.vcfR("/Users/cschwitalla/Documents/WES_data/sarek3_annotation/manta/NEC/filtered/Manta_QATLV067AW_vs_QATLV069AE.diploidSV_VEP.ann_VEP.ann_filtered.vcf", verbose = FALSE)
# test <- vcf@fix
# 
# 
# #MY VARIANT TEST----------------------------------------------------------------
# 
# 
# 
# library("readxl")
# library("myvariant")
# stuf <- read_xlsx("/Users/cschwitalla/Downloads/Binders Medulloblastoma.xlsx")
# test_vcf <- readVcf("/Users/cschwitalla/Documents/exome_data_indexed/INF/raw_data/Strelka/somatic/filtered/concat_snvs_filtered.vcf")
# hgvs <- formatHgvs(test_vcf, variant_type = "snp")
# 
# head(hgvs)
# 
# 
# test_clinvar <- getVariants(hgvs, fields = "clinvar.type")
# test_cadd <- getVariants(hgvs, fields = "cadd.consequence")
# 
# # include lollipopplot maybe
# # circouise plots
# #circos.initializeWithIdeogram(species = "hg19")
# 
# # circos plits ---------------
# bed_file <- read.table("/Users/cschwitalla/Documents/exome_data_indexed/INF/raw_data/Strelka/somatic/filtered/all_unsorted.bed", sep = "\t", header = FALSE)
# bed_file <- bed_file[c(1:3, 6:7)]
# names(bed_file) <- c("chr", "start", "end", "ref", "alt")
# circos.initializeWithIdeogram(bed_file)
# circos.genomicTrack(test,)
# test_split <- split(bed_file,bed_file$chr)
# 
# 
# chr1 <- test_split$chr1
# test_chr1 <- table(chr1)
# 
# test_freq <- apply(chr1[c("start", "end", "ref", "alt")], 2, table)
# 
# 
# 
# circos.initializeWithIdeogram(plotType = NULL)
# circos.track(ylim = c(0,1), panel.fun = function(x,y){
#   chr = CELL_META$sector.index
#   xlim = CELL_META$xlim
#   ylim = CELL_META$ylim
#   circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
#   circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
#               facing = "inside", niceFacing = TRUE)
# }, track.height = 0.15, bg.border = NA)  
# 
# 
# 
# 
# circos.genomicTrackPlotRegion(bed_file, panel.fun = function(region, value , ...) {
#   if(CELL_META$sector.index == "chr1") {
#     print(head(region, n = 2))
#     print(head(value, n = 2))
#   }
# })
#   
# 
# 
# 
# 
