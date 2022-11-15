################################################################################
###                                Purpose                                   ###
################################################################################
## Title: Functions for WES Analysis
## Author: Carolin Schwitalla
##
## Description: This script contains all functions that are used for the WES
##              analysis for my master thesis
################################################################################
###                             Hard coded Variables                         ###
################################################################################
# colors for oncoplot itself
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
# plot summary 
plot_summary <- function(output_dir,dir,maf){
  pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)],"_maf_summary.pdf"),
      width = 8,
      height = 6)
  plotmafSummary(maf = maf)
  dev.off()
}

#drug interactions
drug_interactions <- function(output_dir,dir,maf) {
  pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)],"_drug_interactions.pdf"),
      width = 7,
      height = 8)
  drugInteractions(maf, fontSize = 0.75)
  dev.off()
}


#oncogenic pathways
onco_pathways <- function(maf, output_dir, dir) {
  pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)],"_onco_pathway.pdf"),
      width = 7,
      height = 8)
  OncogenicPathways(maf = maf)
  dev.off()
}


# tcga compare
tcga_comapre <- function(maf, dir, kit_size, output_dir, name) {
  if (is.null(dir)){
    print("hello ")
    pdf(paste0(output_dir, name,"_mutload.pdf"),
        width = 10,
        height = 5)
    tcgaCompare(maf = maf, cohortName = name, logscale = TRUE, capture_size = kit_size/1000000)
    dev.off()
  }else {
    pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)],"_mutload.pdf"),
        width = 10,
        height = 5)
    tcgaCompare(maf = maf, cohortName = names(region_dir)[match(dir, region_dir)], logscale = TRUE, capture_size = kit_size/1000000)
    dev.off()
  }

}




# #make oncoplot pdf for novel stuff
# make_oncoplot_novel_snvs <- function(output_dir, region_dir, maf, anno_color, title, genes, sample_order) {
#   
#   pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)], title,"_oncoplot.pdf"),
#       width = 10,
#       height = 8)
#   # make oncoplot for current maf that will exported as pdf in the output dir
#   oncoplot(
#     maf = maf,
#     genes = genes,
#     draw_titv = FALSE,
#     clinicalFeatures = c("Tumor_region", "Patient_ID"),
#     leftBarData = compute_vaf(maf), # see function script
#     annotationColor = anno_color, # see function script
#     colors = vc_cols,
#     #fontSize = 0.4,
#     sampleOrder = sample_order,
#     sortByAnnotation = TRUE,
#     gene_mar = 6
#   )# see function script
#   dev.off()
#   
#   
#   
# }


#make oncoplot pdf
make_oncoplot <- function(output_dir, region_dir, maf, anno_color, title, genes, sample_order, top, minMut, rightBarData) {
    vaf <- compute_vaf(maf)
    # filter out genes that have vaf < 0.02
    if(!is.null(genes)){
      genes <- genes[genes %in% vaf$Hugo_Symbol]
    }
    pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)], title,"_oncoplot.pdf"),
        width = 10,
        height = 8)
    # make oncoplot for current maf that will exported as pdf in the output dir
    oncoplot(
      maf = subsetMaf(maf,genes = vaf$Hugo_Symbol),
      genes = genes,
      draw_titv = FALSE,
      clinicalFeatures = c("Tumor_region", "Patient_ID"),
      leftBarData = vaf, # see function script
      annotationColor = anno_color, # see function script
      colors = vc_cols,
      sampleOrder = sample_order,
      sortByAnnotation = TRUE,
      gene_mar = 6,
      
      #top = top,
      fontSize = 0.6
      #rightBarData = rightBarData
    )# see function script
    dev.off()
    return(vaf)

    
  
}
compute_vaf2 <- function(maf_obj) {
  vaf <- maf_obj@data[c("Hugo_Symbol", "t_alt_count", "t_depth")]
  vaf_genes$Hugo_Symbol <- unique(vaf$Hugo_Symbol)
  vaf$vaf <- vaf$t_alt_count/vaf_t_depth
  vaf
}


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



## DESCRIPTION: mapping meta data to colors for annotation color in oncoplots
## PARAMETERS:
##             - dregion_column: matadata df column were tumor region / conditions
##                               are listed
##             - patient_column: metadata df column were patien ids are listed
## OUTPUT: returns a list of lists
##
create_annotation_color <- function(patient_column, region_column) {
  patient_scale <- colorRampPalette(c("#543005", "#f5f5f5", "#003c30"))
  region_color <- c("#4B6C22", "#74add1", "#9e0142", "#fdae61")
  #meth_color <- c("#4B6C22", "#74add1")
  #get the number of unique patient ids to extract colors from color scale
  sample_num <- length(unique(patient_column))
  # get a specific color palette with num of patients
  patient_color <- patient_scale(sample_num)
  # set names to asign for each level the right color
  names(patient_color) <- unique(patient_column)
  #names(meth_color) <- sort(unique(meth_column))
  names(region_color) <- sort(unique(region_column))
  annotation_color <- list(Tumor_region = region_color,
                           Patient_ID = patient_color)
  return(annotation_color)
}


## DESCRIPTION: lollipop plots to visualize amino acid changes in protein sequence
##              from maf file. The code from the function lollipopPlot from the 
##              package maftools was used as a basis and modified to fit maffiles
##              without columns like HGVSp_Short, protein_Change, AAChange, and 
##              using ggplot for plotting
##              
## PARAMETERS:
##             - genename: name of gene which will be plotted [string]
##             - maffile_subset: maffile@data 
## OUTPUT: returns a list of lists
##
Lollipop_plot <- function(genename, maffile_subset){
  # adapted from package maftools --> lollipopPlot()
  gene_mafsubset = maffile_subset[which(maffile_subset$Hugo_Symbol == genename),]
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
  
  geneID = genename
  #Protein domain source.
  gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
  gff = readRDS(file = gff)
  data.table::setDT(x = gff)
  prot = gff[HGNC %in% geneID]
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
  
  # Plot
  ggplot()+#mut_occurance, mapping=aes(x=as.numeric(mut_occurance$AA_mut), y=mut_occurance$V1, label = paste(mut_occurance$Amino_acids, mut_occurance$AA_mut, sep = ": "))) +
    annotate("rect", xmin = 1, xmax = as.numeric(total_len), ymin = -0.1, ymax = 0.1, alpha = 0.4, fill = "grey") +
    geom_segment( mut_occurance,mapping=aes(x=as.numeric(AA_mut), xend=as.numeric(AA_mut), y=0, yend=V1), color ="grey") +
    geom_point(mut_occurance,mapping=aes(colour = factor(Variant_Classification), x=as.numeric(AA_mut),y=V1), size=4) + 
    #geom_text(nudge_x = 0, nudge_y = 0.4) +
    geom_text_repel(mut_occurance,mapping=aes(x=as.numeric(AA_mut),y=V1+0.15,label = paste(Amino_acids, AA_mut, sep = ": ")))+
    geom_rect(prot,mapping=aes(NULL,NULL,xmin = prot[,Start], xmax= prot[,End], fill = prot[,Label]), ymin = -0.2, ymax = 0.2, alpha = 0.7 ) +
    #geom_segment( mut_occurance,mapping=aes(x=as.numeric(AA_mut), xend=as.numeric(AA_mut), y=0, yend=V1), color ="grey") +
    scale_x_continuous(limits = c(1,as.numeric(total_len)))+
    #annotate("rect", xmin = 1, xmax = as.numeric(total_len), ymin = -0.25, ymax = 0.25, alpha = 0.5, fill = "grey") +
    #annotate("rect", xmin = prot[,Start], xmax = prot[,End], ymin = -0.4, ymax = 0.4, fill= prot[,domainCol], alpha = 0.75, color = "lightsteelblue4")+
    #annotate("text", x = ((prot[,End]-prot[,Start])/2) + prot[,Start], y = 0, label = prot[,Label])+
    labs(x = "protein length [amino acids]", y = "number of mutated samples")+
    guides(color = guide_legend(title = "type of mutation")) +
    guides(fill = guide_legend(title = "protein domain")) +
    theme_classic() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  
}









