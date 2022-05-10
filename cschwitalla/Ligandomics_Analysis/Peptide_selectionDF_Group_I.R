setwd("/Users/cschwitalla/Documents/Immunopeptidomics/")

library("stringr")
library("tidyr")
library("dplyr")
library("tibble")
library("tuple")
library("org.Hs.eg.db") 
library("readxl")
library("EnsDb.Hsapiens.v86")
library("ensembldb")
library("AnnotationDbi")
library("readxl")



#read in file from inhouse databank-----------------------
inhouseDB_search = read_xlsx("/Users/cschwitalla/Downloads/number_of_observations.xlsx", col_names = TRUE)

#HUMAN PROTEIN ATLAS-------------
normal_tissue_proteom = read.csv("/Users/cschwitalla/Downloads/normal_tissue.tsv", header = TRUE, sep = "\t")
cancer_tissue_proteom = read.csv("/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/pathology.tsv", header = TRUE, sep = "\t")


#read in data for gioele---------------------------
NEC_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_frequend_peptides.tsv")
T1_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_frequend_peptides.tsv")
INF_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_frequend_peptides.tsv")
#shared_without_ben_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/shared_I_high_f.tsv")
Tumor_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_I_frequend_peptides.tsv")


NEC_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_frequend_peptides.tsv")
T1_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_frequend_peptides.tsv")
INF_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_frequend_peptides.tsv")
#shared_without_ben_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/shared_II_high_f.tsv")
Tumor_II_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_II_frequend_peptides.tsv")



#tumor region peptides 
NEC_I_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_seq.tsv")
T1_I_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_seq.tsv")
INF_I_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_seq.tsv")
BEN_I_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/BEN_I_seq.tsv")


NEC_II_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_seq.tsv")
T1_II_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_seq.tsv")
INF_II_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_seq.tsv")
BEN_II_seq = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/BEN_II_seq.tsv")






#MAP UNIPROT ACCESSIONS TO GENE NAMES----------------------------------
#genename to protein acc
edb <- EnsDb.Hsapiens.v86
hasProteinData(edb)
listTables(edb)
accesionList = c(Tumor_acc_pep$Accessions, Tumor_II_acc_pep$Accessions)

GeneNames_list = sapply(accesionList, function(gen){
  genname = AnnotationDbi::select(edb, keys = gen, keytype = "UNIPROTID", columns = "GENENAME")[1,][2]
  
})
Acc = c()
Gene = c()
for(i in 1:length(GeneNames_list)){
  
  acc = names(GeneNames_list[i])
  gene = GeneNames_list[[i]]
  Acc = append(Acc,acc)
  if(is.na(gene)){
    Gene = append(Gene, "NO")
  }
  else{
    Gene = append(Gene, gene)
  }
  
}


accessions = c()
for(i in 1: length(Acc)){
  str = strsplit(Acc[i], ".", fixed = TRUE)[[1]][1]
  accessions = append(accessions, str)
}

mapping_df = as.data.frame(accessions)
mapping_df$GeneNames = Gene
mapping_df = unique(mapping_df)
mapping_df = as.data.frame(mapping_df)
names(mapping_df)[1] = "Accessions"
names(mapping_df)[2] = "GeneName"

#include to df 
Tumor_acc_pep = merge(Tumor_acc_pep,mapping_df, by = "Accessions", all.x = TRUE)
NEC_excl_acc_pep = merge(NEC_excl_acc_pep, mapping_df, by = "Accessions", all.x = TRUE)
T1_excl_acc_pep = merge(T1_excl_acc_pep, mapping_df, by = "Accessions", all.x = TRUE)
INF_excl_acc_pep = merge(INF_excl_acc_pep, mapping_df, by = "Accessions", all.x = TRUE)

Tumor_II_acc_pep = merge(Tumor_II_acc_pep,mapping_df, by = "Accessions", all.x = TRUE)
NEC_excl_acc_pep_II = merge(NEC_excl_acc_pep_II, mapping_df, by = "Accessions", all.x = TRUE)
T1_excl_acc_pep_II = merge(T1_excl_acc_pep_II, mapping_df, by = "Accessions", all.x = TRUE)
INF_excl_acc_pep_II = merge(INF_excl_acc_pep_II, mapping_df, by = "Accessions", all.x = TRUE)


#PEPTIDE SELECTION FOR GIOELE -----------------------------------------------------------------------
#Tumor-exclusive peptides, presented on at least 2 patients --> xy-excl_acc_pep dfs enthalten das DONE



#Sorted based on ascending order of protein frequency in healthy tissue database --> in wie viele healthy tissues kommt es vor ? ====================
#https://www.proteinatlas.org/about/download-->normal tissue data 
#Expression profiles for proteins in human tissues based on immunohistochemisty using tissue micro arrays

# filter for level--> all not detected will be removed
geneNames_list = c(Tumor_acc_pep$GeneName, Tumor_II_acc_pep$GeneName)
normal_tissue_proteom_filter = normal_tissue_proteom[-which(normal_tissue_proteom$Level == "Not detected"),] #begründung
normal_tissue_proteom_filter = normal_tissue_proteom_filter[-which(normal_tissue_proteom_filter$Reliability == "Uncertain"),] #begründung
normal_tissue_proteom_filter = normal_tissue_proteom_filter %>% dplyr::filter(normal_tissue_proteom_filter$Gene.name %in% geneNames_list)

# abundace in tissue 
uniq_genes = unique(normal_tissue_proteom_filter$Gene.name)
normal_tissue_proteinf_df = data.frame(unique(normal_tissue_proteom_filter$Gene.name))
normal_tissue_proteom_freq = normal_tissue_proteom_filter %>% group_by(Gene.name) %>% summarise(Tissue = unique(Tissue))

normal_tissue_proteinf_df$Frequnecy = sapply(uniq_genes, function(gen){
  mean(ifelse(is.na(matchAll(gen, normal_tissue_proteom_freq$Gene.name)), 0,
              (length(matchAll(gen, normal_tissue_proteom_freq$Gene.name))/length(unique(normal_tissue_proteom_freq$Tissue)))*100))
})

length(unique(normal_tissue_proteom_freq$Tissue))
normal_tissue_df_protein_abundance = normal_tissue_proteinf_df
names(normal_tissue_df_protein_abundance)[1] = "GeneName"
names(normal_tissue_df_protein_abundance)[2] = "NormalTissueFrequency"

#combine normal tissue frequency 

Tumor_acc_pep = merge(Tumor_acc_pep, normal_tissue_df_protein_abundance, by ="GeneName", all.x = TRUE)
NEC_excl_acc_pep = merge(NEC_excl_acc_pep, normal_tissue_df_protein_abundance, by = "GeneName", all.x = TRUE)
T1_excl_acc_pep = merge(T1_excl_acc_pep, normal_tissue_df_protein_abundance, by = "GeneName", all.x = TRUE)
INF_excl_acc_pep = merge(INF_excl_acc_pep, normal_tissue_df_protein_abundance, by = "GeneName", all.x = TRUE)

Tumor_II_acc_pep = merge(Tumor_II_acc_pep, normal_tissue_df_protein_abundance, by ="GeneName", all.x = TRUE)
NEC_excl_acc_pep_II = merge(NEC_excl_acc_pep_II, normal_tissue_df_protein_abundance, by = "GeneName", all.x = TRUE)
T1_excl_acc_pep_II = merge(T1_excl_acc_pep_II, normal_tissue_df_protein_abundance, by = "GeneName", all.x = TRUE)
INF_excl_acc_pep_II = merge(INF_excl_acc_pep_II, normal_tissue_df_protein_abundance, by = "GeneName", all.x = TRUE)

#Descending order of peptide frequency in glioblastoma tissue samples ===================================000

inhouseDB_search$total = inhouseDB_search$benign + inhouseDB_search$glioblastoma + inhouseDB_search$`other cancer entities`
inhouseDB_search$GlioblastomaFreq = inhouseDB_search$glioblastoma / inhouseDB_search$total *100
inhouseDB_search$BenignFreq = inhouseDB_search$benign / inhouseDB_search$total *100
inhouseDB_search$OtherCancerFreq = inhouseDB_search$`other cancer entities`/ inhouseDB_search$total *100
Immuno_DB_pep_f = inhouseDB_search[c(2,7,8,9)]

Tumor_acc_pep = merge(Tumor_acc_pep, Immuno_DB_pep_f, by ="Sequence", all.x = TRUE)
NEC_excl_acc_pep = merge(NEC_excl_acc_pep, Immuno_DB_pep_f, by = "Sequence", all.x = TRUE)
T1_excl_acc_pep = merge(T1_excl_acc_pep, Immuno_DB_pep_f, by = "Sequence", all.x = TRUE)
INF_excl_acc_pep = merge(INF_excl_acc_pep, Immuno_DB_pep_f, by = "Sequence", all.x = TRUE)

Tumor_II_acc_pep = merge(Tumor_II_acc_pep, Immuno_DB_pep_f, by ="Sequence", all.x = TRUE)
NEC_excl_acc_pep_II = merge(NEC_excl_acc_pep_II, Immuno_DB_pep_f, by = "Sequence", all.x = TRUE)
T1_excl_acc_pep_II = merge(T1_excl_acc_pep_II, Immuno_DB_pep_f, by = "Sequence", all.x = TRUE)
INF_excl_acc_pep_II = merge(INF_excl_acc_pep_II, Immuno_DB_pep_f, by = "Sequence", all.x = TRUE)

#Descending order of protein source frequency in glioblastoma tissue samples ==================================================
#hier pep frequencies from own data set
#protein frequencies werdenoben schon eingeladen deshlab brauch ich das nicht mehr 


#protein frequency in non-glioblastoma malignant samples in descending order =============
cancer_tissue_proteom$totalPatients = rowSums(cancer_tissue_proteom[4:7])
cancer_tissue_proteom$PatientNumberDetected = rowSums(cancer_tissue_proteom[4:6])
cancer_tissue_proteom$detectionPercentage = cancer_tissue_proteom$PatientNumberDetected / cancer_tissue_proteom$totalPatients *100

#filtern nach proteinen die in meiner selection drin sind 
cancer_tissue_proteom_filtered = cancer_tissue_proteom %>% dplyr::filter(cancer_tissue_proteom$Gene.name %in% geneNames_list)
cancer_tissue_proteom_filtered = cancer_tissue_proteom_filtered[which(cancer_tissue_proteom_filtered$Cancer != "glioma"), ]

cancer_tissue_proteom_filtered_small = cancer_tissue_proteom_filtered[c(2,3,14)]
cancer_tissue_proteom_filtered_small_means = cancer_tissue_proteom_filtered_small %>% group_by(Gene.name) %>% summarise_at(vars(detectionPercentage), list(Mean = mean))
names(cancer_tissue_proteom_filtered_small_means)[1] = "GeneName"
names(cancer_tissue_proteom_filtered_small_means)[2] = "Cancer_tissue_protein_frequency"




Tumor_I_peptide_selection = merge(Tumor_acc_pep, cancer_tissue_proteom_filtered_small_means, by ="GeneName", all.x = TRUE)
NEC_excl_I_peptide_selection = merge(NEC_excl_acc_pep, cancer_tissue_proteom_filtered_small_means, by = "GeneName", all.x = TRUE)
T1_excl_I_peptide_selection = merge(T1_excl_acc_pep, cancer_tissue_proteom_filtered_small_means, by = "GeneName", all.x = TRUE)
INF_excl_I_peptide_selection = merge(INF_excl_acc_pep, cancer_tissue_proteom_filtered_small_means, by = "GeneName", all.x = TRUE)

Tumor_II_peptide_selection = merge(Tumor_II_acc_pep, cancer_tissue_proteom_filtered_small_means, by ="GeneName", all.x = TRUE)
NEC_excl_II_peptide_selection = merge(NEC_excl_acc_pep_II, cancer_tissue_proteom_filtered_small_means, by = "GeneName", all.x = TRUE)
T1_excl_II_peptide_selection = merge(T1_excl_acc_pep_II, cancer_tissue_proteom_filtered_small_means, by = "GeneName", all.x = TRUE)
INF_excl_II_peptide_selection = merge(INF_excl_acc_pep_II,cancer_tissue_proteom_filtered_small_means, by = "GeneName", all.x = TRUE)






# final df
names(Tumor_I_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")
names(NEC_excl_I_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")
names(T1_excl_I_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")
names(INF_excl_I_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")


names(Tumor_II_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")
names(NEC_excl_II_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")
names(T1_excl_II_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")
names(INF_excl_II_peptide_selection) = c("GeneName", "Sequence", "Accessions", "Peptide_frequency_own_Dataset","Protein_frequency_own_Dataset","HealthyTissue_ProteinFrequency", "Peptide_Freq_Glioblastoma_external",
                                     "Peptide_Freq_Benign_external","Peptide_Freq_otherCancer_external", "CancerTissue_protein_frequency")

Tumor_I_peptide_selection = Tumor_I_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]
NEC_excl_I_peptide_selection = NEC_excl_I_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]
T1_excl_I_peptide_selection = T1_excl_I_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]
INF_excl_I_peptide_selection = INF_excl_I_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]

Tumor_II_peptide_selection = Tumor_II_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]
NEC_excl_II_peptide_selection = NEC_excl_II_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]
T1_excl_II_peptide_selection = T1_excl_II_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]
INF_excl_II_peptide_selection = INF_excl_II_peptide_selection[,c(1,2,3,4,5,7,8,9,6,10)]



#MERGE THE BINDERS PREDICTION !!! U DUMB ASS -----------


NEC_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_high_f.tsv")
T1_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_high_f.tsv")
INF_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_high_f.tsv")
Tumor_I_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_I_high_f.tsv")


NEC_II_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_high_f.tsv")
T1_II_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_high_f.tsv")
INF_II_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_high_f.tsv")
Tumor_II_binder_pred = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_II_high_f.tsv")


#merge

Tumor_I_peptide_selection = merge(Tumor_I_peptide_selection, Tumor_I_binder_pred[c(1,4,5,6)], by ="Sequence", all.x = TRUE)
NEC_excl_I_peptide_selection = merge(NEC_excl_acc_pep, NEC_binder_pred[c(1,4,5,6)], by = "Sequence", all.x = TRUE)
T1_excl_I_peptide_selection = merge(T1_excl_acc_pep, T1_binder_pred[c(1,4,5,6)], by = "Sequence", all.x = TRUE)
INF_excl_I_peptide_selection = merge(INF_excl_acc_pep, INF_binder_pred[c(1,4,5,6)], by = "Sequence", all.x = TRUE)

Tumor_II_peptide_selection = merge(Tumor_II_peptide_selection, Tumor_II_binder_pred[c(1,4,5,6)], by ="Sequence", all.x = TRUE)
NEC_excl_II_peptide_selection = merge(NEC_excl_acc_pep_II, NEC_II_binder_pred[c(1,4,5,6)], by = "Sequence", all.x = TRUE)
T1_excl_II_peptide_selection = merge(T1_excl_acc_pep_II, T1_II_binder_pred[c(1,4,5,6)], by = "Sequence", all.x = TRUE)
INF_excl_II_peptide_selection = merge(INF_excl_acc_pep_II,INF_II_binder_pred[c(1,4,5,6)], by = "Sequence", all.x = TRUE)

#merge tumors 

Tumor_I_peptide_selection_test= unique(list(Tumor_I_peptide_selection,NEC_I_seq,T1_I_seq,INF_I_seq,BEN_I_seq) %>% reduce(left_join, by = "Sequence"))
Tumor_II_peptide_selection_test = unique(list(Tumor_II_peptide_selection,NEC_II_seq,T1_II_seq,INF_II_seq,BEN_II_seq) %>% reduce(left_join, by = "Sequence"))

names(Tumor_I_peptide_selection_test)[14:17] = c("NEC_region", "T1_region", "INF_region", "BEN_region")
names(Tumor_II_peptide_selection_test)[14:17] = c("NEC_region", "T1_region", "INF_region", "BEN_region")
#write df 


write.table(Tumor_I_peptide_selection_test,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_I_peptide_selection.tsv")
write.table(NEC_excl_I_peptide_selection,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_excl_I_peptide_selection.tsv")
write.table(T1_excl_I_peptide_selection,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_excl_I_peptide_selection.tsv")
write.table(INF_excl_I_peptide_selection,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_excl_I_peptide_selection.tsv")
write.table(Tumor_II_peptide_selection_test,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_II_peptide_selection.tsv")
write.table(NEC_excl_II_peptide_selection,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_excl_II_peptide_selection.tsv")
write.table(T1_excl_II_peptide_selection,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_excl_II_peptide_selection.tsv")
write.table(INF_excl_II_peptide_selection,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_excl_II_peptide_selection.tsv")











#PEPTIDE SELECTION GROUP II TRANSCRIPTOMIC APPROACH----------------------------------------------











