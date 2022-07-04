################################################################################
###                                Purpose                                   ###
################################################################################

################################################################################
###                               Manual steps                               ###
################################################################################
# clear R enviroment
rm(list = ls())
## TODO: set workig dir to file source dir 
# set directory fpr the input files
input_dir <- ("/Users/cschwitalla/Documents/Immunopeptidomics/")
# set directory for the output 
output_dir <-  ("/Users/cschwitalla/Documents/Ligandomics_analysis/")
# import functions for the analysis 
source("/Users/cschwitalla/git/students/cschwitalla/Ligandomics_Analysis/functions_ligandomics.R")
################################################################################
###                             Load libraries                               ###
################################################################################
required_Libs <- c("tidyr","readxl", "ggVennDiagram", "dplyr", "stringr",
                   "tibble", "ggplot2", "org.Hs.eg.db")

suppressMessages(invisible(lapply(required_Libs, library, character.only = T)))


################################################################################
###                            Load Data                                     ###
################################################################################
# load HLA-typing dataframe form Marcel Immunology------------------------------
GB_HLA_types <- read_xlsx(paste0(input_dir, "HLA-Typisierung_GBM.xlsx"), col_names = TRUE)

# get list of unique HLA types
uniqe_HLA_types <- unique(c(as.matrix(GB_HLA_types[2:16, 2:7])))
# rewrite HLA types that NetMHCpan can use them and save them as vector
new_HLA_types <- rewrite_HLA_types(uniqe_HLA_types, as_string = FALSE)

# Benign data Immunology -------------------------------------------------------
# more specific
# less hits
benign_pep_I <- read.csv(paste0(input_dir, "newBenignmorespecific/Benign_class1.csv"),
  header = FALSE,
  sep = ","
)[, 1:2]
benign_pep_II <- read.csv(paste0(input_dir, "newBenignmorespecific/Benign_class2.csv"),
  header = FALSE,
  sep = ","
)[, 1:2]


#HLA-Ligand Atlas data ---------------------------------------------------------
# read in dataframes from HLA ligand antlas
HLA_ligand_atlas_pep <- read.csv(paste0(input_dir, "hla_2020.12/HLA_aggregated.tsv"),
  header = TRUE,
  sep = "\t"
)
HLA_ligand_atlas_acc <- read.csv(paste0(input_dir, "hla_2020.12/HLA_protein_map.tsv"),
  header = TRUE,
  sep = "\t"
)
# aggregate HLA ligand atlas dataframe bevor merging 
HLA_ligand_atlas_acc <- HLA_ligand_atlas_acc %>%
  group_by(peptide_sequence_id ) %>% summarise(acc = toString(uniprot_id))

# merge dataframes by peptide sequence id to add protein accession numbers to the dataframe
HLA_atlas_data <- merge(HLA_ligand_atlas_pep, HLA_ligand_atlas_acc, by = "peptide_sequence_id")
# split the ligand atlas according to the hla class ( 1 or 2 )
HLA_atlas_data <- split(HLA_atlas_data, HLA_atlas_data$hla_class)

# TODO maybe--------------------------------------------------------------------
# #get mapping frequencie for every peptide id 
# HLA_Ligand_mapperfreq = as.data.frame(table(HLA_Ligand_acc_df$peptide_sequence_id))
# # get peptide ids with freq 1
# Uniqe_mappers = HLA_Ligand_mapperfreq[which(HLA_Ligand_mapperfreq$Freq == 1),]
# # get all acc that are uniqe mappers
# HLA_Ligand_atlas_unique_acc = HLA_Ligand_acc_df %>% dplyr::filter(HLA_Ligand_acc_df$peptide_sequence_id  %in% Uniqe_mappers$Var1)
# #split in class I nd class II 
# HLA_Atlas_classI = HLA_Ligand_atlas_df[which((HLA_Ligand_atlas_df$hla_class == "HLA-I")|(HLA_Ligand_atlas_df$hla_class == "HLA-I+II")),]
# HLA_Atlas_classII = HLA_Ligand_atlas_df[which((HLA_Ligand_atlas_df$hla_class == "HLA-II")|(HLA_Ligand_atlas_df$hla_class == "HLA-I+II")),]
# # split protein acc in classI and classII 
# HLA_Atlas_uniqueacc_classI = HLA_Atlas_classI %>% dplyr::filter(HLA_Atlas_classI$peptide_sequence_id %in% HLA_Ligand_atlas_unique_acc$peptide_sequence_id)
# HLA_Atlas_uniqueacc_classII  = HLA_Atlas_classII %>% dplyr::filter(HLA_Atlas_classII$peptide_sequence_id %in%HLA_Ligand_atlas_unique_acc$peptide_sequence_id)
# 
# HLA_Atlas_I_uniq_acc = HLA_Ligand_atlas_unique_acc %>% dplyr::filter(HLA_Ligand_atlas_unique_acc$peptide_sequence_id %in% HLA_Atlas_uniqueacc_classI$peptide_sequence_id)
# HLA_Atlas_II_uniq_acc = HLA_Ligand_atlas_unique_acc %>% dplyr::filter(HLA_Ligand_atlas_unique_acc$peptide_sequence_id %in% HLA_Atlas_uniqueacc_classII$peptide_sequence_id)
# 
# # add source protein acc to hla ligand atlas df 
# HLA_Atlas_I_uniqaccs = merge(HLA_Atlas_uniqueacc_classI, HLA_Atlas_I_uniq_acc, by= "peptide_sequence_id")
# HLA_Atlas_II_uniqaccs = merge(HLA_Atlas_uniqueacc_classII, HLA_Atlas_II_uniq_acc, by= "peptide_sequence_id")


# TUMOR REGION LIGANDOMICS DATA-------------------------------------------------
# get all filenames with path 
files_cI <- dir(paste0(input_dir, "ClassI"),
  recursive = TRUE,
  pattern = ".csv",
  full.names = TRUE
)
files_cII <- dir(paste0(input_dir, "ClassII"),
  recursive = TRUE,
  pattern = ".csv",
  full.names = TRUE
)

#read in data 
classI_df = create_data_frame(files_cI)
classII_df = create_data_frame(files_cII)

# DATA PREPERATION -------------------------------------------------------------
# combine benign data 

combi_benign_pep_I <- c(HLA_atlas_data$`HLA-I`$peptide_sequence,
                        HLA_atlas_data$`HLA-I+II`$peptide_sequence,
                        benign_pep_I$V1)
combi_benign_pep_II <- c(HLA_atlas_data$`HLA-II`$peptide_sequence,
                         HLA_atlas_data$`HLA-I+II`$peptide_sequence,
                         benign_pep_II$V1)
# filter my datasets so that the known benign peptides frombenign pep db are not 
# there anymore 
classI_df <- subset(classI_df, !(Sequence %in% combi_benign_pep_I))
classII_df <- subset(classII_df, !(Sequence %in% combi_benign_pep_II))

# dataset that hase only peptides predicted to originate only from one protein
# not multiple  
classI_df = classI_df[ -grep(";", classI_df$Accessions), ]
classII_df = classII_df[ -grep(";", classII_df$Accessions), ]

# create than dataframes that are filterd+uniqe origin that are region specific
region_specific_I <- split(classI_df, classI_df$Tumor_region)
region_specific_II <- split(classII_df,classII_df$Tumor_region)
   
# TODO: save dataframes to tsv

################################################################################
##                               Analysis                                     ##
################################################################################


# venn diagram -----------------------------------------------------------------
# extract venn diagrams region exclusive peptides than 
set_I <- setNames(vector("list", length = length(names(region_specific_I))),
                  c(names(region_specific_I)))

set_II <- setNames(vector("list", length = length(names(region_specific_II))),
                  c(names(region_specific_II)))

for (i in region_specific_I) {
  set_I[i$Tumor_region[1]] <- list(i$Sequence)
}
for (i in region_specific_II) {
  set_II[i$Tumor_region[1]] <- list(i$Sequence)
}

venn_data_I <- plot_custom_venn(set_I, "HLA class I Peptides")
venn_data_II <- plot_custom_venn(set_II, "HLA class II Peptides")
#test2@region$item



#length distribution CLASS I + II ================

#TO DO: 
# - color adjust
# - % statt counts

classI_df$length = apply(classI_df,2,nchar)[,4]
table(classI_df$length)
counts_classI = aggregate(classI_df$length, by=list(Tumor_region = classI_df$Tumor_region, len = classI_df$length ), FUN=table )
table(classII_df$length)
counts_classI = counts_classI[order(counts_classI$Tumor_region),]
sum(counts_classI$x)


classII_df$length = apply(classII_df,2,nchar)[,4]
counts_classII = aggregate(classII_df$length, by=list(Tumor_region = classII_df$Tumor_region, len = classII_df$length ), FUN=table )
counts_classII = counts_classII[order(counts_classII$Tumor_region),]
sum(counts_classII$x)



#plot as bar plot 
ggplot(counts_classI, aes(x = len, y = x, group = Tumor_region, fill = Tumor_region)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(4,"Set2")) +
  labs(x="Peptide length", y= "counts of HLA class I peptides") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

ggplot(counts_classII, aes(x = len, y = x, fill = Tumor_region)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(4,"Set2")) +
  labs(x="Peptide length", y= "counts of HLA class II peptides") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))


#plot 

ggplot(counts_classI ,aes(x = len, y = x, group = Tumor_region, color = Tumor_region)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(4,"Set2")) +
  labs(x="Peptide length", y= "counts of HLA class I peptides") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

ggplot(counts_classII ,aes(x = len, y = x, group = Tumor_region, color = Tumor_region)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(4, "Set2")) +
  labs(x="Peptide length", y= "counts of HLA class II peptides") +
  theme_classic() +
  theme( axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

display.brewer.pal(n = 4,name = "Paired")




#saturation analysis ===================

# y axis = num of unique peptides
# x axis = num of samples

sat_classI <- classI_df[match(unique(classI_df$Sequence), classI_df$Sequence), ] %>%
  group_by(Sample_num) %>%
  summarise(num_unique = length(unique(Sequence))) %>%
  arrange(desc(num_unique)) %>%
  mutate(num_unique = cumsum(num_unique)) %>%
  mutate(num_samples = seq(1, length(Sample_num)))

sat_classII <- classII_df[match(unique(classII_df$Sequence), classII_df$Sequence), ] %>%
  group_by(Sample_num) %>%
  summarise(num_unique = length(unique(Sequence))) %>%
  arrange(desc(num_unique)) %>%
  mutate(num_unique = cumsum(num_unique)) %>%
  mutate(num_samples = seq(1, length(Sample_num)))

# fit curve 

fit_param_classI = NLSstAsymptotic(sortedXyData(c(0, sat_classI$num_samples), c(0, sat_classI$num_unique)))
fit_sat_classI = SSasymp(sat_classI$num_samples, abs(fit_param_classI[1]) + fit_param_classI[2], fit_param_classI[1], fit_param_classI[3])

fit_param_classII = NLSstAsymptotic(sortedXyData(c(0, sat_classII$num_samples), c(0, sat_classII$num_unique)))
fit_sat_classII = SSasymp(sat_classII$num_samples, abs(fit_param_classII[1]) + fit_param_classII[2], fit_param_classII[1], fit_param_classII[3])

fit_param_classI

#plot 
ggplot(sat_classI,aes(x = sat_classI$num_samples, y = fit_sat_classI)) +
  geom_smooth(method = "loess", se = F, color = "black",formula = y ~ x) +
  geom_hline(yintercept = abs(fit_param_classI[1])+fit_param_classI[2], linetype = "dashed", alpha = 0.6) +
  labs(x = "Number of Samples", y = "Number of unique Peptides") +
  theme(axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"))

ggplot(sat_classII,aes(x = sat_classII$num_samples, y = fit_sat_classII)) +
  geom_smooth(method = "loess", se = F, color = "black",formula = y ~ x) +
  geom_hline(yintercept = abs(fit_param_classII[1])+fit_param_classII[2], linetype = "dashed", alpha = 0.6) +
  labs(x = "Number of Samples", y = "Number of unique Peptides") +
  theme(axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"))










# Waterfall plots---------------------------------------------------------------
# TODO: plot maybe venn into watferall 
# TODO: adjust color for tumor regions
# TODO: save each plot as pdf 

# for loop mit allen combination wie in rna seq analysis 

for (com in apply(combn(names(region_specific_I),2),2, paste, collapse = "_vs_")) {
  
  comparison <- unlist(strsplit(com, "_vs_"))
  df_1 <- do.call(rbind.data.frame,region_specific_I[comparison[1]])
  df_2 <- do.call(rbind.data.frame,region_specific_I[comparison[2]])
  #make waterfall data frame 
  waterfall_df <- make_waterfall_df(df_1,df_2, with_seq = FALSE)
  plot_wf <- plot_waterfall(waterfall_df,
                 comparison[1],
                 comparison[2],
                 com,
                 with_seq = FALSE)
  plot(plot_wf)
  grid::grid.newpage()
  plot_venn_waterfall(df_1$Accessions, df_2$Accessions)
}


# Multiwaterfall----------------------------------------------------------------



#dataset for region exclusive peptides with highest frequency class I##############################
all_waterfall_withoutBEN = make_multi_waterfall_df(NEC_I_uniqe_acc_filterd, T1_I_uniqe_acc_filterd, INF_I_uniqe_acc_filterd)
NEC = all_waterfall_withoutBEN$Nec
T1 = all_waterfall_withoutBEN$T1
INF = all_waterfall_withoutBEN$INF 

all_waterfall_df_withoutBEn = make_multi_waterfall_acc(NEC_I_uniqe_acc_filterd, T1_I_uniqe_acc_filterd, INF_I_uniqe_acc_filterd)
NEC_acc = all_waterfall_df_withoutBEn$Nec
names(NEC_acc)[1]= "Accessions"
T1_acc = all_waterfall_df_withoutBEn$T1
names(T1_acc)[1] = "Accessions"
INF_acc = all_waterfall_df_withoutBEn$INF
names(INF_acc)[1] = "Accessions"


#classII region exclusive peptides with highest frequency 
all_waterfall_classII_without_BEN = make_multi_waterfall_df(NEC_II_uniqe_acc_filterd, T1_II_uniqe_acc_filterd, INF_II_uniqe_acc_filterd)
NEC_II = all_waterfall_classII_without_BEN$Nec
T1_II = all_waterfall_classII_without_BEN$T1
INF_II = all_waterfall_classII_without_BEN$INF 
shared_without_BEN_II_new = all_waterfall_classII_without_BEN$Shared
all_II = all_waterfall_classII_without_BEN$all

all_waterfall_classII_without_BEN_acc = make_multi_waterfall_acc(NEC_II_uniqe_acc_filterd, T1_II_uniqe_acc_filterd, INF_II_uniqe_acc_filterd)
NEC_acc_II = all_waterfall_classII_without_BEN_acc$Nec
names(NEC_acc_II)[1] = "Accessions"

T1_acc_II = all_waterfall_classII_without_BEN_acc$T1
names(T1_acc_II)[1] = "Accessions"

INF_acc_II = all_waterfall_classII_without_BEN_acc$INF
names(INF_acc_II)[1] = "Accessions"


#get Tumor peptide frequencies

uniqe_I = unique(Tumor_I_uniqe_acc_filterd$Sequence)
uniqe_II = unique(Tumor_II_uniqe_acc_filterd$Sequence)

Tumor_I = Tumor_I_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
Tumor_II = Tumor_II_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))

Tumor_I_freq = data.frame(uniqe_I)
Tumor_II_freq = data.frame(uniqe_II)

Tumor_I_freq$Frequency = sapply(uniqe_I, function(seq){
  mean(ifelse(is.na(matchAll(seq, Tumor_I$Sequence)), 0,
              (length(matchAll(seq, Tumor_I$Sequence))/length(unique(Tumor_I$Patient_ID)))*100))
})

Tumor_II_freq$Frequency = sapply(uniqe_II, function(seq){
  mean(ifelse(is.na(matchAll(seq, Tumor_II$Sequence)), 0,
              (length(matchAll(seq, Tumor_II$Sequence))/ length(unique(Tumor_II$Patient_ID))) *100))
})

# get Tumor protein accesion frequencies
uniqe_I_acc = unique(Tumor_I_uniqe_acc_filterd$Accessions)
uniqe_II_acc = unique(Tumor_II_uniqe_acc_filterd$Accessions)

Tumor_I_acc = Tumor_I_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
Tumor_II_acc = Tumor_II_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))

Tumor_I_freq_acc = data.frame(uniqe_I_acc)
Tumor_II_freq_acc = data.frame(uniqe_II_acc)

Tumor_I_freq_acc$Frequency = sapply(uniqe_I_acc, function(acc){
  mean(ifelse(is.na(matchAll(acc, Tumor_I_acc$Accessions)), 0,
              (length(matchAll(acc, Tumor_I_acc$Accessions))/length(unique(Tumor_I_acc$Patient_ID)))*100))
})

Tumor_II_freq_acc$Frequency = sapply(uniqe_II_acc, function(acc){
  mean(ifelse(is.na(matchAll(acc, Tumor_II_acc$Accessions)), 0,
              (length(matchAll(acc, Tumor_II_acc$Accessions))/ length(unique(Tumor_II_acc$Patient_ID))) *100))
})

# GET PEPTIDES AND PROTEIN ACC and write to tsv###########################################
#class I 
NEC_peptide_list = NEC[which(NEC$set1 > 10),1:2]
names(NEC_peptide_list) = c("Sequence", "Peptide_frequency")
T1_peptide_list = T1[which(T1$set2 >10), c(1,3)]
names(T1_peptide_list) = c("Sequence", "Peptide_frequency")
INF_peptide_list = INF[which(INF$set3 > 10), c(1,4)]
names(INF_peptide_list) = c("Sequence", "Peptide_frequency")
Tumor_I_freq = Tumor_I_freq[which(Tumor_I_freq$Frequency > 10),]
names(Tumor_I_freq) = c("Sequence", "Peptide_frequency")

#class II
NEC_peptide_list_II = NEC_II[which(NEC_II$set1 > 10),1:2]
names(NEC_peptide_list_II) = c("Sequence", "Peptide_frequency")
T1_peptide_list_II = T1_II[which(T1_II$set2 >10), c(1,3)]
names(T1_peptide_list_II) = c("Sequence", "Peptide_frequency")
INF_peptide_list_II = INF_II[which(INF_II$set3 > 10), c(1,4)]
names(INF_peptide_list_II) = c("Sequence", "Peptide_frequency")
Tumor_II_freq = Tumor_II_freq[which(Tumor_II_freq$Frequency > 10),]
names(Tumor_II_freq) = c("Sequence", "Peptide_frequency")


names(Tumor_I_freq_acc) = c("Accessions", "Protein_frequency")
names(Tumor_II_freq_acc) = c("Accessions", "Protein_frequency")
#get protein accesions frequencies

NEC_peptide_list_acc = unique(subset(NEC_I_uniqe_acc_filterd[4:5], (NEC_I_uniqe_acc_filterd$Sequence %in% NEC_peptide_list$Sequence)))
T1_peptide_list_acc = unique(subset(T1_I_uniqe_acc_filterd[4:5], (T1_I_uniqe_acc_filterd$Sequence %in% T1_peptide_list$Sequence)))
INF_peptide_list_acc = unique(subset(INF_I_uniqe_acc_filterd[4:5], (INF_I_uniqe_acc_filterd$Sequence %in% INF_peptide_list$Sequence)))
Tumor_I_peptides_acc = unique(subset(Tumor_I_uniqe_acc_filterd[4:5], (Tumor_I_uniqe_acc_filterd$Sequence %in% Tumor_I_freq$Sequence)))

NEC_peptide_list_II_acc = unique(subset(NEC_II_uniqe_acc_filterd[4:5], (NEC_II_uniqe_acc_filterd$Sequence %in% NEC_peptide_list_II$Sequence)))
T1_peptide_list_II_acc = unique(subset(T1_II_uniqe_acc_filterd[4:5], (T1_II_uniqe_acc_filterd$Sequence %in% T1_peptide_list_II$Sequence)))
INF_peptide_list_II_acc = unique(subset(INF_II_uniqe_acc_filterd[4:5], (INF_II_uniqe_acc_filterd$Sequence %in% INF_peptide_list_II$Sequence)))
Tumor_II_peptides_acc = unique(subset(Tumor_II_uniqe_acc_filterd[4:5], (Tumor_II_uniqe_acc_filterd$Sequence %in% Tumor_II_freq$Sequence)))

#merge peptide frequencies and protein freaquencies 

NEC_I_frequend_peptides = merge(NEC_peptide_list, NEC_peptide_list_acc, by = "Sequence", all.x = TRUE)
NEC_I_frequend_peptides = merge(NEC_I_frequend_peptides, NEC_acc[1:2], by = "Accessions", all.x = TRUE)
NEC_I_frequend_peptides$Accessions = getProteinAcc_uniqemappers(NEC_I_frequend_peptides$Accessions)
names(NEC_I_frequend_peptides)[4] = "Protein_frequency"

T1_I_frequend_peptides = merge(T1_peptide_list, T1_peptide_list_acc, by = "Sequence", all.x = TRUE)
T1_I_frequend_peptides = merge(T1_I_frequend_peptides, T1_acc[c(1,3)], by = "Accessions", all.x = TRUE)
T1_I_frequend_peptides$Accessions = getProteinAcc_uniqemappers(T1_I_frequend_peptides$Accessions)
names(T1_I_frequend_peptides)[4] = "Protein_frequency"

INF_I_frequend_peptides = merge(INF_peptide_list, INF_peptide_list_acc, by = "Sequence", all.x = TRUE)
INF_I_frequend_peptides = merge(INF_I_frequend_peptides, INF_acc[c(1,4)], by = "Accessions", all.x = TRUE)
INF_I_frequend_peptides$Accessions = getProteinAcc_uniqemappers(INF_I_frequend_peptides$Accessions)
names(INF_I_frequend_peptides)[4] = "Protein_frequency"

NEC_II_frequend_peptides = merge(NEC_peptide_list_II, NEC_peptide_list_II_acc, by = "Sequence", all.x = TRUE)
NEC_II_frequend_peptides = merge(NEC_II_frequend_peptides, NEC_acc_II[1:2], by = "Accessions", all.x = TRUE)
NEC_II_frequend_peptides$Accessions = getProteinAcc_uniqemappers(NEC_II_frequend_peptides$Accessions)
names(NEC_II_frequend_peptides)[4] = "Protein_frequency"


T1_II_frequend_peptides = merge(T1_peptide_list_II, T1_peptide_list_II_acc, by = "Sequence", all.x = TRUE)
T1_II_frequend_peptides = merge(T1_II_frequend_peptides, T1_acc_II[c(1,3)], by = "Accessions", all.x = TRUE)
T1_II_frequend_peptides$Accessions = getProteinAcc_uniqemappers(T1_II_frequend_peptides$Accessions)
names(T1_II_frequend_peptides)[4] = "Protein_frequency"

INF_II_frequend_peptides = merge(INF_peptide_list_II, INF_peptide_list_II_acc, by = "Sequence", all.x = TRUE)
INF_II_frequend_peptides = merge(INF_II_frequend_peptides, INF_acc_II[c(1,4)], by = "Accessions", all.x = TRUE)
INF_II_frequend_peptides$Accessions = getProteinAcc_uniqemappers(INF_II_frequend_peptides$Accessions)
names(INF_II_frequend_peptides)[4] = "Protein_frequency"

Tumor_I_frequent_peptides = merge(Tumor_I_freq, unique(Tumor_I_uniqe_acc_filterd[4:5]), by = "Sequence", all.x = TRUE)
Tumor_I_frequent_peptides = merge(Tumor_I_frequent_peptides, Tumor_I_freq_acc, by = "Accessions", all.x = TRUE )
Tumor_I_frequent_peptides$Accessions = getProteinAcc_uniqemappers(Tumor_I_frequent_peptides$Accessions)

Tumor_II_frequent_peptides = merge(Tumor_II_freq, unique(Tumor_II_uniqe_acc_filterd[4:5]), by = "Sequence", all.x = TRUE)
Tumor_II_frequent_peptides = merge(Tumor_II_frequent_peptides, Tumor_II_freq_acc, by = "Accessions", all.x = TRUE )
Tumor_II_frequent_peptides$Accessions = getProteinAcc_uniqemappers(Tumor_II_frequent_peptides$Accessions)


write.table(NEC_I_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_frequend_peptides.tsv")
write.table(T1_I_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_frequend_peptides.tsv")
write.table(INF_I_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_frequend_peptides.tsv")
write.table(Tumor_I_frequent_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_I_frequend_peptides.tsv")
write.table(NEC_II_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_frequend_peptides.tsv")
write.table(T1_II_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_frequend_peptides.tsv")
write.table(INF_II_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_frequend_peptides.tsv")
write.table(Tumor_II_frequent_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_II_frequend_peptides.tsv")





# get protein acc for peptids and export it as csv or tsf 
write_region_pepacc_df <- function(region_peptide_list, region_excl_acc, filename){
  names(region_peptide_list)[1] = "Sequence"
  names(region_peptide_list)[2] = "Frequency"
  subdf = subset(region_excl_acc, (region_excl_acc$Sequence %in% region_peptide_list$Sequence))
  reduced_subdf = unique(subdf[,4:5])
  reduced_subdf$Accessions = getProteinAcc(reduced_subdf$Accessions)
  jointdf = merge(reduced_subdf, region_peptide_list, by = "Sequence")
  write.table(jointdf, sep = "\t", col.names = TRUE, file = paste("/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/", filename, sep = "" ), quote = FALSE, row.names = FALSE) 
  write.table(jointdf, sep = "\t", col.names = TRUE, file = paste("/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/", filename, sep = "" ), quote = FALSE, row.names = FALSE) 
  
  return(jointdf)
}

#classI
NEC_excl_acc_pep = write_region_pepacc_df(NEC_peptide_list, NEC_exclusive_I_uniqeAcc_filterd, "NEC_excl_pepacc.tsv")
T1_excl_acc_pep = write_region_pepacc_df(T1_peptide_list, T1_exclusive_I_uniqeAcc_filterd, "T1_excl_pepacc.tsv")
INF_excl_acc_pep = write_region_pepacc_df(INF_peptide_list, INF_exclusive_I_uniqeAcc_filterd, "INF_excl_pepacc.tsv")
shared_without_ben_acc_pep = write_region_pepacc_df(shared_withoutBEN_list, classI_uniqe_mappers, "shared_without_ben.tsv")
all_acc_pep = write_region_pepacc_df(all_I_list, classI_uniqe_mappers,"all_I_pepacc.tsv")
Tumor_I_acc_pep = write_region_pepacc_df(Tumor_I_freq, Tumor_I_uniqe_acc_filterd, "Tumor_I_pepacc.tsv")

#calssII
NEC_excl_acc_pep_II = write_region_pepacc_df(NEC_peptide_list_II, NEC_exclusive_II_uniqeAcc_filterd, "NEC_excl_pepacc_II.tsv")
T1_excl_acc_pep_II = write_region_pepacc_df(T1_peptide_list_II, T1_exclusive_II_uniqeAcc_filterd, "T1_excl_pepacc_II.tsv")
INF_excl_acc_pep_II = write_region_pepacc_df(INF_peptide_list_II, INF_exclusive_II_uniqeAcc_filterd, "INF_excl_pepacc_II.tsv")
shared_without_ben_acc_pep_II = write_region_pepacc_df(shared_withoutBEN_list_II_new, classII_uniqe_mappers, "shared_without_ben_II.tsv")
all_II_acc_pep = write_region_pepacc_df(all_II_list, classII_uniqe_mappers,"all_II_pepacc.tsv")

new_seq_II = c()
seq_II = Tumor_II_freq$uniqe_II
for(i in seq_II){
  #print(nchar(i))
  if(nchar(i) > 8){
    new_seq_II = append(new_seq_II, i )
  }
}

Tumor_II_freq = Tumor_II_freq %>% dplyr::filter(Tumor_II_freq$uniqe_II %in% new_seq_II)
table(Tumor_II_freq$uniqe_II)
Tumor_II_acc_pep = write_region_pepacc_df(Tumor_II_freq, Tumor_II_uniqe_acc_filterd, "Tumor_II_pepacc.tsv")


for( i in Tumor_II_acc_pep$Sequence){
  print(nchar(i))
}


check_len_II <- function(list){
  new_list = c()
  for(i in list){
    
    if(nchar(i) > 8){
      new_list = append(new_list, i )
    }
  }
  return(new_list)
}

#read in data for gioele 
NEC_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_excl_pepacc.tsv")
T1_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_excl_pepacc.tsv")
INF_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_excl_pepacc.tsv")
shared_without_ben_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/shared_without_ben.tsv")
all_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/all_I_pepacc.tsv")


NEC_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ NEC_excl_pepacc_II.tsv")
T1_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ T1_excl_pepacc_II.tsv")
INF_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ INF_excl_pepacc_II.tsv")
shared_without_ben_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ shared_without_ben_II.tsv")
all_II_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ all_II_pepacc.tsv")


#wirte only sequences to netMHCpan for class I 

write.table(NEC_excl_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/NEC_excl_peptides.tsv")
write.table(T1_excl_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/T1_excl_peptides.tsv")
write.table(INF_excl_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/INF_excl_peptides.tsv")
#write.table(shared_without_ben_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/shared_peptides.tsv")
write.table(Tumor_I_uniqe_acc_filterd$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/TumorI_peptides.tsv")
write.table(Tumor_II_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/Tumor_regions//TumorII_peptides.tsv")

write.table(NEC_excl_acc_pep_II$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/tumor_regions/NEC_excl_peptides_II.tsv")
write.table(check_len_II(T1_excl_acc_pep_II$Sequence),sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/tumor_regions/T1_excl_peptides_II.tsv")
write.table(INF_excl_acc_pep_II$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/tumor_regions/INF_excl_peptides_II.tsv")


#add tumor regions to Tumor_I and Tumor_II dfs 

nec_seq = NEC_I_uniqe_acc_filterd[3:4]
T1_seq = T1_I_uniqe_acc_filterd[3:4]
inf_seq = INF_I_uniqe_acc_filterd[3:4]
ben_seq = BEN_I_uniqe_acc_filterd[3:4]


nec_II_seq = NEC_II_uniqe_acc_filterd[3:4]
T1_II_seq = T1_II_uniqe_acc_filterd[3:4]
inf_II_seq = INF_II_uniqe_acc_filterd[3:4]
ben_II_seq = BEN_II_uniqe_acc_filterd[3:4]



write.table(nec_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_seq.tsv")
write.table(T1_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_seq.tsv")
write.table(inf_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_seq.tsv")
write.table(ben_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/BEN_I_seq.tsv")
write.table(nec_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_seq.tsv")
write.table(T1_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_seq.tsv")
write.table(inf_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_seq.tsv")
write.table(ben_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/BEN_II_seq.tsv")


test_join = unique(list(Tumor_I_acc_pep,nec_seq,T1_seq,inf_seq,ben_seq) %>% reduce(left_join, by = "Sequence"))
test_II_join = unique(list(Tumor_II_acc_pep,nec_II_seq,T1_II_seq,inf_II_seq,ben_II_seq) %>% reduce(left_join, by = "Sequence"))



library("purrr")

Tumor_I_acc_pep_o = merge(Tumor_I_acc_pep, nec_seq, by = "Sequence", all.x=TRUE)

Tumor_I_acc_pep_o  = unique(Tumor_I_acc_pep_o )

#PEPTIDE SELECTION FOR GIOELE -----------------------------------------------------------------------
#Tumor-exclusive peptides, presented on at least 2 patients --> xy-excl_acc_pep dfs enthalten das 



#Sorted based on ascending order of protein frequency in healthy tissue database --> in wie viele healthy tissues kommt es vor ?
#https://www.proteinatlas.org/about/download-->normal tissue data 
library("EnsDb.Hsapiens.v86")
library(ensembldb)

#read in normal tissue data
#Expression profiles for proteins in human tissues based on immunohistochemisty using tissue micro arrays

normal_tissue_proteom = read.csv("/Users/cschwitalla/Downloads/normal_tissue.tsv", header = TRUE, sep = "\t")
# filter for level--> all not detected will be removed

normal_tissue_proteom_filter = normal_tissue_proteom[-which(normal_tissue_proteom$Level == "Not detected"),]

normal_tissue_proteom_filter2 = normal_tissue_proteom_filter[-which(normal_tissue_proteom_filter$Reliability == "Uncertain"),]

normal_tissue_proteom_filter3 = normal_tissue_proteom_filter2 %>% dplyr::filter(normal_tissue_proteom_filter2$Gene.name %in% all_acc_pep$GeneName)
#genename tp protein acc
edb <- EnsDb.Hsapiens.v86
hasProteinData(edb)

#select(edb, keys = "P15502", keytype = "UNIPROTID", columns = "GENENAME")[1,][2]

all_acc_pep$GeneName = sapply(all_acc_pep$Accessions, function(gen){
  genname = select(edb, keys = gen, keytype = "UNIPROTID", columns = "GENENAME")[1,][2]
  
})



#normal_tissue_proteom_filter3 = normal_tissue_proteom_filter2 %>% filter(normal_tissue_proteom_filter2$)
#compute tissue abundance 
#3 = df3 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))

normal_tissue_df = normal_tissue_proteom_filter3 %>% group_by(Tissue) %>% summarise(Gene.name = unique(Gene.name))
uniq_genes = unique(normal_tissue_df$Gene.name)
normal_tissue_proteinf_df = data.frame(unique(normal_tissue_df$Gene.name))
# ratio of peptide occurences für NEC und INF 
normal_tissue_proteinf_df$Frequnecy = sapply(uniq_genes, function(gen){
  mean(ifelse(is.na(matchAll(gen, normal_tissue_df$Gene.name)), 0,
              (length(matchAll(gen, normal_tissue_df$Gene.name))/length(unique(normal_tissue_df$Tissue)))*100))
})

unique(normal_tissue_df$Tissue)
normal_tissue_df_protein_abundance = normal_tissue_proteinf_df

names(normal_tissue_df_protein_abundance)[1] = "GeneName"
names(normal_tissue_df_protein_abundance)[2] = "NormalTissueFrequency"

#combine normal tissue frequency 

all_acc_pep_2 = merge(all_acc_pep, normal_tissue_df_protein_abundance, by = "GeneName")


gNames = do.call(rbind.data.frame, all_acc_pep$GeneName)
names(gNames)[1] = "GeneName"
all_acc_pep$GeneName = gNames$GeneName

#Descending order of peptide frequency in glioblastoma tissue samples

#Descending order of protein source frequency in glioblastoma tissue samples

#protein frequency in non-glioblastoma malignant samples in descending order












#GET REGION EXCLUSIVE DATA FRAMES FOR INTEGRATION---------------------------------------------------------------
# files of tumor region exclusive peptides source protein accessions
# for integration of the different omics types

# veraltet NEU MACHEN UNBEDINGT !!======================
write.csv(Nec_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//NEC_excl_acc.csv", row.names = FALSE, quote = FALSE)
write.csv(T1_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//T1_excl_acc.csv", row.names = FALSE, quote = FALSE)
write.csv(INF_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//INF_excl_acc.csv", row.names = FALSE, quote = FALSE)
write.csv(BEN_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//BEN_excl_acc.csv", row.names = FALSE, quote = FALSE)
#=============================


test = t(as.data.frame(strsplit(files_cI, "_", fixed = TRUE)))
table(test)
