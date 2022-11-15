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
output_dir <- ("/Users/cschwitalla/Documents/Ligandomics_analysis/")
#metadata file
exome_metadata_file <- "/Users/cschwitalla/Documents/WES_analysis/WES_metadata.tsv"
# import functions for the analysis
source("/Users/cschwitalla/git/students/cschwitalla/Ligandomics_Analysis/functions_ligandomics.R")
################################################################################
###                             Load libraries                               ###
################################################################################
required_libs <- c(
  "tidyr", "readxl", "ggVennDiagram", "dplyr", "stringr",
  "tibble", "ggplot2", "org.Hs.eg.db", "RColorBrewer", "EnsDb.Hsapiens.v86"
)

suppressMessages(invisible(lapply(required_libs, library, character.only = TRUE)))


################################################################################
###                            Load Data                                     ###
################################################################################
# load patients HLA-typing data, important for predicting binding affinity
# with netMHCpan later
GB_HLA_types <- read_xlsx(paste0(input_dir, "HLA-Typisierung_GBM.xlsx"),
  col_names = TRUE
)
# get list of unique HLA types
uniqe_HLA_types <- unique(c(as.matrix(GB_HLA_types[2:16, 2:7])))
# process HLA-typing data
new_HLA_types <- rewrite_HLA_types(uniqe_HLA_types, as_string = FALSE)
#load metadata
exome_metadata <- read.table(file = exome_metadata_file, sep = "\t", header = TRUE)


# Benign data from inhouse database of the immunology department tuebingen------
benign_pep_I <- read.csv(paste0(input_dir, "newBenignmorespecific/Benign_class1.csv"),
  header = FALSE,
  sep = ","
)[, 1:2]
benign_pep_II <- read.csv(paste0(input_dir, "newBenignmorespecific/Benign_class2.csv"),
  header = FALSE,
  sep = ","
)[, 1:2]


# Data from the HLA-Ligand Atlas database---------------------------------------
# read in dataframes from HLA ligand atlas
HLA_ligand_atlas_pep <- read.csv(paste0(input_dir, "hla_2020.12/HLA_aggregated.tsv"),
  header = TRUE,
  sep = "\t"
)
HLA_ligand_atlas_acc <- read.csv(paste0(input_dir, "hla_2020.12/HLA_protein_map.tsv"),
  header = TRUE,
  sep = "\t"
)
# aggregate HLA ligand atlas dataframe befor merging
HLA_ligand_atlas_acc <- HLA_ligand_atlas_acc %>%
  dplyr::group_by(peptide_sequence_id) %>%
  dplyr::summarise(acc = toString(uniprot_id))

# merge dataframes by peptide sequence id to add protein accession numbers to the dataframe
HLA_atlas_data <- merge(HLA_ligand_atlas_pep, HLA_ligand_atlas_acc, by = "peptide_sequence_id")
# split the ligand atlas according to the hla class ( 1 or 2 )
HLA_atlas_data <- split(HLA_atlas_data, HLA_atlas_data$hla_class)
# get protein acc frequency of source proteins from the hla ligand atlas for
# integrating the frequencies to the peptide selection
class_I_HLA_atlas_acc <- separate_rows(data.frame(acc = c(HLA_atlas_data$`HLA-I`$acc, HLA_atlas_data$`HLA-I+II`$acc)), acc, sep = ",")
class_II_HLA_atlas_acc <- separate_rows(data.frame(acc = c(HLA_atlas_data$`HLA-II`$acc, HLA_atlas_data$`HLA-I+II`$acc)), acc, sep = ",")

HLA_atlas_protein_f_cI <- data.frame(table(class_I_HLA_atlas_acc))
names(HLA_atlas_protein_f_cI) <- c("Accessions", "abundance")
HLA_atlas_protein_f_cI$percentage <- HLA_atlas_protein_f_cI$abundance / length(unique(HLA_atlas_protein_f_cI$Accessions)) * 100
HLA_atlas_protein_f_cII <- data.frame(table(class_II_HLA_atlas_acc))
names(HLA_atlas_protein_f_cII) <- c("Accessions", "abundance")
HLA_atlas_protein_f_cII$percentage <- HLA_atlas_protein_f_cII$abundance / length(unique(HLA_atlas_protein_f_cII$Accessions)) * 100

# GB tumor region ligandomics data----------------------------------------------
# get all filenames with path of the raw data
files_cI <- dir(paste0(input_dir, "ligandomics_results_classI_8-12"),
  recursive = TRUE,
  pattern = ".tsv",
  full.names = TRUE
)
files_cII <- dir(paste0(input_dir, "ligandomics_results_classII_8-30"),
  recursive = TRUE,
  pattern = ".tsv",
  full.names = TRUE
)

# create dataframes from raw results
classI_df <- create_data_frame(files_cI)
classII_df <- create_data_frame(files_cII)
uniqI <- unique(classI_df)
################################################################################
###                          Data preparation                                ###
################################################################################
# merge HLA-Ligand Atlas data and benign data from inhouse immu. data base
# to create a comined data set of known benign peptides
combi_benign_pep_I <- c(
  HLA_atlas_data$`HLA-I`$peptide_sequence,
  HLA_atlas_data$`HLA-I+II`$peptide_sequence,
  benign_pep_I$V1
)
combi_benign_pep_II <- c(
  HLA_atlas_data$`HLA-II`$peptide_sequence,
  HLA_atlas_data$`HLA-I+II`$peptide_sequence,
  benign_pep_II$V1
)
# filter GB data so that the known benign peptides are excluded
classI_df <- subset(classI_df, !(Sequence %in% combi_benign_pep_I))
#classI_df_f <- subset(classI_df, !(Sequence %in% combi_benign_pep_I))
classII_df <- subset(classII_df, !(Sequence %in% combi_benign_pep_II))

# remove peptides that were predicted to originate from more than one
# source protein in the GB data
classI_df <- classI_df[-grep(";", classI_df$Accessions), ]
#classI_df_ff <- classI_df_f[-grep(";", classI_df_f$Accessions), ]

classII_df <- classII_df[-grep(";", classII_df$Accessions), ]
# rewrite protein accessions in GB data to uniprot-ids
classI_df$Accessions <- get_protein_acc(classI_df$Accessions)
#classI_df_ff$Accessions <- get_protein_acc(classI_df_ff$Accessions)
classII_df$Accessions <- get_protein_acc(classII_df$Accessions)

# create dataframe mappin uniprot ids to gene names
mapping_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
  keys = c(classI_df$Accessions, classII_df$Accessions),
  columns = "GENENAME",
  keytype = "UNIPROTID"
)


test <- unique(classI_df_fff)
test2 <- unique(classI_df)

# add the gene names to the respective uniprot ids
classI_df <- merge(classI_df, mapping_df, by.x = "Accessions", by.y = "UNIPROTID", all.x = TRUE)
#classI_df_fff <- merge(classI_df_ff, mapping_df, by.x = "Accessions", by.y = "UNIPROTID", all.x = TRUE)

classII_df <- merge(classII_df, mapping_df, by = "Accessions", by.y = "UNIPROTID", all.x = TRUE)
# peptide frequency 
classI_df <- calc_peptide_frequency(classI_df)
classII_df <- calc_peptide_frequency(classII_df)
# split GB data to get tumor region specific peptide data
region_specific_I <- split(classI_df, classI_df$Tumor_region)
region_specific_II <- split(classII_df, classII_df$Tumor_region)

#split by patients because for netMHCpan binding prediction 
patient_specific_I <- split(classI_df, classI_df$Patient_ID)

# TODO: save dataframes to tsv
################################################################################
##                               Analysis                                     ##
################################################################################


# Venn diagrams ----------------------------------------------------------------
# extract venn diagrams region exclusive peptides than
set_I <- setNames(
  vector("list", length = length(names(region_specific_I))),
  c(names(region_specific_I))
)

set_II <- setNames(
  vector("list", length = length(names(region_specific_II))),
  c(names(region_specific_II))
)

for (i in region_specific_I) {
  set_I[i$Tumor_region[1]] <- list(i$Sequence)
}
for (i in region_specific_II) {
  set_II[i$Tumor_region[1]] <- list(i$Sequence)
}

venn_data_I <- plot_custom_venn(set_I, "HLA class I Peptides")
venn_data_II <- plot_custom_venn(set_II, "HLA class II Peptides")

# Venn of most frequent pep tides
classI_df_frequent <- classI_df[classI_df$freq >15,]
classII_df_frequent <- classII_df[classII_df$freq >15,]
classI_df_frequent_regions <- split(classI_df_frequent, classI_df_frequent$Tumor_region)
classII_df_frequent_regions <- split(classII_df_frequent, classII_df_frequent$Tumor_region)

set_I_f <- setNames(
  vector("list", length = length(names(classI_df_frequent_regions))),
  c(names(classI_df_frequent_regions))
)
set_II_f <- setNames(
  vector("list", length = length(names(classII_df_frequent_regions))),
  c(names(classII_df_frequent_regions))
)
for (i in classI_df_frequent_regions) {
  set_I_f[i$Tumor_region[1]] <- list(i$Sequence)
}
for (i in classII_df_frequent_regions) {
  set_II_f[i$Tumor_region[1]] <- list(i$Sequence)
}

set_I_II <- list(Class_I = classI_df_frequent$Sequence, ClassII = classII_df_frequent$Sequence)
plot_custom_venn(set_I_II, "HLA class I+II frequent peptides")

venn_data_I_f <- plot_custom_venn(set_I_f, "HLA class I frequent peptides")
venn_data_II_f <- plot_custom_venn(set_II_f, "HLA class II frequent peptides")
grid.arrange(venn_data_I_f, venn_data_II_f,
             ncol = 2)
#length(unique(classI_df_frequent$Sequence))
#length(unique(classII_df_frequent$Sequence))

# Length distribution ----------------------------------------------------------
# TODO: adjust color for tumor regions
# TODO: adjust output that plot is saved as pdf in output dir

par(mfrow = c(2,2))
ln1 <- plot_length_distribution(classI_df_fff,
  as_bar = TRUE,
  "HLA class I peptides"
)

ln2 <- plot_length_distribution(classII_df_fff,
                         as_bar = TRUE,
                         "HLA class II peptides"
)
library(gridExtra)
library(grid)
grid.arrange(ln1, ln2,
             ncol = 2)


# upset region + frequencies----------------------------------------------------
library("UpSetR")
library("ComplexUpset")
frequency_specific_I <- split(classI_df_frequent, classI_df_frequent$freq)
set_pep_I <- list(NEC = classI_df_frequent_regions$NEC$Sequence,
                    T1 = classI_df_frequent_regions$T1$Sequence,
                    INF = classI_df_frequent_regions$INF$Sequence,
                    BEN = classI_df_frequent_regions$BEN$Sequence)
                  #   #"1_Patient" = frequency_specific_I$`7.69230769230769`$Sequence,
                  #   "2_Patients" = frequency_specific_I$`15.3846153846154`$Sequence,
                  #   "3_Patients" =frequency_specific_I$`23.0769230769231`$Sequence,
                  # "4_Patients" = frequency_specific_I$`30.7692307692308`$Sequence,
                  # "5_Patients" =frequency_specific_I$`38.4615384615385`$Sequence,
                  # "6_Patients" =frequency_specific_I$`46.1538461538462`$Sequence,
                  # "7_Patients" =frequency_specific_I$`53.8461538461538`$Sequence,
                  # "9_Patients" =frequency_specific_I$`69.2307692307692`$Sequence
                  # )

classI_df_frequent_regions$NEC[,c("Sequence", freq)]
# "2_Patients" = frequency_specific_I$`15.3846153846154`$Sequence,
# "3_Patients" =frequency_specific_I$`23.0769230769231`$Sequence,
# "4_Patients" = frequency_specific_I$`30.7692307692308`$Sequence,
# "5_Patients" =frequency_specific_I$`38.4615384615385`$Sequence,
# "6_Patients" =frequency_specific_I$`46.1538461538462`$Sequence,
# "7_Patients" =frequency_specific_I$`53.8461538461538`$Sequence,
# "9_Patients" =frequency_specific_I$`69.2307692307692`$Sequence
set_pep_I <- list(NEC = classI_df_frequent_regions$NEC[c(Sequence, freq)],
                  T1 = classI_df_frequent_regions$T1$Sequence,
                  INF = classI_df_frequent_regions$INF$Sequence,
                  BEN = classI_df_frequent_regions$BEN$Sequence)
#length(set_all_all)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), header = T, sep = ";")
#classI
peptides <- unique(c(classI_df_frequent_regions$NEC$Sequence, classI_df_frequent_regions$T1$Sequence, classI_df_frequent_regions$INF$Sequence, classI_df_frequent_regions$BEN$Sequence))

#classII
clII_peptides <- unique(c(classII_df_frequent_regions$NEC$Sequence, classII_df_frequent_regions$T1$Sequence, classII_df_frequent_regions$INF$Sequence, classII_df_frequent_regions$BEN$Sequence))

#class I
upset_df_own <- data_frame(peptides)
upset_df_own$BEN <-ifelse((c(upset_df_own$peptides) %in% classI_df_frequent_regions$BEN$Sequence),1,0)
upset_df_own$NEC <-ifelse((c(upset_df_own$peptides) %in% classI_df_frequent_regions$NEC$Sequence),1,0)
upset_df_own$T1 <-ifelse((c(upset_df_own$peptides) %in% classI_df_frequent_regions$T1$Sequence),1,0)
upset_df_own$INF <-ifelse((c(upset_df_own$peptides) %in% classI_df_frequent_regions$INF$Sequence),1,0)
upset_df_own <- merge(upset_df_own, classI_df_frequent[,c("Sequence", "freq")], by.x = "peptides", by.y = "Sequence", all.x = TRUE)
upset_df_own <- unique(upset_df_own)

#classII 
upset_df_II <- data_frame(clII_peptides)
upset_df_II$BEN <-ifelse((c(upset_df_II$clII_peptides) %in% classII_df_frequent_regions$BEN$Sequence),1,0)
upset_df_II$NEC <-ifelse((c(upset_df_II$clII_peptides) %in% classII_df_frequent_regions$NEC$Sequence),1,0)
upset_df_II$T1 <-ifelse((c(upset_df_II$clII_peptides) %in% classII_df_frequent_regions$T1$Sequence),1,0)
upset_df_II$INF <-ifelse((c(upset_df_II$clII_peptides) %in% classII_df_frequent_regions$INF$Sequence),1,0)
upset_df_II <- merge(upset_df_II, classII_df_frequent[,c("Sequence", "freq")], by.x = "clII_peptides", by.y = "Sequence", all.x = TRUE)
upset_df_II <- unique(upset_df_II)




#row.names(upset_df_own) <- upset_df_own$peptides
upset_df_own2 <- upset_df_own
upset_df_own2 <- upset_df_own2[order(upset_df_own2$freq,decreasing = TRUE),]
#upset_df_own2$peptides <- NULL
#colnames(upset_df_own2) <- c("NEC", "BEN", "T1", "INF", "Frequency")
# mach die richtige column class I 
upset_df_own2$freq <- round(upset_df_own2$freq, digits = 2)
upset_df_own2$freq <- as.character(upset_df_own2$freq)
upset_df_own2$patients <- NA
upset_df_own2$patients <- ifelse((c(upset_df_own2$freq) == "15.38"),"2 Patients",upset_df_own2$patients)
upset_df_own2$patients <- ifelse((c(upset_df_own2$freq) == "23.08"),"3 Patients",upset_df_own2$patients)
upset_df_own2$patients <- ifelse((c(upset_df_own2$freq) == "30.77"),"4 Patients",upset_df_own2$patients)
upset_df_own2$patients <- ifelse((c(upset_df_own2$freq) == "38.46"),"5 Patients",upset_df_own2$patients)
upset_df_own2$patients <- ifelse((c(upset_df_own2$freq) == "46.15"),"6 Patients",upset_df_own2$patients)
upset_df_own2$patients <- ifelse((c(upset_df_own2$freq) == "53.85"),"7 Patients",upset_df_own2$patients)
upset_df_own2$patients <- ifelse((c(upset_df_own2$freq) == "69.23"),"9 Patients",upset_df_own2$patients)

#class II 
upset_df_II$freq <- round(upset_df_II$freq, digits = 2)
upset_df_II$freq <- as.character(upset_df_II$freq)
upset_df_II$patients <- NA
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "15.38"),"2 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "23.08"),"3 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "30.77"),"4 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "38.46"),"5 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "46.15"),"6 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "53.85"),"7 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "61.54"),"8 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "69.23"),"9 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "76.92"),"10 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "84.62"),"11 Patients",upset_df_II$patients)
upset_df_II$patients <- ifelse((c(upset_df_II$freq) == "92.31"),"12 Patients",upset_df_II$patients)


# UpSetR::upset(classI_df_frequent, order.by = "freq", nsets = length(set_pep_I),
#               sets = c("NEC", "T1", "INF", "BEN" ),#,"2_Patients","3_Patients","4_Patients","5_Patients","6_Patients","7_Patients","9_Patients"),
#               keep.order = TRUE,
#               sets.bar.color = c("#9e0142", "#fdae61", "#74add1" , "#4B6C22"),
#               queries =list(list(query = intersects, params = list("NEC"), color = "#9e0142", active = F)
#                             ,list(query = elements, params=list(frequency_specific_I$`15.3846153846154`$Sequence), color = "red")))#,"black","black", "black","black", "black" ,"black","black"))
# 

test_f <- function (row, freq) {
  data <- (row["freq"] == freq)
  
}

ComplexUpset::upset(upset_df_own2, c("NEC", "T1", "INF", "BEN"),
                    base_annotations = list(
                      "Intersection size" = intersection_size(
                        mapping = aes(fill = patients)
                      ) + scale_fill_manual(values = c(
                        "2 Patients" = "#F49B67",
                        "3 Patients" = "#F28140",
                        "4 Patients" = "#EA5D4D",
                        "5 Patients" = "#E3436B",
                        "6 Patients" = "#B3337A",
                        "7 Patients" =  "#832388",
                        "9 Patients" =  "#5E1961"))),
                    queries = list(   upset_query(set = "NEC", fill = "#9e0142"),
                      upset_query(set = "T1", fill = "#fdae61"),
                      upset_query(set = "INF", fill = "#74add1"),
                      upset_query(set = "BEN", fill = "#4B6C22")),
                    
                    )

ComplexUpset::upset(upset_df_own2, c("NEC", "T1", "INF", "BEN"),
                    base_annotations = list(
                      "Intersection size" = intersection_size(
                        mapping = aes(fill = patients)
                      ) + scale_fill_manual(values = c(
                        "9 Patients" =  "#5E1961",
                        "7 Patients" =  "#832388",
                        "6 Patients" = "#B3337A",
                        "5 Patients" = "#E3436B",
                        "4 Patients" = "#EA5D4D",
                        "3 Patients" = "#F28140",
                        "2 Patients" = "#F49B67"
                    
                        ))),
                    queries = list(   upset_query(set = "NEC", fill = "#9e0142"),
                                      upset_query(set = "T1", fill = "#fdae61"),
                                      upset_query(set = "INF", fill = "#74add1"),
                                      upset_query(set = "BEN", fill = "#4B6C22")),
                    matrix = 
                      intersection_matrix(geom = geom_point(shape ="circle filled", size =3))+
                                            scale_color_manual(
                                              values = c("NEC" = "#9e0142", 
                                                         "T1" = "#fdae61",
                                                         "INF"= "#74add1",
                                                         "BEN" = "#4B6C22")
                                          
                                            ))
#classII
ComplexUpset::upset(upset_df_II, c("NEC", "T1", "INF", "BEN"),
                    base_annotations = list(
                      "Intersection size" = intersection_size(
                        mapping = aes(fill = patients)
                      ) + scale_fill_manual(values = c(
                       "12 Patients" =  "#000000",
                        "11 Patients" =  "#311864",
                        "10 Patients" =  "#4F1C60",
                        "9 Patients" =  "#6F1D72",
                        "8 Patients" =  "#832388",
                        "7 Patients" =  "#B3337A",
                        "6 Patients" = "#E3436B",
                        "5 Patients" = "#EA5D4D",
                        "4 Patients" = "#F28140",
                        "3 Patients" = "#F49B67",
                        "2 Patients" = "#F1B676"
                        
                      ))),
                    queries = list(   upset_query(set = "NEC", fill = "#9e0142"),
                                      upset_query(set = "T1", fill = "#fdae61"),
                                      upset_query(set = "INF", fill = "#74add1"),
                                      upset_query(set = "BEN", fill = "#4B6C22")),
                    matrix = 
                      intersection_matrix(geom = geom_point(shape ="circle filled", size =3))+
                      scale_color_manual(
                        values = c("NEC" = "#9e0142", 
                                   "T1" = "#fdae61",
                                   "INF"= "#74add1",
                                   "BEN" = "#4B6C22")
                        
                      ))


plot_custom_venn(set_pep_I, "t")
 UpSetR::upset(upset_df_own2, order.by = "freq",queries = list(
  list(query = test_f, params = list( c(15.38462)), color = "#F28E54", active= T),
  list(query = test_f, params = list( c(23.07692)), color = "#F28140", active= T),
  list(query = test_f, params = list( c(30.76923)), color = "#EA5D4D", active= T),
  list(query = test_f, params = list( c(38.46154)), color = "#E3436B", active= T),
  list(query = test_f, params = list( c(46.15385)), color = "#B3337A", active= T),
  list(query = test_f, params = list( c(53.84615)), color = "#832388", active= T),
  list(query = test_f, params = list( c(69.23077)), color = "#5E1961", active= F)))

UpSetR::upset(upset_df_own2, order.by = "freq")

UpSetR::upset(fromList(set_pep_I), order.by = "freq", nsets = length(set_pep_I),
              sets = c("NEC", "T1", "INF", "BEN" ),#,"2_Patients","3_Patients","4_Patients","5_Patients","6_Patients","7_Patients","9_Patients"),
              keep.order = TRUE,
              sets.bar.color = c("#9e0142", "#fdae61", "#74add1" , "#4B6C22"),
              queries =list(list(query = intersects, params = list("NEC"), color = "#9e0142", active = F)
              ,list(query = elements, params=list(frequency_specific_I$`15.3846153846154`$Sequence), color = "red")))#,"black","black", "black","black", "black" ,"black","black"))


# queries = list(
#   upset_query(set = "NEC", fill = "#9e0142"),
#   upset_query(set = "T1", fill = "#fdae61"),
#   upset_query(set = "INF", fill = "#74add1"),
#   upset_query(set = "BEN", fill = "#4B6C22")

upset(set_pep_I,)

# saturation analysis-----------------------------------------------------------
# TODO: make functions here
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

fit_param_classI <- NLSstAsymptotic(sortedXyData(c(0, sat_classI$num_samples), c(0, sat_classI$num_unique)))
fit_sat_classI <- SSasymp(sat_classI$num_samples, abs(fit_param_classI[1]) + fit_param_classI[2], fit_param_classI[1], fit_param_classI[3])

fit_param_classII <- NLSstAsymptotic(sortedXyData(c(0, sat_classII$num_samples), c(0, sat_classII$num_unique)))
fit_sat_classII <- SSasymp(sat_classII$num_samples, abs(fit_param_classII[1]) + fit_param_classII[2], fit_param_classII[1], fit_param_classII[3])

fit_param_classI

# plot
ggplot(sat_classI, aes(x = sat_classI$num_samples, y = fit_sat_classI)) +
  geom_smooth(method = "loess", se = F, color = "black", formula = y ~ x) +
  geom_hline(yintercept = abs(fit_param_classI[1]) + fit_param_classI[2], linetype = "dashed", alpha = 0.6) +
  labs(x = "Number of Samples", y = "Number of unique Peptides") +
  theme(
    axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent")
  )

ggplot(sat_classII, aes(x = sat_classII$num_samples, y = fit_sat_classII)) +
  geom_smooth(method = "loess", se = F, color = "black", formula = y ~ x) +
  geom_hline(yintercept = abs(fit_param_classII[1]) + fit_param_classII[2], linetype = "dashed", alpha = 0.6) +
  labs(x = "Number of Samples", y = "Number of unique Peptides") +
  theme(
    axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent")
  )

# Waterfall plots---------------------------------------------------------------
# TODO: plot maybe venn into watferall
# TODO: save each plot as pdf
wf_list <- list()
com_list <- list()
for (df in list(region_specific_I,region_specific_II)){
  # loop over all pairwise comparisons to get for each a waterfall plot + venn diagram
  for (com in apply(combn(names(df), 2), 2, paste, collapse = "_vs_")) {
    com_list <- append(com_list,com)
    comparison <- unlist(strsplit(com, "_vs_"))
    df_1 <- do.call(rbind.data.frame, df[comparison[1]])
    df_2 <- do.call(rbind.data.frame, df[comparison[2]])
    # make waterfall data frame
    waterfall_df <- make_waterfall_df(df_1, df_2, with_seq = TRUE)
    plot_wf <- plot_waterfall(waterfall_df,
                              comparison[1],
                              comparison[2],
                              com,
                              with_seq = TRUE
    )
    wf_list <- append(wf_list, list(plot_wf$data))
    plot(plot_wf,)
    grid::grid.newpage()
    plot_venn_waterfall(df_1$Sequence, df_2$Sequence,comparison[1],comparison[2])
  }
}
names(wf_list) <- com_list
tumor_ligands <- rbind(region_specific_I$INF,region_specific_I$T1, region_specific_I$NEC)
region_color["BEN"][[1]]
#################
comparison[1]
######
waterfall_tumor_vs_be <- make_waterfall_df(region_specific_I$T1, region_specific_I$BEN, with_seq = TRUE)
plot_waterfall_2(waterfall_tumor_vs_be, "T1", "BEN", "NEC regions vs BEN", with_seq = TRUE)
region_color["BEN"][[1]]
make_waterfall_df()

# Multiwaterfall----------------------------------------------------------------
# TODO: get whole tumor frequencies using total columns of all dfs
multi_wf_I <- make_multi_waterfall_df(region_specific_I[["NEC"]],
  region_specific_I[["T1"]],
  region_specific_I[["INF"]],
  region_specific_I[["BEN"]],
  with_seq = TRUE
)
multi_wf_II <- make_multi_waterfall_df(region_specific_II[["NEC"]],
  region_specific_II[["T1"]],
  region_specific_II[["INF"]],
  region_specific_II[["BEN"]],
  with_seq = TRUE
)
multi_wf_I_acc <- make_multi_waterfall_df(region_specific_I[["NEC"]],
  region_specific_I[["T1"]],
  region_specific_I[["INF"]],
  region_specific_I[["BEN"]],
  with_seq = FALSE
)
multi_wf_II_acc <- make_multi_waterfall_df(region_specific_II[["NEC"]],
  region_specific_II[["T1"]],
  region_specific_II[["INF"]],
  region_specific_II[["BEN"]],
  with_seq = FALSE
)

# export results for netMHCpan--------------------------------------------------
peptide_list_I <- sapply(c("NEC", "T1", "INF", "BEN", "Tumor", "Shared_all", "Tumor_shared"),function(x) NULL)
for (i in names(multi_wf_I)) {
  subset <- subset(multi_wf_I[[i]], Total > 5) # hier sollte es 5 sein wegen 37 samples insgesamt 
  peplist <- subset$uniq_seq
  new_peplist <- check_len_I(peplist)
  peptide_list_I[[i]] <- new_peplist
  write.table(new_peplist,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE,
    row.names = FALSE,
    file = paste0("/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/", i, "_peptide_list_I.tsv")
  )
}
peptide_list_II <- sapply(c("NEC", "T1", "INF", "BEN", "Tumor", "Shared_all", "Tumor_shared"),function(x) NULL)
for (i in names(multi_wf_II)) {
  subset <- subset(multi_wf_II[[i]], Total > 5)# hier sollte 5 sein weil 37 samples insgesamt 
  peplist <- subset$uniq_seq
  new_peplist <- check_len_II(peplist)
  peptide_list_II[[i]] <- new_peplist
  write.table(new_peplist,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE,
    row.names = FALSE,
    file = paste0("/Users/cschwitalla/Documents/netMHCIIpan-4.0/Tumor_regions/", i, "_peptide_list_II.tsv")
  )
}

# write patients peptides in tsv for netMHC netMHCpan-4.1/patients to predict on
# patient specific hla types

patient_peptide_list <- sapply(c(names(patient_specific_I)), function(x) NULL)
for (i in names(patient_specific_I) ){
  patient_peptide_list[[i]] <- patient_specific_I[[i]]$Sequence
  write.table(patient_specific_I[[i]]$Sequence,
              sep = "\t",
              col.names = FALSE,
              quote = FALSE,
              row.names = FALSE,
              file = paste0("/Users/cschwitalla/Documents/netMHCpan-4.1/patients/", i, "_peptide_list_I.tsv")
  )
}



################################################################################
###                              Peptide selection                           ###
################################################################################
# TODO: better refactoring still needed
# TODO: better comments

# read in file from inhouse databank--------------------------------------------
inhouse_DB_search <- read_xlsx(paste0(input_dir, "number_of_observations.xlsx"),
  col_names = TRUE
)

inhouse_DB_search$total <- rowSums(inhouse_DB_search[3:5])
inhouse_DB_search$peptide_GB_frequency_immunology_database <- (inhouse_DB_search$glioblastoma / inhouse_DB_search$total) * 100
inhouse_DB_search$peptide_benign_frequency_immunology_database <- (inhouse_DB_search$benign / inhouse_DB_search$total) * 100
inhouse_DB_search$peptide_other_cancer_frequency_immunology_database <- (inhouse_DB_search$`other cancer entities` / inhouse_DB_search$total) * 100

# HUMAN PROTEIN ATLAS-----------------------------------------------------------
# healthy tissue data set preparation
healthy_tissue_proteom <- read.csv(paste0(input_dir, "normal_tissue.tsv"),
  header = TRUE,
  sep = "\t"
)
healthy_tissue_proteom <- healthy_tissue_proteom[-which(healthy_tissue_proteom$Level == "Not detected" | healthy_tissue_proteom$Reliability == "Uncertain"), ]
healthy_tissue_proteom <- unique(healthy_tissue_proteom[c("Gene.name", "Tissue")])
# get protein abundace in healthy tissue proteom
healthy_df <- as.data.frame(table(healthy_tissue_proteom$Gene.name))
healthy_df$total <- rep.int(length(unique(healthy_tissue_proteom$Tissue)), length(healthy_df))
healthy_df$abundance <- (healthy_df$Freq / healthy_df$total) * 100

# cancer tissue data set preparation
cancer_tissue_proteom <- read.csv(paste0(input_dir, "pathology.tsv"),
  header = TRUE,
  sep = "\t"
)
cancer_tissue_proteom <- cancer_tissue_proteom[, 2:7]
cancer_tissue_proteom$sum <- rowSums(cancer_tissue_proteom[3:6])
cancer_tissue_proteom$num_detected <- rowSums(cancer_tissue_proteom[3:5])
cancer_tissue_proteom$abundance <- (cancer_tissue_proteom$num_detected / cancer_tissue_proteom$sum) * 100
cancer_protein_abundance <- cancer_tissue_proteom %>%
  dplyr::filter(cancer_tissue_proteom$Cancer != "glioma") %>%
  group_by(Gene.name) %>%
  summarise_at(vars(abundance), list(Mean = mean))

glioma_protein_abundance <- cancer_tissue_proteom %>%
  dplyr::filter(cancer_tissue_proteom$Cancer == "glioma") %>%
  group_by(Gene.name) %>%
  summarise_at(vars(abundance), list(Mean = mean))


# netMHCpan results-------------------------------------------------------------
patients_netMHC_results_I <- dir(input_dir,
                                 pattern = "ZH",
                                 full.names = TRUE)

netMHCpan_I_results <- dir(input_dir,
  recursive = TRUE,
  pattern = "_I_netMHCpan_res.txt",
  full.names = TRUE
)
netMHCpan_I_results
netMHCpan_II_results <- dir(input_dir,
  recursive = TRUE,
  pattern = "_II_netMHCpan_res.txt",
  full.names = TRUE
)
#patient dataframe
patients_netMHC_I_df <- setNames(vector("list", length = length(patients_netMHC_results_I)),
                                 c(names(patient_specific_I)))



netMHCpan_I_df <- setNames(
  vector("list", length = length(netMHCpan_I_results)),
  c("BEN","INF", "NEC","Shared_all","T1", "Tumor", "Tumor_shared")
)
netMHCpan_II_df <- setNames(
  vector("list", length = length(netMHCpan_II_results)),
  c("BEN","INF", "NEC", "Shared_all" ,"T1", "Tumor", "Tumor_shared")
)


# for patients
for (name in names(patients_netMHC_I_df)) {
  print(name)
  temp_file <- file(toString(grep(pattern = name, x = patients_netMHC_results_I, value = TRUE)))
  res_df <- netMHCpan_results_to_df_I(temp_file)
  patients_netMHC_I_df[[name]] <- res_df
}

for (name in names(netMHCpan_I_df[1:5])) {
  print(name)
  temp_file <- file(toString(grep(pattern = name, x = netMHCpan_I_results, value = TRUE)))
  res_df <- netMHCpan_results_to_df_I(temp_file)
  netMHCpan_I_df[[name]] <- res_df
}
netMHCpan_I_df[["Tumor"]] <- netMHCpan_results_to_df_I(file(toString(netMHCpan_I_results[6])))
netMHCpan_I_df[["Tumor_shared"]] <- netMHCpan_results_to_df_I(file(toString(netMHCpan_I_results[7])))

for (name in names(netMHCpan_II_df[1:5])) {
  temp_file <- file(grep(pattern = name, x = netMHCpan_II_results, value = TRUE))
  res_df <- netMHCpan_results_to_df_II(temp_file)
  netMHCpan_II_df[[name]] <- res_df
}
netMHCpan_II_df[["Tumor"]] <- netMHCpan_results_to_df_II(file(netMHCpan_II_results[6]))
netMHCpan_II_df[["Tumor_shared"]] <- netMHCpan_results_to_df_II(file(netMHCpan_II_results[7]))


for (region in names(netMHCpan_II_df)) {
  netMHCpan_II_df[region][[1]]$BindLevel <-  str_replace_all(netMHCpan_II_df[region][[1]]$BindLevel, "<=", "")
}
#-------------------------------------------------------------------------------
all_patients_netMHC_results <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c( "MHC","Peptide", "BindLevel", "Patient"))

for (name in names(patients_netMHC_I_df)) {
  new_df <- patients_netMHC_I_df[[name]][c("MHC","Peptide", "BindLevel")]
  new_df <- cbind(new_df, rep(name, times = length(new_df$MHC)))
  all_patients_netMHC_results <- rbind(all_patients_netMHC_results,new_df)
}
names(all_patients_netMHC_results)[4] <- "Patient"

all_patients_binderprediction <- summarise_binder_pred(all_patients_netMHC_results)
patients_summary <- all_patients_netMHC_results %>% group_by(Peptide) %>% summarise(toString(c(unique(Patient))))
names(patients_summary) <- c("Peptide", "Patient")
all_patients_binderprediction <- merge(all_patients_binderprediction,patients_summary[c("Peptide", "Patient")], by.x = "Peptide", by.y = "Peptide", all.x = TRUE)

#-------------------------------------------------------------------------------

# ligandom data set per patient machen + region information 
# dann nur net mhc pan mit patient specifischen hla typen 
# filter out benign und keine binder
# ligandomat --> anschauen  
#netmhcpan score 2 threshold
# sifpeiti scroe 60 

names(netMHCpan_II_df)
# create and export peptide selection dataframes for class I + II peptides
for (class in c("class_I", "class_II")) {
  if (class == "class_I") {
    peptide_list <- peptide_list_I
    binding_pred_df <- all_patients_binderprediction
    class_df <- classI_df
    multi_wf_df <- multi_wf_I
    HLA_atlas <- HLA_atlas_protein_f_cI
    peptide_class <- "_I_peptide_selection.tsv"
  } else {
    peptide_list <- peptide_list_II
    netMHCpan <- netMHCpan_II_df
    class_df <- classII_df
    multi_wf_df <- multi_wf_II
    HLA_atlas <- HLA_atlas_protein_f_cII
    peptide_class <- "_II_peptide_selection.tsv"
    binding_pred_df <- summarise_binder_pred(netMHCpan[region][[1]])
  }
  for (region in names(peptide_list)) {
  print(region)
  # use selected peptides for mergeing all other dataframes
  peptide_selection <- data.frame(Peptide = peptide_list[[region]])
  # get ppetides from netmhcpan results + binding pred +on which hla class + add region info
  #binding_pred_df <- summarise_binder_pred(netMHCpan[region][[1]])
  # mmerge peptide, acc, binder prediction,gene name, patient id, sample num
  comp_data <- class_df %>%
    group_by(Sequence) %>%
    summarise(
      Tumor_region = toString(unique(c(Tumor_region))),
      Patient_ID = toString(unique(c(Patient_ID))),
      Accessions = toString(unique(c(Accessions))),
      Gene_Name = toString(unique(c(GENENAME)))
    )
  print("1 durch")
  
  # merge patient id, acc, gene names and tumor region information
  peptide_selection <- merge(peptide_selection, 
                             comp_data,
                             by.x = "Peptide",
                             by.y = "Sequence",
                             all.x = TRUE
  )
  # merge binding prediction 
  peptide_selection <- merge(peptide_selection,
                             binding_pred_df,
                             by.x = "Peptide",
                             by.y = "Peptide",
                             all.x = TRUE
  )
  print("2 durch")
  
  # merge own data set frequncies for peptide and protein --> sind proteine wirklich da ?? 
  peptide_selection <- merge(peptide_selection,
                             multi_wf_df[region][[1]][c("uniq_seq", "Total")],
                             by.x = "Peptide",
                             by.y = "uniq_seq",
                             all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "Total"] <- "peptide_frequency - own_data"
  print("3 durch")
  
  # merge healthy tissue frequency
  peptide_selection <- merge(peptide_selection,
                             healthy_df[c("Var1", "abundance")],
                             by.x = "Gene_Name",
                             by.y = "Var1",
                             all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "abundance"] <- "healthy_tissue_protein_frequency - protein-atlas"
  print("4 durch")
  
  # merge HLA atlas protein abundaces 
  peptide_selection <- merge(peptide_selection,
                             HLA_atlas,
                             by.x = "Accessions",
                             by.y = "Accessions",
                             all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "abundance"] <- "HLA_atlas_protein_abundace"
  names(peptide_selection)[names(peptide_selection) == "percentage"] <- "HLA_atlas_protein_frequency"
  print("5 durch")
  
  
  # merge tumor tissue frequency
  peptide_selection <- merge(peptide_selection,
                             cancer_protein_abundance,
                             by.x = "Gene_Name",
                             by.y = "Gene.name",
                             all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "Mean"] <- "cancer_protein_frequency - protein-atlas"
  print("6 durch")
  
  peptide_selection <- merge(peptide_selection,
                             glioma_protein_abundance,
                             by.x = "Gene_Name",
                             by.y = "Gene.name",
                             all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "Mean"] <- "glioma_protein_frequency - protein_atlas"
  print("7 durch")
  
  # merge  inhouse frequencies
  peptide_selection <- merge(peptide_selection,
                             inhouse_DB_search[c(1:2, 7:9)],
                             by.x = "Peptide",
                             by.y = "Sequence",
                             all.x = TRUE
  )
  
  print("8 durch")
  # filter out every column that has benign in it
  
  # exclude peptides with no binding prediction 
  # peptide_selection <- peptide_selection[!is.na(peptide_selection$Weak_binders) & 
  #                                          !is.na(peptide_selection$Strong_binders) &
  #                                          is.na(peptide_selection$peptide_benign_frequency_immunology_database) |
  #                                          peptide_selection$peptide_benign_frequency_immunology_database == 0,]
  #peptide_selection <- peptide_selection[-grep("BEN", peptide_selection$Tumor_region), ]
  #peptide_selection <- peptide_selection %>% dplyr::filter(grepl(",",Patient_ID))
  # write peptide selection to output directory
  peptide_selection <- peptide_selection %>% dplyr::filter(!(is.na(Weak_binders) & is.na(Strong_binders)) & 
                                                    (is.na(peptide_benign_frequency_immunology_database) | peptide_benign_frequency_immunology_database == 0))
  
  
  write.table(peptide_selection,
              sep = "\t",
              col.names = TRUE,
              quote = FALSE,
              row.names = FALSE,
              file = paste0(output_dir, region, peptide_class)
  )
  print(paste0(output_dir, region, "_peptide_selection.tsv"))
  }
}
################################################################################
###                           Peptide selecteion 2                          ###
################################################################################

# read in peptide selection
files_peptideselection_I <- dir("/Users/cschwitalla/Documents/Ligandomics_analysis/peptide_selection_300822/",
  pattern = "_peptide_selection.tsv",
  full.names = TRUE
)
# make list of peptide selection data frames
peptide_selection_list <- list()
for (file in files_peptideselection_I) {
  df <- read.csv(file, header = TRUE, sep = "\t")
  peptide_selection_list <- append(peptide_selection_list, list(df))
}
names(peptide_selection_list) <- c(
  "BEN", "INF", "NEC", "Shared_all", "T1",
  "Tumor", "Tumor_shared"
)

peptide_selection_list$Tumor <- peptide_selection_list$Tumor[peptide_selection_list$Tumor$freq >10,]
peptide_selection_list$Shared_all <- peptide_selection_list$Shared_all[peptide_selection_list$Shared_all$freq >10,]
peptide_selection_list$Tumor_shared <- peptide_selection_list$Tumor_shared[peptide_selection_list$Tumor_shared$freq >10,]
write.table(peptide_selection_list$Shared_all,
            sep = "\t",
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE,
            file = paste0(output_dir, "Shared_all_class_I_patientid_column.tsv"))




# include peptide frequencies according to patients 
for (name in names(peptide_selection_list)) {
  new <- merge(peptide_selection_list[[name]],
    classI_df[c("Sequence", "freq")],
    by.x = "Peptide",
    by.y = "Sequence",
    all.x = TRUE
  )
  
  new$peptide_frequency...own_data <- NULL # exlcude old frequency that is on sample num 
  peptide_selection_list[[name]] <- new
}

# peptide list filtering after binding prediction 
nec_binder <- dplyr::select(dplyr::filter(peptide_selection_list$NEC, grepl(",",peptide_selection_list$NEC$Patient)),c("Peptide", "freq.x","Tumor_region"))
t1_binder <- dplyr::select(dplyr::filter(peptide_selection_list$T1, grepl(",",peptide_selection_list$T1$Patient)),c("Peptide", "freq.x","Tumor_region"))
inf_binder <- dplyr::select(dplyr::filter(peptide_selection_list$INF, grepl(",",peptide_selection_list$INF$Patient)),c("Peptide", "freq.x","Tumor_region"))
ben_binder <- dplyr::select(dplyr::filter(peptide_selection_list$BEN, grepl(",",peptide_selection_list$BEN$Patient)),c("Peptide", "freq.x","Tumor_region"))
tumor_binder <- dplyr::select(dplyr::filter(peptide_selection_list$Tumor, grepl(",",peptide_selection_list$Tumor$Patient)),c("Peptide", "freq.x","Tumor_region", "Patient"))
tumor_binder2 <- calc_peptide_selection_frequency(tumor_binder)
length(unique(tumor_binder2$Peptide))
#tumor upset 
upset_tumor <- upset_df_own2[upset_df_own2$peptides %in% c(tumor_binder$Peptide),]
upset_tumor$freq<- NULL
upset_tumor$patients <- NULL
upset_tumor <- merge(upset_tumor, tumor_binder2[,c("Peptide", "freq")], by.x = "peptides", by.y = "Peptide", all.x = TRUE)
upset_tumor <- unique(upset_tumor)

upset_tumor$freq <- round(upset_tumor$freq, digits = 2)
upset_tumor$freq <- as.character(upset_tumor$freq)
upset_tumor$patients <- NA
upset_tumor$patients <- ifelse((c(upset_tumor$freq) == "15.38"),"2 Patients",upset_tumor$patients)
upset_tumor$patients <- ifelse((c(upset_tumor$freq) == "23.08"),"3 Patients",upset_tumor$patients)
upset_tumor$patients <- ifelse((c(upset_tumor$freq) == "30.77"),"4 Patients",upset_tumor$patients)
upset_tumor$patients <- ifelse((c(upset_tumor$freq) == "38.46"),"5 Patients",upset_tumor$patients)


ComplexUpset::upset(upset_tumor, c("NEC", "T1", "INF"),
                    base_annotations = list(
                      "Intersection size" = intersection_size(
                        mapping = aes(fill = patients)
                      ) + scale_fill_manual(values = c(
                        "5 Patients" = "#E3436B",
                        "4 Patients" = "#EA5D4D",
                        "3 Patients" = "#F28140",
                        "2 Patients" = "#F49B67"
                        
                      ))),
                    queries = list(   upset_query(set = "NEC", fill = "#9e0142"),
                                      upset_query(set = "T1", fill = "#fdae61"),
                                      upset_query(set = "INF", fill = "#74add1")),
                    matrix = 
                      intersection_matrix(geom = geom_point(shape ="circle filled", size =3))+
                      scale_color_manual(
                        values = c("NEC" = "#9e0142", 
                                   "T1" = "#fdae61",
                                   "INF"= "#74add1")
                        
                      ))


#upset_tumor <- data_frame(clII_peptides)



# include expression values in TPM ---------------------------------------------
merged_genes_tpm <- read.csv("/Users/cschwitalla/Documents/transcriptomics_results/salmon.merged.gene_tpm.tsv",
                             sep = "\t",
                             header = TRUE
)
# TODO: mak it in a loop !
ben_colname <- grep("BEN", colnames(merged_genes_tpm))
nec_colname <- grep("NEC", colnames(merged_genes_tpm))
inf_colname <- grep("INF", colnames(merged_genes_tpm))
t1_colname <- grep("T1", colnames(merged_genes_tpm))
tumor_colname <- grep("T1|INF|NEC",colnames(merged_genes_tpm))

merged_genes_tpm$overall_mean_TPM <- rowMeans(merged_genes_tpm[3:51])
merged_genes_tpm$BEN_mean_TPM <- rowMeans(merged_genes_tpm[ben_colname])
merged_genes_tpm$INF_mean_TPM <- rowMeans(merged_genes_tpm[inf_colname])
merged_genes_tpm$T1_mean_TPM <- rowMeans(merged_genes_tpm[t1_colname])
merged_genes_tpm$NEC_mean_TPM <- rowMeans(merged_genes_tpm[nec_colname])
merged_genes_tpm$tumor_mean_TPM <- rowMeans(merged_genes_tpm[tumor_colname])

# merge with peptide selection 
peptide_selection_list$BEN <- merge(peptide_selection_list$BEN, merged_genes_tpm[c("gene_name", "BEN_mean_TPM")], by.x = "Gene_Name", by.y ="gene_name", all.x = TRUE)
peptide_selection_list$INF <- merge(peptide_selection_list$INF, merged_genes_tpm[c("gene_name", "INF_mean_TPM")], by.x = "Gene_Name", by.y ="gene_name", all.x = TRUE)
peptide_selection_list$NEC <- merge(peptide_selection_list$NEC, merged_genes_tpm[c("gene_name", "NEC_mean_TPM")], by.x = "Gene_Name", by.y ="gene_name", all.x = TRUE)
peptide_selection_list$T1 <- merge(peptide_selection_list$T1, merged_genes_tpm[c("gene_name", "T1_mean_TPM")], by.x = "Gene_Name", by.y ="gene_name", all.x = TRUE)
peptide_selection_list$Tumor <- merge(peptide_selection_list$Tumor, merged_genes_tpm[c("gene_name", "tumor_mean_TPM")], by.x = "Gene_Name", by.y ="gene_name", all.x = TRUE)
peptide_selection_list$Tumor_shared <- merge(peptide_selection_list$Tumor_shared, merged_genes_tpm[c("gene_name", "tumor_mean_TPM")], by.x = "Gene_Name", by.y ="gene_name", all.x = TRUE)
peptide_selection_list$Shared_all <- merge(peptide_selection_list$Shared_all, merged_genes_tpm[c("gene_name", "overall_mean_TPM")], by.x = "Gene_Name", by.y ="gene_name", all.x = TRUE)


peptide_selection_list$Tumor <- peptide_selection_list$Tumor[peptide_selection_list$Tumor$freq >10,]
peptide_selection_list$Shared_all <- peptide_selection_list$Shared_all[peptide_selection_list$Shared_all$freq >10,]
peptide_selection_list$Tumor_shared <- peptide_selection_list$Tumor_shared[peptide_selection_list$Tumor_shared$freq >10,]


# for(name in names(peptide_selection_list)) {
#   write.table(unique(peptide_selection_list[[name]]), paste0(output_dir,name, "_peptide_selection.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# }

for (name in names(peptide_selection_list)) {
  peptide_selection_list[[name]][, c("Patient_ID", "Protein.Group.Accessions", "peptide_GB_frequency_immunology_database.1")] <- NULL
  peptide_selection_list[[name]] <- unique(peptide_selection_list[[name]])
}

# CALC ALLOTYPE FREQUENCY-------------------------------------------------------
#test with inf 
test <- peptide_selection_list$INF[c("Peptide", "Patient", "Strong_binders", "Weak_binders")]
# allotype patient number 
allotype_df <- data.frame(uniqe_HLA_types)
library(data.table)
allotype_df <- melt(as.data.table(GB_HLA_types[,c(1:7)],), id.vars = "Sample Name", na.rm = TRUE)
allotype_df <- allotype_df %>% group_by(value) %>% summarise(Patient_ID = unique(`Sample Name`))
allotype_df <- allotype_df %>% group_by(value) %>% summarise(toString(c(Patient_ID)))
allotype_df$len <- stringr::str_count(allotype_df$`toString(c(Patient_ID))`,"Z")
allotype_df <- data.frame(allotype_df)

for (i in allotype_df$value) {
  allotype_df[allotype_df$value == i,1] <- paste0("HLA-",i) 
}

# include allotype frequency 
for (name in names(peptide_selection_list)) {
  print(name)
  df <- peptide_selection_list[[name]]
  df$allo_f <- NA
  for(nrow in 1:nrow(df)) {
    df[nrow, "allo_f"] <- calc_allotype_f(df[nrow,], allotype_df)
  }
  peptide_selection_list[[name]] <- df
}
# patient frequency after binding prediction
for (name in names(peptide_selection_list)){
  peptide_selection_list[[name]] <- calc_peptide_selection_frequency(peptide_selection_list[[name]])
}
# reorder columns 
column_names <- colnames(peptide_selection_list$BEN)
colnames_new <- c("Peptide","Gene_Name", "Accessions","freq.x","freq.y","allo_f","NEC_mean_TPM",
                  "Patient", "Tumor_region", "Strong_binders", "Weak_binders",
                  "HLA_atlas_protein_frequency","HLA_atlas_protein_abundace",
                  "peptide_GB_frequency_immunology_database",
                  "peptide_other_cancer_frequency_immunology_database",
                  "healthy_tissue_protein_frequency...protein.atlas",
                  "cancer_protein_frequency...protein.atlas",
                  "glioma_protein_frequency...protein_atlas")


peptide_selection_list[["BEN"]] <- peptide_selection_list[["BEN"]][colnames_new]
peptide_selection_list[["INF"]] <- peptide_selection_list[["INF"]][colnames_new]
peptide_selection_list[["T1"]] <- peptide_selection_list[["T1"]][colnames_new]
peptide_selection_list[["NEC"]] <- peptide_selection_list[["NEC"]][colnames_new]
peptide_selection_list[["Shared_all"]] <- peptide_selection_list[["Shared_all"]][colnames_new]
peptide_selection_list[["Tumor"]] <- peptide_selection_list[["Tumor"]][colnames_new]
peptide_selection_list[["Tumor_shared"]] <- peptide_selection_list[["Tumor_shared"]][colnames_new]






# write peptide selections to tsv
for (name in names(peptide_selection_list)){
  
  write.table(peptide_selection_list[[name]],
              sep = "\t",
              col.names = TRUE,
              quote = FALSE,
              row.names = FALSE,
              file = paste0(output_dir, "/peptide_selection200922/",name, "_class_I.tsv"))
  
}




# calc_allotype_f <- function(x, allotype_df) {
#   if (!is.na(x$Strong_binders)) {  # use strong binders as default
#     if (!grepl(",", x$Strong_binders)) { # if mutliple strong binders then skip
#       allotype <- x$Strong_binders
#       #check which patients have allotypes
#       pattern <- paste(as.list(strsplit(allotype_df[allotype_df$value == allotype,]$toString.c.Patient_ID.., ", " )[[1]]), collapse = "|")
#       test_patients <- c(strsplit(x$Patient,", ")[[1]])
#       patients <- grep(pattern, test_patients, value = TRUE )
#       n_patiens <- stringr::str_count(paste(patients, collapse = ","), "Z")
#       totl_patiens <- allotype_df[allotype_df$value == allotype, ]$len
#       allo_f <- n_patiens / totl_patiens * 100
#       print(allo_f)
#     } else {
#       return(NA)
#     }
#   } else {
#     if (!grepl(",", x$Weak_binders, fixed = TRUE)) { # if mutliple strong binders then skip
#       allotype <- x$Weak_binders
#       print("weak_b")
#       pattern <- paste(as.list(strsplit(allotype_df[allotype_df$value == allotype,]$toString.c.Patient_ID.., ", " )[[1]]), collapse = "|")
#       test_patients <- c(strsplit(x$Patient,", ")[[1]])
#       patients <- grep(pattern, test_patients, value = TRUE )
#       n_patiens <- stringr::str_count(paste(patients, collapse = ","), "Z")
#       totl_patiens <- allotype_df[allotype_df$value == allotype, ]$len
#       allo_f <- n_patiens / totl_patiens * 100
#       return(allo_f)
#     } else {
#       return(NA)
#     }
#   }
# }


################################################################################
###                           Neoepitopes                                    ###
################################################################################
filter_inhouse_db <- subset(inhouse_DB_search, benign != 0)

# load neoepitope data 
neo_epitopes <- read.csv(paste0(input_dir, "neoepitopes_I.tsv"),
                                         header = TRUE,
                                         sep = "\t"
)
# filter neoepitopes with benign db  -------------------------------------------
neo_epitopes_filtered <- subset(neo_epitopes, !(sequence %in% combi_benign_pep_II |sequence %in% combi_benign_pep_I | sequence %in% filter_inhouse_db$Sequence))
num_genes <- length(unique(neo_epitopes_filtered$gene))
num_peptides <- length(unique(neo_epitopes_filtered$sequence))
num_proteins <- length(unique(neo_epitopes_filtered$proteins))

# calc peptide frequency 


neo_epitopes_filtered <- calc_neoepitopes_frequency(neo_epitopes_filtered)
neo_sum_patients <- neo_epitopes_filtered %>% dplyr::group_by(sequence) %>% dplyr::summarise(patient_id = unique(patient_id))
neo_sum_patients <- neo_sum_patients %>% dplyr::group_by(sequence) %>% dplyr::summarise(patient_id = toString(c(patient_id)))
neo_sum_regions <- neo_epitopes_filtered %>% dplyr::group_by(sequence) %>% dplyr::summarise(region = unique(region))
neo_sum_regions <- neo_sum_regions %>% dplyr::group_by(sequence) %>% dplyr::summarise(region = toString(c(region)))
neo_epitopes_filtered_final <- merge(neo_sum_patients,neo_sum_regions, by = "sequence", all.x = TRUE)
neo_epitopes_filtered_final <- merge(neo_epitopes_filtered_final,unique(neo_epitopes_filtered[c("sequence","GENENAME.y", "freq")]), by = "sequence", all.x = TRUE)
neo_epitopes_filtered_final <- unique(neo_epitopes_filtered_final)

#--replace na
na_epitipe <- neo_epitopes_filtered_final[is.na(neo_epitopes_filtered_final$GENENAME.y),]
merge_neoepitope <- neo_epitopes_filtered[,c("sequence", "GENENAME.x")]
new <- merge(na_epitipe,merge_neoepitope,by = "sequence", all.x = TRUE)
new$GENENAME.y <- NULL
names(new)[names(new) == "GENENAME.x"] <- "GENENAME.y"
new <- new[,c(1,2,3,5,4)]
neo_epitopes_filtered_final <- neo_epitopes_filtered_final[!(neo_epitopes_filtered_final$sequence %in% na_epitipe$sequence),]
neo_epitopes_filtered_final <- rbind(neo_epitopes_filtered_final, new)
neo_epitopes_filtered_final <- unique(neo_epitopes_filtered_final)
# neoepitopes that occure in at least 2 patients
frequent_neoepitopes <- neo_epitopes_filtered_final[neo_epitopes_filtered_final$freq >10,]
# epitopes only occuring in tumor regions NEC, T1, INF
tumor_neoepitopes <- neo_epitopes_filtered_final[-grep("BEN", neo_epitopes_filtered_final$region),]

#length(unique(neo_epitopes_filtered_final$sequence))
split_epitopes <- split(neo_epitopes_filtered_final, neo_epitopes_filtered_final$region)

region_epitopes <- split(neo_epitopes_filtered, neo_epitopes_filtered$region)  
set_region_distrib <- list(NEC = region_epitopes$NEC$sequence,
                           T1 = region_epitopes$T1$sequence,
                           INF = region_epitopes$INF$sequence,
                           BEN = region_epitopes$BEN$sequence)
neo1 <- plot_custom_venn(set_region_distrib, "")



set_region_distrib_genes <- list(NEC = region_epitopes$NEC$GENENAME.x,
                           T1 = region_epitopes$T1$GENENAME.x,
                           INF = region_epitopes$INF$GENENAME.x,
                           BEN = region_epitopes$BEN$GENENAME.x)
neo2 <- plot_custom_venn(set_region_distrib_genes, "")

grid.arrange(neo1, neo2,
             ncol = 2)
################################################################################
###                           Integration                                    ###
################################################################################
# #venn for testing the filtering ------------------------------------------------
# set_neoepitopes <- list(benign = combi_benign_pep_I, own_data = classI_df$Sequence, neoepitopes = neo_epitopes_filterd_benign$sequence, class2 = classII_df$Sequence)
# 
# set_neo_prot <- list(class_I = classI_df$Accessions, neoepitopes = neo_epitopes_filterd_benign$UNIPROTID, classs_II = classII_df$Accessions)
# set_neo_prot2 <- list(class_I = classI_df$Accessions, neoepitopes = neo_epitopes_filterd_benign$UNIPROTID)
# 
# #plotting venn 
# venn <- plot_custom_venn(set_neoepitopes, "filtering")
# venn2 <- plot_custom_venn(set_neo_prot2, "immunopeptidome vs neoepitopes")
# 
# neo_epitopes_filterd_benign <- unique(neo_epitopes_filterd_benign)
# 
# #region_neoepitopes <- split(neo_epitopes_filterd,neo_epitopes_filterd_benign$region)
# set_neo_regions <- list(NEC = region_neoepitopes$NEC$sequence, T1 = region_neoepitopes$T1$sequence, INF =region_neoepitopes$INF$sequence , BEN = region_neoepitopes$BEN$sequence)
# venn_neo_regions <- plot_custom_venn(set_neo_regions, "region distribution neoepitopes")

#-------------------------------------------------------------------------------


# venn_INF <-  plot_custom_venn(list(neoepitopes = region_neoepitopes$INF$UNIPROTID,
#                                    INF_excl = peptide_selection_list$INF$Accessions,
#                                    Tumor = peptide_selection_list$Tumor$Accessions,
#                                    ALL = classI_df$Accessions),
#                               "inf neoepitopes")
# venn_INF@region$item
# test <- classI_df %>% filter(classI_df$Accessions == "Q9H3M7")
# test <- data.frame(HLA_atlas_data) %>% filter(HLA_atlas_data$`HLA-I`$acc == "Q9H3M7")


# VENN INTEGRATION--------------------------------------------------------------

# 1: comparing tumor from each data set 
set_tumor <- list(exome = unique(exome_data$Tumor$Hugo_Symbol),
                  peptidome = peptide_selection_list$Tumor$Gene_Name ,
                  neoepitopes = tumor_neoepitopes$GENENAME)
plot_custom_venn(set_tumor, "tumor")
unique(neo_epitopes$GENENAME)
# 2: comparing all from all 
set_all <- list(exome = exome_data$ALL$Hugo_Symbol,
                peptidome = classI_df$GENENAME,
                neoepitopes = neo_epitopes_filtered_final$GENENAME )
set_all_degenes <- list(exome = exome_data$ALL$Hugo_Symbol,
                peptidome = classI_df$GENENAME,
                neoepitopes = neo_epitopes_filtered_final$GENENAME,
                de_genes = de_gene_list$V1)
integration_all <- plot_custom_venn(set_all, "all")
integration_all_degenes <- plot_custom_venn(set_all_degenes, "all + de genes ")
# 3: compraing high frequent from all 
set_high_f <- list(exome = exome_high_f$Hugo_Symbol,
                   peptidome = classI_df[classI_df$freq >15,]$GENENAME,
                   neoepitopes = neoepitopes_final[neoepitopes_final$freq > 13,]$GENENAME)

plot_custom_venn(set_high_f, "High frequency")

# import de genes
small_tpm <- merged_genes_tpm[c(1,2,52,53,54,55,56,57)]

de_gene_list <- read.table("/Users/cschwitalla/Documents/RNAseq_analysis/DE_gene_list.tsv")

#import data from ca atlas
ca_antigens <- read.csv("/Users/cschwitalla/Documents/Immunopeptidomics/CA_antigens_ca_atlas.tsv",
                          sep = "\t")
ct_antigens <- read.csv("/Users/cschwitalla/Documents/Immunopeptidomics/CT_antigens_ca_atlas.tsv",
                          sep = "\t")
ptm_antigens <- read.csv("/Users/cschwitalla/Documents/Immunopeptidomics/PTM_antigens_ca_atalas.tsv",
                           sep = "\t")




set_all__ca<- list(exome = exome_data$ALL$Hugo_Symbol,
                        peptidome = classI_df$GENENAME,
                        neoepitopes = neo_epitopes_filtered_final$GENENAME,
                        de_genes = de_gene_list$V1,
                    ca_antigens = ca_antigens$GeneName)

integration_all_ct <- plot_custom_venn(set_all__ca , "all + ca")


#upset plot 
library("UpSetR")

set_all_all <- list(exome = exome_data$ALL$Hugo_Symbol,
                    immunopeptidome = classI_df$GENENAME,
                    neoepitopes = neo_epitopes_filtered_final$GENENAME,
                    de_genes = de_gene_list$V1,
                   ca_antigens = ca_antigens$GeneName,
                    ct_antigens = ct_antigens$GeneName,
                    ptm_antigens = ptm_antigens$GeneName)

length(set_all_all)
UpSetR::upset(fromList(set_all_all), order.by = "freq", empty.intersections = "on", nsets = length(set_all_all) )

# list to matrix

# frequently observed peptides in tumor regions NEC, T1, INF 
# in recurrent exome mutation genes 
#in neoepitopes 
# set_petidome_neo_exomehighf <- list(neo = neoepitopes_final$GENENAME,
#                                peptidome = peptide_selection_list$Tumor$Gene_Name,
#                                exome = exome_high_f_$Hugo_Symbol)
# 
# venn_petidome_neo_exome <- plot_custom_venn(set_petidome_neo_exomehighf, "exome_high_f, peptidome, neo")
# 
# 
# 
# set_petidome_neo_exomehighf_tumor <- list(neo = neoepitopes_final$GENENAME,
#                                     peptidome = peptide_selection_list$Tumor$Gene_Name,
#                                     exome = exome_high_f_tumor$Hugo_Symbol)
# venn_petidome_neo_exome <- plot_custom_venn(set_petidome_neo_exomehighf_tumor, "exome_high_f_tumor, peptidome, neo")
# 
# 

#--> was ist die tpm 

# DE_results <- list()
# for ( file in files_rna) {
#   df <- read.csv(file, header = TRUE, sep = "\t")
#   DE_results <- append(DE_results, list(df))
# }
# names(DE_results) <- c("BEN_vs_INF", "BEN_vs_NEC", "BEN_vs_T1", "INF_vs_NEC", "INF_vs_T1", "NEC_vs_T1")
# 
# mapping_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
#                                     keys = c(DE_results$BEN_vs_INF$Names,
#                                              DE_results$BEN_vs_NEC$Names,
#                                              DE_results$BEN_vs_T1$Names,
#                                              DE_results$INF_vs_NEC$Names,
#                                              DE_results$INF_vs_T1$Names,
#                                              DE_results$NEC_vs_T1$Names),
#                                     columns = "UNIPROTID",
#                                     keytype = "GENENAME"
# )

# for (i in names(DE_results)) {
#   new_df <- merge(DE_results[[i]], mapping_df, by.x = "Names", by.y = "GENENAME", all.x = TRUE )
#   DE_results[[i]] <- new_df
# }
# de_genes <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
#                                   keys = DE_genes,
#                                   columns = "UNIPROTID",
#                                   keytype = "GENENAME")

#--EXOME DATA-------------------------------------------------------------------
files_exome <- dir("/Users/cschwitalla/Documents/WES_analysis/",
                   pattern = "_variants.tsv",
                   full.names = TRUE)


exome_data <- list()
for (file in files_exome) {
  df <- read.csv(file, header = TRUE, sep = ",")
  exome_data <- append(exome_data, list(df))
}
names(exome_data) <- c("BEN", "INF", "NEC", "T1")
# include region info as column for each 
exome_data$BEN$region <- "BEN"
exome_data$NEC$region <- "NEC"
exome_data$T1$region <- "T1"
exome_data$INF$region <- "INF"


#mapp proteinids to acc
createprotein2acc <- function(){
  edb <- EnsDb.Hsapiens.v86
  k <- keys(edb, keytype = "PROTEINID")
  gene2gene <- AnnotationDbi::select(edb, k, "UNIPROTID", "PROTEINID")
  return(gene2gene)
  
}
protein2acc <- createprotein2acc()


all_variants <- c()
for (i in names(exome_data)) {
  new_df <- merge(exome_data[[i]], protein2acc, by.x = "ENSP", by.y = "PROTEINID", all.x = TRUE )
  new_df <- merge(new_df, exome_metadata[c("Tumor_Sample_Barcode", "Patient_ID")], by = "Tumor_Sample_Barcode", all.x = TRUE)
  exome_data[[i]] <- new_df
  all_variants <- append(all_variants,new_df$UNIPROTID)
}
# get dataframe with all region variants
exome_data$ALL <- rbind(exome_data$INF, exome_data$T1, exome_data$NEC, exome_data$BEN)
#get dataset with only tumor region NEC T1 INF
exome_data$Tumor <- rbind(exome_data$INF, exome_data$T1, exome_data$NEC)
# exclude all genes that are present in benign 
exome_data$Tumor <- exome_data$Tumor[!(exome_data$Tumor$Hugo_Symbol %in% c(exome_data$BEN$Hugo_Symbol)),]


# get frequency of genes mutated in patients
exome_frequency <- exome_data$Tumor[c("Patient_ID", "Hugo_Symbol")]
exome_frequency <- exome_frequency %>% dplyr::group_by(Hugo_Symbol) %>% dplyr::summarise(Patient_ID =unique(Patient_ID))
exome_frequency <- exome_frequency %>% dplyr::group_by(Hugo_Symbol) %>% dplyr::summarise(toString(c(Patient_ID)))
exome_frequency$len <- stringr::str_count(exome_frequency$`toString(c(Patient_ID))`,"Z")
exome_frequency$freq <- exome_frequency$len / length(unique(exome_data$ALL$Patient_ID)) *100

# summarise 

for (i in names(exome_data)) {
  new_df <- merge(exome_data[[i]], exome_frequency, bx = "Hugo_Symbol", all.x = TRUE)
  new_df[c("Tumor_Sample_Barcode", "Patient_ID")] <-  NULL
  exome_data[[i]] <- unique(new_df)
}

for(i in names(exome_data)) {
  new_df <- exome_data[[i]] %>%
    group_by(Hugo_Symbol) %>%
    summarise(frequency = freq, patients = `toString(c(Patient_ID))`, region = toString(c(unique(region)))) %>%
    unique()
  exome_data[[i]] <- unique(new_df)
}



# exome_high_f <- exome_data$ALL[exome_data$ALL$freq > 15,]
# exome_high_f <- exome_high_f %>% group_by(Hugo_Symbol) %>% summarise(frequency = freq, patients = `toString(c(Patient_ID))`, region = toString(c(unique(region)))) %>% unique()
# exome_high_f <- exome_high_f %>% group_by(Hugo_Symbol) %>% summarise(toString(c(region)))
# exome_high_f_tumor <- exome_data$Tumor[exome_data$Tumor$freq > 15,]
# exome_high_f_tumor <- exome_high_f_tumor %>% group_by(Hugo_Symbol) %>% summarise(frequency = freq, patients = `toString(c(Patient_ID))`, region = toString(c(unique(region)))) %>% unique()
# exome_high_f_tumor <- exome_high_f_tumor %>% group_by(Hugo_Symbol) %>% summarise(toString(c(region)))




# ven_exom_immunop <- list( exome_frequent = exome_frequency$Hugo_Symbol, immunopeptidome = tumor_ligands$GENENAME, neo_epitopes = neoepitopes_final$GENENAME.x )
# exom_vs_immunopeptidome <- plot_custom_venn(ven_exom_immunop, "exom vs immunopeptidome")
# 
# 
# selection <- c(exom_vs_immunopeptidome@region$item[[7]])
# selection_df <- exome_high_f[exome_high_f$Hugo_Symbol %in% selection,]

#DE results---------------------------------------------------------------------
de_genes <- read.csv("/Users/cschwitalla/Documents/RNAseq_analysis/DE_gene_list.tsv", header = FALSE)
#-------------------------------------------------------------------------------
exome_frequent_mutations <- exome_frequency[exome_frequency$len>=2,]
#exome_frequent_mutations_xBEN <- exome_frequency[exome_frequency$]
new_tumor_binder <- merge(tumor_binder2, peptide_selection_list$Tumor[,c("Peptide", "Gene_Name")], by = "Peptide", all.x = TRUE )
new_tumor_binder <- unique(new_tumor_binder)

neo_epitopes_filtered2 <- unique(neo_epitopes_filtered)
venn_set_all <- list(
  DE_genes = de_genes$V1,
  exome = exome_data$ALL$Hugo_Symbol,
  ligandome = peptide_selection_list$Tumor$Gene_Name,
  neo_epitope = neo_epitopes_filtered_final$GENENAME.y
)
all <- plot_custom_venn(venn_set_all, "all")

# venn with frequent exome mutations
venn_frequent_exome <- list(
 "DE genes" = de_genes$V1,
  "somatic mutations" = exome_frequent_mutations$Hugo_Symbol,
  "immunopeptidome" = new_tumor_binder$Gene_Name,
  "neo-epitope" = neo_epitopes_filtered_final$GENENAME.y
)
freq_exome_venndata <- plot_custom_venn(venn_frequent_exome, "")


neo_epitope_ligandome <- list(
  "neo-epitopes" = neo_epitopes_filtered_final$sequence,
  "immunopeptidome" = tumor_binder2$Peptide
)
neo_vs_binder <- plot_custom_venn(neo_epitope_ligandome,"")




acc2gene <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                  keys = c("Q460N5" ,"Q9NQT8" ,"Q9H1A4"),
                                  columns = "GENENAME",
                                  keytype = "UNIPROTID")
all@region$item


#-------------------------------------------------------------------------------
# ven1 <- plot_custom_venn(list(neo_epitopes = neo_epitopes$GENENAME.y, exome_INF = exome_data$INF, DE_genes = DE_results$INF_vs_T1$Names), "ven1")
# 
# test <- ven1@region$item[[7]]
# ( c("TXNIP", "DNAH8") %in% test)
# 
# inf_vs_nec <- data.frame(DE_results$INF_vs_NEC) %>% filter(DE_results$INF_vs_NEC$Names %in% c("TXNIP", "DNAH8"))
# inf_vs_t1 <- data.frame(DE_results$INF_vs_T1) %>% filter(DE_results$INF_vs_T1$Names %in% c("TXNIP", "DNAH8"))
# ben_vs_inf <- data.frame(DE_results$BEN_vs_INF) %>% filter(DE_results$BEN_vs_INF$Names %in% c("TXNIP", "DNAH8"))
# test2 <- data.frame(counts_data) %>% filter(row.names(data.frame(counts_data)) %in% c("TXNIP", "DNAH8"))
# names(test2) <- dds_default@colData$Tumor_region
# test2 <- t(test2)
# 
# test4 <- neo_epitopes %>% filter(neo_epitopes$GENENAME.y %in% c("TXNIP", "DNAH8"))
# 
