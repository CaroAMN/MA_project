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
# import functions for the analysis
source("/Users/cschwitalla/git/students/cschwitalla/Ligandomics_Analysis/functions_ligandomics.R")
################################################################################
###                             Load libraries                               ###
################################################################################
required_libs <- c(
  "tidyr", "readxl", "ggVennDiagram", "dplyr", "stringr",
  "tibble", "ggplot2", "org.Hs.eg.db", "RColorBrewer", "EnsDb.Hsapiens.v86"
)

suppressMessages(invisible(lapply(required_Libs, library, character.only = TRUE)))


################################################################################
###                            Load Data                                     ###
################################################################################
# load HLA-typing dataframe form Marcel Immunology------------------------------
GB_HLA_types <- read_xlsx(paste0(input_dir, "HLA-Typisierung_GBM.xlsx"),
  col_names = TRUE
)

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


# HLA-Ligand Atlas data ---------------------------------------------------------
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
  group_by(peptide_sequence_id) %>%
  summarise(acc = toString(uniprot_id))

# merge dataframes by peptide sequence id to add protein accession numbers to the dataframe
HLA_atlas_data <- merge(HLA_ligand_atlas_pep, HLA_ligand_atlas_acc, by = "peptide_sequence_id")
# split the ligand atlas according to the hla class ( 1 or 2 )
HLA_atlas_data <- split(HLA_atlas_data, HLA_atlas_data$hla_class)

# TUMOR REGION LIGANDOMICS DATA-------------------------------------------------
# dataframe mappin uniprot ids to gene names
mapping_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
  keys = accesionList,
  columns = "GENENAME",
  keytype = "UNIPROTID"
)


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

# read in data
classI_df <- create_data_frame(files_cI)
classII_df <- create_data_frame(files_cII)

# DATA PREPERATION -------------------------------------------------------------
# combine benign data

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
# filter my datasets so that the known benign peptides frombenign pep db are not
# there anymore
classI_df <- subset(classI_df, !(Sequence %in% combi_benign_pep_I))
classII_df <- subset(classII_df, !(Sequence %in% combi_benign_pep_II))

# dataset that hase only peptides predicted to originate only from one protein
# not multiple
classI_df <- classI_df[-grep(";", classI_df$Accessions), ]
classI_df$Accessions <- get_protein_acc(classI_df$Accessions)

classII_df <- classII_df[-grep(";", classII_df$Accessions), ]
classII_df$Accessions <- get_protein_acc(classII_df$Accessions)

classI_df <- merge(classI_df, mapping_df, by = "Accessions", all.x = TRUE)
classII_df <- merge(classII_df, mapping_df, by = "Accessions", all.x = TRUE)
# create than dataframes that are filterd+uniqe origin that are region specific
region_specific_I <- split(classI_df, classI_df$Tumor_region)
region_specific_II <- split(classII_df, classII_df$Tumor_region)

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


# Length distribution ----------------------------------------------------------
# TODO: adjust color for tumor regions
# TODO: adjust output that plot is saved as pdf in output dir


plot_length_distribution(classI_df,
  as_bar = TRUE,
  "length distributeion of class I HLA peptides"
)


# saturation analysis ===================
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
# TODO: adjust color for tumor regions
# TODO: save each plot as pdf

# for loop mit allen combination wie in rna seq analysis

for (com in apply(combn(names(region_specific_I), 2), 2, paste, collapse = "_vs_")) {
  comparison <- unlist(strsplit(com, "_vs_"))
  df_1 <- do.call(rbind.data.frame, region_specific_I[comparison[1]])
  df_2 <- do.call(rbind.data.frame, region_specific_I[comparison[2]])
  # make waterfall data frame
  waterfall_df <- make_waterfall_df(df_1, df_2, with_seq = FALSE)
  plot_wf <- plot_waterfall(waterfall_df,
    comparison[1],
    comparison[2],
    com,
    with_seq = FALSE
  )
  plot(plot_wf)
  grid::grid.newpage()
  plot_venn_waterfall(df_1$Accessions, df_2$Accessions)
}


# Multiwaterfall----------------------------------------------------------------
# TODO: get whole tumor frequencies using total columns of all dfs
multi_wf_I <- make_multi_waterfall_df(region_specific_I[["NEC"]],
  region_specific_I[["T1"]],
  region_specific_I[["INF"]],
  with_seq = TRUE
)
multi_wf_II <- make_multi_waterfall_df(region_specific_II[["NEC"]],
  region_specific_II[["T1"]],
  region_specific_II[["INF"]],
  with_seq = TRUE
)
multi_wf_I_acc <- make_multi_waterfall_df(region_specific_I[["NEC"]],
  region_specific_I[["T1"]],
  region_specific_I[["INF"]],
  with_seq = FALSE
)
multi_wf_II_acc <- make_multi_waterfall_df(region_specific_II[["NEC"]],
  region_specific_II[["T1"]],
  region_specific_II[["INF"]],
  with_seq = FALSE
)


# export results for netMHCpan--------------------------------------------------

for (i in names(multi_wf_I)) {
  subset <- subset(multi_wf_I[[i]], Total > 4)
  peplist <- subset$uniq_seq
  new_peplist <- check_len_I(peplist)
  write.table(new_peplist,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE,
    row.names = FALSE,
    file = paste0("/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/", i, "_peptide_list_I.tsv")
  )
}
for (i in names(multi_wf_II)) {
  subset <- subset(multi_wf_II[[i]], Total > 4)
  peplist <- subset$uniq_seq
  new_peplist <- check_len_II(peplist)
  write.table(new_peplist,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE,
    row.names = FALSE,
    file = paste0("/Users/cschwitalla/Documents/netMHCIIpan-4.0/Tumor_regions/", i, "_peptide_list_II.tsv")
  )
}


################################################################################
###                             Peptide selection                           ###
################################################################################
# TODO: better refactoring still needed
# TODO: better comments

# read in file from inhouse databank--------------------------------------------
inhouse_DB_search <- read_xlsx(paste0(input_dir, "number_of_observations.xlsx"),
  col_names = TRUE
)

inhouse_DB_search$total <- rowSums(inhouse_DB_search[3:5])
inhouse_DB_search$pep_GB_f <- (inhouse_DB_search$glioblastoma / inhouse_DB_search$total) * 100
inhouse_DB_search$pep_benign_f <- (inhouse_DB_search$benign / inhouse_DB_search$total) * 100
inhouse_DB_search$pep_other_cancer_f <- (inhouse_DB_search$`other cancer entities` / inhouse_DB_search$total) * 100

# HUMAN PROTEIN ATLAS-----------------------------------------------------------
# healthy tissue data set preparation
healthy_tissue_proteom <- read.csv(paste0(input_dir, "normal_tissue.tsv"),
  header = TRUE,
  sep = "\t"
)
healthy_tissue_proteom <- healthy_tissue_proteom[-which(healthy_tissue_proteom$Level == "Not detected" | healthy_tissue_proteom$Reliability == "Uncertain"), ]
# get protein abundace in healthy tissue proteom
healthy_df <- as.data.frame(table(healthy_tissue_proteom_f$Gene.name))
healthy_df$total <- rep.int(length(unique(healthy_tissue_proteom_f$Tissue)), length(healthy_df))
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
netMHCpan_I_results <- dir(input_dir,
  recursive = TRUE,
  pattern = "_I_netMHCpan_res.txt",
  full.names = TRUE
)

netMHCpan_II_results <- dir(input_dir,
  recursive = TRUE,
  pattern = "_II_netMHCpan_res.txt",
  full.names = TRUE
)


netMHCpan_I_df <- setNames(
  vector("list", length = length(netMHCpan_I_results)),
  c("INF", "Nec", "T1", "Tumor")
)
netMHCpan_II_df <- setNames(
  vector("list", length = length(netMHCpan_II_results)),
  c("INF", "Nec", "T1", "Tumor")
)

for (name in names(netMHCpan_I_df)) {
  temp_file <- file(grep(pattern = name, x = netMHCpan_I_results, value = TRUE))
  res_df <- netMHCpan_results_to_df_I(temp_file)
  netMHCpan_I_df[[name]] <- res_df
}
for (name in names(netMHCpan_II_df)) {
  temp_file <- file(grep(pattern = name, x = netMHCpan_II_results, value = TRUE))
  res_df <- netMHCpan_results_to_df_II(temp_file)
  netMHCpan_II_df[[name]] <- res_df
}

#-------------------------------------------------------------------------------


################################################################################
# loop over all region in net MHCp results
for (region in names(netMHCpan_I_df)) {
  # get ppetides from netmhcpan results + binding pred +on which hla class + add region info
  binding_pred_df <- summarise_binder_pred(netMHCpan_I_df[region][[1]])
  # mmerge peptide, acc, binder prediction,gene name, patient id, sample num
  comp_data <- classI_df %>%
    group_by(Sequence) %>%
    summarise(
      Tumor_region = toString(unique(c(Tumor_region))),
      Patient_ID = toString(unique(c(Patient_ID))),
      Accessions = toString(unique(c(Accessions))),
      Gene_Name = toString(unique(c(GeneName)))
    )

  peptide_selection <- merge(binding_pred_df,
    comp_data,
    by.x = "Peptide",
    by.y = "Sequence",
    all.x = TRUE
  )

  # merge own data set frequncies for peptide and protein
  peptide_selection <- merge(peptide_selection,
    multi_wf_I[region][[1]][c("uniq_seq", "Total")],
    by.x = "Peptide",
    by.y = "uniq_seq",
    all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "Total"] <- "peptide_f_own_data"

  # merge healthy tissue frequency
  peptide_selection <- merge(peptide_selection,
    healthy_df[c("Var1", "abundance")],
    by.x = "Gene_Name",
    by.y = "Var1",
    all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "abundance"] <- "healthy_tissue_protein_f"

  # merge tumor tissue frequency
  peptide_selection <- merge(peptide_selection,
    cancer_protein_abundance,
    by.x = "Gene_Name",
    by.y = "Gene.name",
    all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "Mean"] <- "cancer_protein_f"

  peptide_selection <- merge(peptide_selection,
    glioma_protein_abundance,
    by.x = "Gene_Name",
    by.y = "Gene.name",
    all.x = TRUE
  )
  names(peptide_selection)[names(peptide_selection) == "Mean"] <- "glioma_protein_f"

  # merge  inhouse frequencies
  peptide_selection <- merge(peptide_selection,
    inhouse_DB_search[c(1:2, 7:9)],
    by.x = "Peptide",
    by.y = "Sequence",
    all.x = TRUE
  )

  # filter out every column that has benign in it
  peptide_selection <- peptide_selection[-grep("BEN",peptide_selection$Tumor_region), ]

  # write peptide selection to output directory
  write.table(peptide_selection,
    sep = "\t",
    col.names = TRUE,
    quote = FALSE,
    row.names = FALSE,
    file = paste0(output_dir, region, "_peptide_selection.tsv")
  )
}
