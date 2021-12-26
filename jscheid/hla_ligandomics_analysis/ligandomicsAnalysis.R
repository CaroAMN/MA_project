# ------------------------
#         PURPOSE
# ------------------------
# Ligandomic analysis includes:
# - Class1 and Class2 peptide yields for all samples
# - Length distribution of Class1 and Class2 peptides
# - Waterfallplot between primary + recurrent VS benign data base
# - Waterfallplot between primary VS recurrent (in general and in autologous samples)
# - Saturation analysis
# - Vulcano plot over all patients
# - Correlation mass VS no. peptides
# - Population coverage
# - Peptide-derived protein coverage
# - Binder-non binder visualization
# - FDR calculation of (tumor) exclusive peptides
# - Comparison of cancer testis antigene derived from proteins ⇒ CTDatabase
# - Crosspresentation

# ----------------------------------
#       Manual preparation
# ----------------------------------
# - For Class1 peptides: Predict according to their HLA-type using NetMHCpan 4.0 and SYPFEITHI
#                        Filter NetMHCpan 4.0 Rank < 2 or SYFPEITHI Score > 0.6/60%
#                        Those peptides are considered binders
# - For Class2 peptides: Take every detected peptide and once interesting targets occur. 
#                        Predict only those, which are relevant after the analysis using NetMHCIIpan 4
# - The excel tables have been transformed into tsv files for easier downstream analysis
# - The purity has been manually extracted from 20200204_Sample Overview_GBM_P-R.xlsx into a tsv file
# 

# ----------------------------------
#       Mandatory libraries
# ----------------------------------
setwd("~/Documents/Uni/Master/Thesis/Ligandomics")

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library("RColorBrewer")
library(foreach)
library(tuple)
library(doParallel)
library(readxl)
library(aomisc)
library(VennDiagram)

# ----------------------------------
#       Read in necessary files
# ----------------------------------
binders_classI <- read.csv("Class1_binders.tsv", sep='\t')[-1]
peptides_classII <- read.csv("Class2_peptides.tsv", sep='\t')[-1]
purity_classI <- read.csv("Class1_purity.tsv", sep='\t')
benign_DB_classI <- read_xlsx('Class1_benign_DB.xlsx', sheet = 1)
benign_DB_classII <- read_xlsx('Class2_benign_DB.xlsx', sheet = 1)
# ----------------------------------
#     Peptide-Sample yields
# ---------------------------------------------------------------------------------------------------
# This part aims to visualize the number of measured peptides as well as the purity calculated by: 
# binders / (binders + non-binders)
# ---------------------------------------------------------------------------------------------------
# CLASS 1
# ---------
# Prepare class1 data for plotting
yields_classI <- binders_classI %>%
      group_by(Sample_ID) %>%
      summarize(n=n()) %>%
      mutate(gbm_status = ifelse(grepl("prim", Sample_ID), "P", "R")) %>%
      mutate(sample_num = gsub("GBM(.*)_.*","\\1", Sample_ID)) %>%
      arrange(match(sample_num, sample_num[order(strtoi(sample_num))])) %>%
      mutate(Sample_ID = factor(Sample_ID, levels = unique(Sample_ID)))
# Rich barplot of class1 peptide yields
yields_classI %>% 
  ggplot(., aes(x = Sample_ID, y = n)) +
  geom_bar(stat="identity", fill = ifelse(yields_classI$gbm_status == "P", brewer.pal(3, 'RdBu')[3], brewer.pal(3, 'RdBu')[1])) +
  geom_point(aes(x = purity_classI$Sample_ID, y = purity_classI$Purity * 80), shape = 18, size = 3, alpha = 0.7) +
  geom_hline(yintercept = mean(yields_classI$n), linetype= "dashed", alpha = 0.5) +
  annotate(geom = "text", x = seq_len(nrow(yields_classI)), y = -300, label = yields_classI$gbm_status, size = 4) + 
  annotate(geom = "text", x = 2 * (1:(nrow(yields_classI)/2)) - 0.5, y = -700, label = unique(yields_classI$sample_num), size = 4) + 
  annotate(geom = "text", x = -1, y = -700, label = "GBM", size = 4) + 
  coord_cartesian(ylim=c(0,8000), clip="off") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),
                     sec.axis = sec_axis(~./80, name = "Purity [%]")) +
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("Number of binders")

# ---------
# CLASS 2
# ---------
# Since class 2 peptides cannot be predicted at the moment, we simply visualize the number of measured peptides
# Prepare class2 data for plotting
yields_classII <- peptides_classII %>%
  group_by(Sample_ID) %>%
  summarize(n=n()) %>%
  mutate(gbm_status = ifelse(grepl("prim", Sample_ID), "P", "R")) %>%
  mutate(sample_num = gsub("GBM(.*)_.*","\\1", Sample_ID)) %>%
  arrange(match(sample_num, sample_num[order(strtoi(sample_num))])) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = unique(Sample_ID)))

# Rich barplot of class2 peptide yields
yields_classII %>% 
  ggplot(., aes(x = Sample_ID, y = n)) +
  geom_bar(stat="identity", fill = ifelse(yields_classII$gbm_status == "P", brewer.pal(3, 'RdBu')[3], brewer.pal(3, 'RdBu')[1])) +
  geom_hline(yintercept = mean(yields_classII$n), linetype= "dashed", alpha = 0.5) +
  annotate(geom = "text", x = seq_len(nrow(yields_classII)), y = -300, label = yields_classII$gbm_status, size = 4) + 
  annotate(geom = "text", x = 2 * (1:(nrow(yields_classII)/2)) - 0.5, y = -700, label = unique(yields_classII$sample_num), size = 4) + 
  annotate(geom = "text", x = -1, y = -700, label = "GBM", size = 4) + 
  coord_cartesian(ylim=c(0,10000), clip="off") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("Number of peptides")

# ----------------------------------
#       Length distribution
# --------------------------------------------------------------------------------
# TODO: The column creation can be refactorized towards more appealing code
# --------------------------------------------------------------------------------
# CLASS 1
# ---------
# Prepare data for plotting
len_dist_classI <- aggregate(binders_classI$Sequence, by=list(Sample_ID = binders_classI$Sample_ID), FUN=nchar) %>%
  mutate(sample_num = gsub("GBM(.*)_.*","\\1", Sample_ID)) %>%
  arrange(match(sample_num, sample_num[order(strtoi(sample_num))])) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = unique(Sample_ID))) %>%
  add_column(eight = NA) %>%
  add_column(nine = NA) %>%
  add_column(ten = NA) %>%
  add_column(eleven = NA) %>%
  add_column(twelve = NA)

# Calculate relative abundances of peptide length
for (i in 1:nrow(len_dist_classI)) {
  counts = table(len_dist_classI$x[i])
  len_dist_classI$eight[i] = counts[1] / sum(counts)
  len_dist_classI$nine[i] = counts[2] / sum(counts)
  len_dist_classI$ten[i] = counts[3] / sum(counts)
  len_dist_classI$eleven[i] = counts[4] / sum(counts)
  len_dist_classI$twelve[i] = counts[5] / sum(counts)
}

# Drop not needed cols
len_dist_classI$x <- NULL

# Transform dataframe into long format such that it can be used in ggplot
len_dist_classI <- len_dist_classI %>% 
  gather("len", "ratio", -c(Sample_ID, sample_num)) %>%
  mutate(len = replace(len, which(len == "eight"), "8-mers")) %>%
  mutate(len = replace(len, which(len == "nine"), "9-mers")) %>%
  mutate(len = replace(len, which(len == "ten"), "10-mers")) %>%
  mutate(len = replace(len, which(len == "eleven"), "11-mers")) %>%
  mutate(len = replace(len, which(len == "twelve"), "12-mers")) %>%
  mutate(len = factor(len, levels = unique(len))) %>%
  arrange(match(sample_num, sample_num[order(strtoi(sample_num))])) 

# Plot Class1 length distribution
len_dist_classI %>%
  ggplot(., aes(x = len, y = ratio, group = Sample_ID , color = Sample_ID)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(11, 'RdBu')[sample.int(11,44, replace = TRUE)]) +
  coord_cartesian(ylim=c(0,1)) +
  labs(x = "", y = "% of HLA class I ligands") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

# ---------
# CLASS 2
# ---------
# Prepare data for plotting
len_dist_classII <- aggregate(peptides_classII$Sequence, by=list(Sample_ID = peptides_classII$Sample_ID), FUN=nchar) %>%
  mutate(sample_num = gsub("GBM(.*)_.*","\\1", Sample_ID)) %>%
  arrange(match(sample_num, sample_num[order(strtoi(sample_num))])) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = unique(Sample_ID))) %>%
  add_column(eight = NA) %>%
  add_column(nine = NA) %>%
  add_column(ten = NA) %>%
  add_column(eleven = NA) %>%
  add_column(twelve = NA) %>%
  add_column(thirteen = NA) %>%
  add_column(fourteen = NA) %>%
  add_column(fiveteen = NA) %>%
  add_column(sixteen = NA) %>%
  add_column(seventeen = NA) %>%
  add_column(eighteen = NA) %>%
  add_column(nineteen = NA) %>%
  add_column(twenty = NA) %>%
  add_column(twentyone = NA) %>%
  add_column(twentytwo = NA) %>%
  add_column(twentythree = NA) %>%
  add_column(twentyfour = NA) %>%
  add_column(twentyfive = NA) 

# Calculate relative abundances of peptide length
for (i in 1:nrow(len_dist_classII)) {
  counts = table(len_dist_classII$x[i])
  len_dist_classII$eight[i] = counts[1] / sum(counts)
  len_dist_classII$nine[i] = counts[2] / sum(counts)
  len_dist_classII$ten[i] = counts[3] / sum(counts)
  len_dist_classII$eleven[i] = counts[4] / sum(counts)
  len_dist_classII$twelve[i] = counts[5] / sum(counts)
  len_dist_classII$thirteen[i] = counts[6] / sum(counts)
  len_dist_classII$fourteen[i] = counts[7] / sum(counts)
  len_dist_classII$fiveteen[i] = counts[8] / sum(counts)
  len_dist_classII$sixteen[i] = counts[9] / sum(counts)
  len_dist_classII$seventeen[i] = counts[10] / sum(counts)
  len_dist_classII$eighteen[i] = counts[11] / sum(counts)
  len_dist_classII$nineteen[i] = counts[12] / sum(counts)
  len_dist_classII$twenty[i] = counts[13] / sum(counts)
  len_dist_classII$twentyone[i] = counts[14] / sum(counts)
  len_dist_classII$twentytwo[i] = counts[15] / sum(counts)
  len_dist_classII$twentythree[i] = counts[16] / sum(counts)
  len_dist_classII$twentyfour[i] = counts[17] / sum(counts)
  len_dist_classII$twentyfive[i] = counts[18] / sum(counts)
}

# Drop not needed cols
len_dist_classII$x <- NULL

# Transform dataframe into long format such that it can be used in ggplot
len_dist_classII <- len_dist_classII %>% 
  gather("len", "ratio", -c(Sample_ID, sample_num)) %>%
  mutate(len = replace(len, which(len == "eight"), "8-mers")) %>%
  mutate(len = replace(len, which(len == "nine"), "9-mers")) %>%
  mutate(len = replace(len, which(len == "ten"), "10-mers")) %>%
  mutate(len = replace(len, which(len == "eleven"), "11-mers")) %>%
  mutate(len = replace(len, which(len == "twelve"), "12-mers")) %>%
  mutate(len = replace(len, which(len == "thirteen"), "13-mers")) %>%
  mutate(len = replace(len, which(len == "fourteen"), "14-mers")) %>%
  mutate(len = replace(len, which(len == "fiveteen"), "15-mers")) %>%
  mutate(len = replace(len, which(len == "sixteen"), "16-mers")) %>%
  mutate(len = replace(len, which(len == "seventeen"), "17-mers")) %>%
  mutate(len = replace(len, which(len == "eighteen"), "18-mers")) %>%
  mutate(len = replace(len, which(len == "nineteen"), "19-mers")) %>%
  mutate(len = replace(len, which(len == "twenty"), "20-mers")) %>%
  mutate(len = replace(len, which(len == "twentyone"), "21-mers")) %>%
  mutate(len = replace(len, which(len == "twentytwo"), "22-mers")) %>%
  mutate(len = replace(len, which(len == "twentythree"), "23-mers")) %>%
  mutate(len = replace(len, which(len == "twentyfour"), "24-mers")) %>%
  mutate(len = replace(len, which(len == "twentyfive"), "25-mers")) %>%
  mutate(len = factor(len, levels = unique(len))) %>%
  arrange(match(sample_num, sample_num[order(strtoi(sample_num))])) 

# Plot Class2 length distribution
len_dist_classII %>%
  ggplot(., aes(x = len, y = ratio, group = Sample_ID , color = Sample_ID)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(11, 'RdBu')[sample.int(11,44, replace = TRUE)]) +
  coord_cartesian(ylim=c(0,0.3)) +
  labs(x = "", y = "% of HLA class II ligands") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

# ------------------------------------
#         Waterfall plots 
# --------------------------------------------------------------------------------------------------------
# A waterfall plot is an ordered barplot of each peptide with respect to their ocurrence in 2 conditions
# On both ends idealy one can identify the exlcusive peptides of each group
# --------------------------------------------------------------------------------------------------------
# Function which creates the dataframe for the waterfall plot. This step might take a few minutes
compute_waterfall <- function(df1, df2){
  # Get union of unique peptides
  sequences <- unique(c(df1$Sequence, df2$Sequence))
  # dataframes need to have unique sequences per sample since the frequency of positive ligandomes is observed
  df1 <- df1 %>%
    group_by(Sample_ID) %>%
    summarise(Sequence = unique(Sequence))
  df2 <- df2 %>%
    group_by(Sample_ID) %>%
    summarise(Sequence = unique(Sequence))
  # Create data frame to fill accordingly
  waterfall <- data.frame(sequences)
  # Calculate the ratio of peptide occurences per sample. This step might take a while
  waterfall$group1 <- sapply(sequences, function(seq){
    mean(ifelse(is.na(matchAll(seq, df1$Sequence)), 0,
                (length(matchAll(seq, df1$Sequence))/length(unique(df1$Sample_ID)))*100))
  })
  waterfall$group2 <- sapply(sequences, function(seq){
    mean(ifelse(is.na(matchAll(seq, df2$Sequence)), 0,
                (length(matchAll(seq, df2$Sequence))/length(unique(df2$Sample_ID)))*100))
  })
  waterfall$total <- sapply(sequences, function(seq){
    (length(matchAll(seq, c(df2$Sequence, df1$Sequence)))/
       length(unique(c(df2$Sample_ID, df1$Sample_ID))))*100
  })
  
  # Order data according to the desired input of a waterfall plot
  exclusive_group1 <- waterfall[which(waterfall$group2 == 0), ]
  exclusive_group1 <- exclusive_group1[order(exclusive_group1$group1, decreasing = T), ]
  exclusive_group2 <- waterfall[which(waterfall$group1 == 0), ]
  exclusive_group2 <- exclusive_group2[order(exclusive_group2$group2), ]
  other <- waterfall[which(waterfall$group1 != 0 & waterfall$group2 != 0), ]
  other <- other[order(other$group1, other$group2), ]
  
  waterfall <- rbind(exclusive_group1, other, exclusive_group2)
  # Filter out one-hit-wonders
  sample_size <- length(unique(df1$Sample_ID)) + length(unique(df2$Sample_ID)) 
  waterfall <- waterfall[which((waterfall$total/100) * sample_size > 1), ]
  
  return(waterfall)
}

# Function which takes the output of compute_waterfall and plots the data accordingly
plot_waterfall <- function(wtf, pth){
  blue <- brewer.pal(5, 'RdBu')[5]
  red <- brewer.pal(5, 'RdBu')[1]
  color_legend <- c("GBM" = blue, "Benign" = red)
  w <- wtf %>%
    # Manually determine the order of sequences on the x-axis since ggplot orders lexiographically by default
    mutate(sequences = factor(sequences, levels = unique(sequences))) %>%
    ggplot(aes(x=sequences)) +
    ## Retrieve all of the different fasta that are in the directory 
    geom_bar( aes(y = group1, fill = "GBM"), stat="identity", width = 1) +
    geom_bar( aes(y = -group2, fill = "Benign"), stat="identity", width = 1) +
    # Manually determine the y ticks
    scale_y_continuous(limits = c(-70,70), breaks = seq(-100, 100, 10)) +
    scale_color_manual(values = color_legend) +
    scale_fill_discrete(breaks = c("GBM", "Benign")) +
    ## Change the text on the y axis 
    ylab("Frequency of positive ligandomes (%)") +
    xlab("Sequences") +
    theme_bw() + 
    ## Make sure that the background is white, and that the text is reduced to hardly visible
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.title = element_blank()) 
    
  
  ggsave(pth, w, width = 16, height = 8)
}

# Function plotting a venn diagram of 2 or 3 sets
plot_venny <- function(ls, ls_names){
  a1 <- unlist(ls[1])
  a2 <- unlist(ls[2])
  n12 <- length(intersect(a1, a2))
  if (length(ls) == 2){
    v <- draw.pairwise.venn(length(a1), length(a2), n12,
                            fill = c(brewer.pal(9, 'RdBu')[2], brewer.pal(9, 'YlGnBu')[5]),
                            lty = "blank", scale = F, euler.d = F, cex = 2.5)
    return(v)
  } else {
    a3 <- unlist(ls[3])
    n13 <- length(intersect(a1, a3))
    n23 <- length(intersect(a2, a3))
    n123 <- length(Reduce(intersect, ls))
    v <- draw.triple.venn(length(a1), length(a2), length(a3), n12, n13, n23, n123,
                          fill = c(brewer.pal(7, 'RdBu')[7], brewer.pal(7, 'YlGnBu')[3], brewer.pal(5, 'RdBu')[2]),
                            category = ls_names, lty = "blank",
                          cat.dist = c(.05,.05,.03), cat.pos = c(340,20, 180))
    return(v)
  }
}


# ------------------------------------
#     Waterfall GBM vs benign DB
# ------------------------------------
# ---------
# CLASS 1
# ---------
waterfall_gbm_vs_benign_classI <- compute_waterfall(binders_classI, benign_DB_classI)
# or load DF since it took hours:
#load(file = "waterfall_gbm_vs_benign_classI.Rdata")
plot_waterfall(waterfall_gbm_vs_benign_classI, "figures/classI_gbm_vs_benign_waterfallplot.jpeg")
gbm_exclusives_classI <- waterfall_gbm_vs_benign_classI$sequences[which(waterfall_gbm_vs_benign_classI$group2 == 0)]

# Plot Vennie (and export manually as pdf with 4x4 width and height since saving w/ this package is shit)
binders_classI_prim <- binders_classI[which(binders_classI$GBM_status == "Primary GBM"), ]
binders_classI_rec <- binders_classI[which(binders_classI$GBM_status == "Recurrent GBM"), ]
plot_venny(list(unique(binders_classI_prim$Sequence), unique(binders_classI_rec$Sequence), unique(benign_DB_classI$Sequence)),
           c("Primary GBM", "Recurrent GBM", "Benign"))
dev.off()

# ---------
# CLASS 2
# ---------
waterfall_gbm_vs_benign_classII <- compute_waterfall(peptides_classII, benign_DB_classII)
# or load DF since it took hours:
#load(file = "waterfall_gbm_vs_benign_classII.Rdata")
plot_waterfall(waterfall_gbm_vs_benign_classII, "figures/classII_gbm_vs_benign_waterfallplot.jpeg")
gbm_exclusives_classII <- waterfall_gbm_vs_benign_classII$sequences[which(waterfall_gbm_vs_benign_classII$group2 == 0)]

# Plot Vennie (and export manually as pdf with 4x4 width and height since saving w/ this package is shit)
peptides_classII_prim <- peptides_classII[which(peptides_classII$GBM_status == "Primary GBM"), ]
peptides_classII_rec <- peptides_classII[which(peptides_classII$GBM_status == "Recurrent GBM"), ]
plot_venny(list(unique(peptides_classII_prim$Sequence), unique(peptides_classII_rec$Sequence),
                unique(benign_DB_classII$Sequence)),
           c("Primary GBM", "Recurrent GBM", "Benign"))
dev.off()

# ------------------------------------------------------
# Waterfall plot exclusive GBMs: primary vs recurrent
# ------------------------------------------------------
# ---------
# CLASS 1
# ---------
excl_binders_classI_prim <- binders_classI[which(binders_classI[match(gbm_exclusives_classI, binders_classI$Sequence), ]$GBM_status == "Primary GBM"), ]
excl_binders_classI_rec <- binders_classI[which(binders_classI[match(gbm_exclusives_classI, binders_classI$Sequence), ]$GBM_status == "Recurrent GBM"), ]

waterfall_prim_vs_rec_classI <- compute_waterfall(excl_binders_classI_prim, excl_binders_classI_rec)
plot_waterfall(waterfall_prim_vs_rec_classI, "figures/classI_gbm_exclusives_p_vs_r_waterfallplot.jpeg")

# Plot Vennie (and export manually as pdf with 4x4 width and height since saving w/ this package is shit)
plot_venny(list(unique(excl_binders_classI_prim$Sequence), unique(excl_binders_classI_rec$Sequence)),
           c("Primary GBM", "Recurrent GBM"))
dev.off()

# ---------
# CLASS 2
# ---------
excl_peptides_classII_prim <- peptides_classII[which(peptides_classII[match(gbm_exclusives_classII, peptides_classII$Sequence), ]$GBM_status == "Primary GBM"), ]
excl_peptides_classII_rec <- peptides_classII[which(peptides_classII[match(gbm_exclusives_classII, peptides_classII$Sequence), ]$GBM_status == "Recurrent GBM"), ]

waterfall_prim_vs_rec_class2 <- compute_waterfall(excl_peptides_classII_prim, excl_peptides_classII_rec)
plot_waterfall(waterfall_prim_vs_rec_class2, "figures/classII_gbm_exclusives_p_vs_r_waterfallplot.jpeg")

# Plot Vennie (and export manually as pdf with 4x4 width and height since saving w/ this package is shit)
plot_venny(list(unique(excl_peptides_classII_prim$Sequence), unique(excl_peptides_classII_rec$Sequence)),
           c("Primary GBM", "Recurrent GBM"))
dev.off()

# ------------------------------------
#       Saturation analysis
# ------------------------------------------------------------------------------------------------------------
# At which sample size does the increase in number of proteins or unique peptides stagnate?
# This is archieved through a non-linear least-squared (nls) fit on the relation of #samples and #unique peptides
# ------------------------------------------------------------------------------------------------------------
# Plot saturation curve
plot_saturation <- function(df, fit_data, fit_params){
  p <- df %>%
    ggplot(aes(x = n_samples, y = fit_data)) +
    geom_smooth(method = "loess", se = F, color = "black", formula = y ~ x) +
    geom_hline(yintercept = abs(fit_params[1])+fit_params[2], linetype = "dashed", alpha = 0.6) + 
    labs(x = "Number of samples", y = "Number of unique peptides") +
    theme(axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
          axis.text.y = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent")) # get rid of legend b
  return(p)
}
# ---------------
# Unique Peptides
# ---------------
# CLASS 1
# ---------
saturation_classI <- binders_classI[match(unique(binders_classI$Sequence), binders_classI$Sequence), ] %>%
  group_by(Sample_ID) %>%
  summarise(n_unique = length(unique(Sequence)))  %>%
  arrange(desc(n_unique))  %>%
  mutate(n_unique = cumsum(n_unique)) %>%
  mutate(n_samples = seq(1, length(Sample_ID)))

fit_params_saturation_classI <- NLSstAsymptotic(sortedXyData(c(0,saturation_classI$n_samples), c(0,saturation_classI$n_unique)))
fit_saturation_classI <- SSasymp(saturation_classI$n_samples, abs(fit_params_saturation_classI[1])+fit_params_saturation_classI[2], fit_params_saturation_classI[1], fit_params_saturation_classI[3])

ggsave(filename = "figures/classI_saturation_peptides.jpeg",
       plot_saturation(saturation_classI, fit_saturation_classI, fit_params_saturation_classI),
       width = 5, height = 4)

# ---------
# CLASS 2
# ---------
saturation_classII <- peptides_classII[match(unique(peptides_classII$Sequence), peptides_classII$Sequence), ] %>%
  group_by(Sample_ID) %>%
  summarise(n_unique = length(unique(Sequence)))  %>%
  arrange(desc(n_unique))  %>%
  mutate(n_unique = cumsum(n_unique)) %>%
  mutate(n_samples = seq(1, length(Sample_ID)))

fit_params_saturation_classII <- NLSstAsymptotic(sortedXyData(c(0,saturation_classII$n_samples), c(0,saturation_classII$n_unique)))
fit_saturation_classII <- SSasymp(saturation_classII$n_samples, abs(fit_params_saturation_classII[1])+fit_params_saturation_classII[2], fit_params_saturation_classII[1], fit_params_saturation_classII[3])

ggsave(filename = "figures/classII_saturation_peptides.jpeg",
       plot_saturation(saturation_classII, fit_saturation_classII, fit_params_saturation_classII),
       width = 5, height = 4)


# --------------------------------------------------------
#    Unique peptides on primary and recurrent level
# --------------------------------------------------------
# ----------------
# CLASS 1 primary
# ----------------
binders_classI_prim <- binders_classI[which(binders_classI$GBM_status == "Primary GBM"), ]
                                     
saturation_classI_prim <- binders_classI_prim[match(unique(binders_classI_prim$Sequence), binders_classI_prim$Sequence), ] %>%
  group_by(Sample_ID) %>%
  summarise(n_unique = length(unique(Sequence)))  %>%
  arrange(desc(n_unique))  %>%
  mutate(n_unique = cumsum(n_unique)) %>%
  mutate(n_samples = seq(1, length(Sample_ID)))

fit_params_saturation_classI_prim <- NLSstAsymptotic(sortedXyData(c(0,saturation_classI_prim$n_samples), c(0,saturation_classI_prim$n_unique)))
fit_saturation_classI_prim <- SSasymp(saturation_classI_prim$n_samples, abs(fit_params_saturation_classI_prim[1])+fit_params_saturation_classI_prim[2], fit_params_saturation_classI_prim[1], fit_params_saturation_classI_prim[3])

ggsave(filename = "figures/classI_saturation_prim_peptides.jpeg",
       plot_saturation(saturation_classI_prim, fit_saturation_classI_prim, fit_params_saturation_classI_prim),
       width = 5, height = 4)

# ----------------
# CLASS 1 recurrent
# ----------------
binders_classI_rec <- binders_classI[which(binders_classI$GBM_status == "Recurrent GBM"), ]

saturation_classI_rec <- binders_classI_rec[match(unique(binders_classI_rec$Sequence), binders_classI_rec$Sequence), ] %>%
  group_by(Sample_ID) %>%
  summarise(n_unique = length(unique(Sequence)))  %>%
  arrange(desc(n_unique))  %>%
  mutate(n_unique = cumsum(n_unique)) %>%
  mutate(n_samples = seq(1, length(Sample_ID)))

fit_params_saturation_classI_rec <- NLSstAsymptotic(sortedXyData(c(0,saturation_classI_rec$n_samples), c(0,saturation_classI_rec$n_unique)))
fit_saturation_classI_rec <- SSasymp(saturation_classI_rec$n_samples, abs(fit_params_saturation_classI_rec[1])+fit_params_saturation_classI_rec[2], fit_params_saturation_classI_rec[1], fit_params_saturation_classI_rec[3])

ggsave(filename = "figures/classI_saturation_rec_peptides.jpeg",
       plot_saturation(saturation_classI_rec, fit_saturation_classI_rec, fit_params_saturation_classI_rec),
       width = 5, height = 4)


# ----------------
# CLASS 2 primary
# ----------------
peptides_classII_prim <- peptides_classII[which(peptides_classII$GBM_status == "Primary GBM"), ]

saturation_classII_prim <- peptides_classII_prim[match(unique(peptides_classII_prim$Sequence), peptides_classII_prim$Sequence), ] %>%
  group_by(Sample_ID) %>%
  summarise(n_unique = length(unique(Sequence)))  %>%
  arrange(desc(n_unique))  %>%
  mutate(n_unique = cumsum(n_unique)) %>%
  mutate(n_samples = seq(1, length(Sample_ID)))

fit_params_saturation_classII_prim <- NLSstAsymptotic(sortedXyData(c(0,saturation_classII_prim$n_samples), c(0,saturation_classII_prim$n_unique)))
fit_saturation_classII_prim <- SSasymp(saturation_classII_prim$n_samples, abs(fit_params_saturation_classII_prim[1])+fit_params_saturation_classII_prim[2], fit_params_saturation_classII_prim[1], fit_params_saturation_classII_prim[3])

ggsave(filename = "figures/classII_saturation_prim_peptides.jpeg",
       plot_saturation(saturation_classII_prim, fit_saturation_classII_prim, fit_params_saturation_classII_prim),
       width = 5, height = 4)

# ----------------
# CLASS 2 recurrent
# ----------------
peptides_classII_rec <- peptides_classII[which(peptides_classII$GBM_status == "Recurrent GBM"), ]

saturation_classII_rec <- peptides_classII_rec[match(unique(peptides_classII_rec$Sequence), peptides_classII_rec$Sequence), ] %>%
  group_by(Sample_ID) %>%
  summarise(n_unique = length(unique(Sequence)))  %>%
  arrange(desc(n_unique))  %>%
  mutate(n_unique = cumsum(n_unique)) %>%
  mutate(n_samples = seq(1, length(Sample_ID)))

fit_params_saturation_classII_rec <- NLSstAsymptotic(sortedXyData(c(0,saturation_classII_rec$n_samples), c(0,saturation_classII_rec$n_unique)))
fit_saturation_classII_prim <- SSasymp(saturation_classII_rec$n_samples, abs(fit_params_saturation_classII_rec[1])+fit_params_saturation_classII_rec[2], fit_params_saturation_classII_rec[1], fit_params_saturation_classII_rec[3])

ggsave(filename = "figures/classII_saturation_rec_peptides.jpeg",
       plot_saturation(saturation_classII_rec, fit_saturation_classII_prim, fit_params_saturation_classII_rec),
       width = 5, height = 4)


# ----------------------------------------
#              Vulcano plots
# -------------------------------------------------------------------------------------------------------
# Differential expression and fold change between primary and recurrent per donor can be computed using
# the Peptide Spectrum Matches (PSMs) of the intersection of primary and recurrent peptides per donor
# -------------------------------------------------------------------------------------------------------
## A customized intensity scaling function using median robust scaling
## REF: QBiC_workflow_result_parser.R
scaleMaps <- function(maps){
  ## Scale intensities
  scaling_factors <- max(maps)/colSums(maps)
  maps <- t(t(maps) * scaling_factors)
  ## Return the dataframe
  return(maps)
}

## p-value computation
## REF: QBiC_workflow_result_parser.R 
compute_pval <- function(con1, con2) {
  ## Always apply t.test (other tests will not work)
  pvals <- t.test(con1,con2)$p.value
  ## NA values are changed to a value of 1
  pvals[is.na(pvals)] <- 1
  ## Transform to values as -log10
  pvals_adj <- -log10(p.adjust(pvals, method="BH"))
  ## Return the data frame
  return(pvals_adj)
}

## Compute log2 foldchange 
## REF: QBiC_workflow_result_parser.R 
compute_fch <- function(con1, con2) {
  ## Get the mean of the different conditions 
  fch <- log2(10**(median(as.numeric(con1)) - median(as.numeric(con2))))
  ## Return the log2foldchange
  return(fch)
}

