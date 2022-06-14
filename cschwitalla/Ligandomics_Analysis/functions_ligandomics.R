################################################################################
###                                Purpose                                   ###
################################################################################
## Title: Functions for Immunopeptidomics Analysis
## Author: Carolin Schwitalla
##
## Description: This script contains all functions that are used for the
##              immunopeptidomics analysis for my master thesis
################################################################################
###                                                                          ###
################################################################################

# function to read in ligandomics data from csv files
createDataFrame <- function(filelist) {
  # create empty classI dataframe
  file_count <- 0
  df <- data.frame(
    Patient_ID = character(),
    Tumor_region = character(),
    Sequence = character(),
    Accessions = character()
  )
  # loop over all files
  for (file in filelist) {
    file_count <- file_count + 1
    tempfile <- read.csv(file, header = TRUE, sep = "\t")
    patientid <- strsplit(strsplit(file, split = "/")[[1]][3], split = "_")[[1]][1]
    tregion <- strsplit(strsplit(strsplit(file, split = "/")[[1]][3], split = "_")[[1]][2], split = ".", fixed = TRUE)[[1]][1]
    seq <- tempfile$sequence
    acc <- tempfile$accessions
    # temporaray df
    temp_df <- data.frame(
      Patient_ID = rep(patientid, times = length(seq)),
      Sample_num = rep(paste("Sample", as.character(file_count), sep = "_"), times = length(seq)),
      Tumor_region = rep(tregion, times = length(seq)),
      Sequence = seq,
      Accessions = acc
    )
    # append to output df
    df <- rbind(df, temp_df)
  }
  df$Sequence <- str_replace_all(df$Sequence, "\\(Oxidation\\)", "") # get rid of the (Oxidation) string in sequences
  df$Sequence <- toupper(df$Sequence) # sequence all upper case
  return(df)
}



getProteinAcc <- function(list){
  Pacc = c()
  acc_only = c()
  
  for( i in list){
    str = strsplit(i, ";")
    Pacc = append(Pacc, str[[1]])
  }
  for (i in Pacc){
    str = strsplit(i, "|", fixed = TRUE)
    acc_only = append(acc_only, str[[1]][2])
  }
  return(acc_only)
}

getProteinAcc_uniqemappers <- function(list){
  acc_only = c()
  for (i in list){
    str = strsplit(i, "|", fixed = TRUE)
    acc_only = append(acc_only, str[[1]][2])
  }
  return(acc_only)
  
}

