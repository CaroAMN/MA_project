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



## DESCRIPTION: function that rewrites HLA types to a format that is compatible
##              with netMHCpan HLA type input
## PARAMETERS:
##             - uniqe_HLA_types: vector with unique HLA types containing * [vector]
##             - as_string: set to FALSE will return new HLA types as vector
##                           set to TRUE will return a single string of new HLA
##                           types which can be used directly as netMHCpan HLA
##                           type input[boolean]
##             
## OUTPUT: either vector or string of reformatted HLA types from A*01:01
##         to HLA-A01:01
##
rewrite_HLA_types <- function(uniqe_HLA_types,as_string) {
  
  new_HLA_types <- c()
  for (i in uniqe_HLA_types) {
    str2 = strsplit(i, "*", fixed = TRUE)
    new_HLA_types = append(rewritten_hlatypes,
                           paste0("HLA-", str2[[1]][1],str2[[1]][2]))
  }
  if (as_string == TRUE) {
    HLA_str = ""
    for (i in unique(new_HLA_types)) {
      HLA_str = paste0(HLA_str, i, ",")
    }
    return(HLA_str)

  } else
    return(unique(new_HLA_types))
}





## DESCRIPTION: function to read in ligandomics data from csv files
##            
## PARAMETERS: - 
##             
## OUTPUT: 
##
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
    file_list <- strsplit(file,split = "/")[[1]]
    patient_id <- strsplit(grep("ZH",file_list, value = TRUE), split = "_")[[1]][1]
    t_region <- strsplit(grep("ZH",file_list, value = TRUE), split = "_")[[1]][2]
    seq <- tempfile$sequence
    acc <- tempfile$accessions
    # temporaray df
    temp_df <- data.frame(
      Patient_ID = rep(patient_id, times = length(seq)),
      Sample_num = rep(paste("Sample", as.character(file_count), sep = "_"), times = length(seq)),
      Tumor_region = rep(t_region, times = length(seq)),
      Sequence = seq,
      Accessions = acc
    )
    # append to output df
    df <- rbind(df, temp_df)
  }
  df$Sequence <- str_replace_all(df$Sequence, "\\(Oxidation\\)", "") # get rid of the (Oxidation) string in sequences
  df$Tumor_region <- str_replace_all(df$Tumor_region,"\\.csv", "" )
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

