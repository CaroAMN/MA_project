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
## TITLE: make multi Waterfall df ----------------------------------------------------
## DESCRIPTION:
##
## PARAMETERS: -
##
## OUTPUT:
##
make_multi_waterfall_df <- function(df1, df2, df3){
  
  if (with_seq == TRUE) {
    uniq_seq <- unique(c(df1$Sequence, df2$Sequence, df3$Sequence))
    
    #unique sequence per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample 
    df1 = df1 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
    df2 = df2 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
    df3 = df3 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
    # neues df 
    waterfall_df = data.frame(uniq_seq)
    # ratio of peptide occurences für NEC und INF 
    waterfall_df$set1 = sapply(uniq_seq, function(seq){
      mean(ifelse(is.na(matchAll(seq, df1$Sequence)), 0,
                  (length(matchAll(seq, df1$Sequence))/length(unique(df1$Patient_ID)))*100))
    })
    
    waterfall_df$set2 = sapply(uniq_seq, function(seq){
      mean(ifelse(is.na(matchAll(seq, df2$Sequence)), 0,
                  (length(matchAll(seq, df2$Sequence))/length(unique(df2$Patient_ID)))*100))
    })
    
    waterfall_df$set3 = sapply(uniq_seq, function(seq){
      mean(ifelse(is.na(matchAll(seq, df3$Sequence)), 0,
                  (length(matchAll(seq, df3$Sequence))/length(unique(df3$Patient_ID)))*100))
    })
    
    
    
    waterfall_df$Total = sapply(uniq_seq, function(seq){
      (length(matchAll(seq, c(df1$Sequence, df2$Sequence, df3$Sequence)))/
         length(c(unique(df1$Patient_ID), unique(df2$Patient_ID), unique(df3$Patient_ID))))*100
    })
    
    
    
  } else (with_seq == FALSE)
  
  uniq_acc <- unique(c(df1$Accessions, df2$Accessions, df3$Accessions))
  
  #unique sequence per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample 
  df1 = df1 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  df2 = df2 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  df3 = df3 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  # neues df 
  waterfall_df = data.frame(uniq_acc)
  # ratio of peptide occurences für NEC und INF 
  waterfall_df$set1 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df1$Accessions)), 0,
                (length(matchAll(acc, df1$Accessions))/length(unique(df1$Patient_ID)))*100))
  })
  
  waterfall_df$set2 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df2$Accessions)), 0,
                (length(matchAll(acc, df2$Accessions))/length(unique(df2$Patient_ID)))*100))
  })
  
  waterfall_df$set3 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df3$Accessions)), 0,
                (length(matchAll(acc, df3$Accessions))/length(unique(df3$Patient_ID)))*100))
  })
  
  
  
  waterfall_df$Total = sapply(uniq_acc, function(acc){
    (length(matchAll(acc, c(df1$Accessions, df2$Accessions, df3$Accessions)))/
       length(c(unique(df1$Patient_ID), unique(df2$Patient_ID), unique(df3$Patient_ID))))*100
  })
  
  # filtern des df 
  exclusive_1 = waterfall_df[which(waterfall_df$set2 == 0 & waterfall_df$set3 == 0 ),]
  exclusive_2 = waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set3 == 0 ),]
  exclusive_3 = waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set2 == 0 ),]
  #exclusive_4 = waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set2 == 0 & waterfall_df$set3 == 0),]
  
  shared = waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0 & waterfall_df$set3 != 0),]
  all = waterfall_df

  #create new waterfall df waterfall_df = rbind(exclusive_1, exclusive_2, exclusive_3, exclusive_4, shared)
  waterfall_df2 = list(Nec = exclusive_1, T1 = exclusive_2, INF = exclusive_3, Shared = shared, all = all)
  return(waterfall_df2)
  
}




## TITLE: make Waterfall df ----------------------------------------------------
## DESCRIPTION:
##
## PARAMETERS: -
##
## OUTPUT:
##
# plots with peptide sequences
make_waterfall_df <- function(df1, df2, with_seq) {
  if (with_seq == TRUE) {
    # union of unique sequencs --> x achse
    uniq_seq <- unique(c(df1$Sequence, df2$Sequence))
    # unique sequence per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample
    df1 <- df1 %>%
      group_by(Patient_ID) %>%
      summarise(Sequence = unique(Sequence))
    df2 <- df2 %>%
      group_by(Patient_ID) %>%
      summarise(Sequence = unique(Sequence))
    # neues df
    waterfall_df <- data.frame(uniq_seq)
    # ratio of peptide occurences für NEC und INF
    waterfall_df$set1 <- sapply(uniq_seq, function(seq) {
      mean(ifelse(is.na(tuple::matchAll(seq, df1$Sequence)), 0,
        (length(tuple::matchAll(seq, df1$Sequence)) / length(unique(df1$Patient_ID))) * 100
      ))
    })

    waterfall_df$set2 <- sapply(uniq_seq, function(seq) {
      mean(ifelse(is.na(tuple::matchAll(seq, df2$Sequence)), 0,
        (length(tuple::matchAll(seq, df2$Sequence)) / length(unique(df2$Patient_ID))) * 100
      ))
    })

    waterfall_df$Total <- sapply(uniq_seq, function(seq) {
      (length(tuple::matchAll(seq, c(df1$Sequence, df2$Sequence))) /
        length(c(unique(df1$Patient_ID), unique(df2$Patient_ID)))) * 100
    })
  } else {
    (with_seq == FALSE)
  }
  # unique acc per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample
  df1 <- df1 %>%
    group_by(Patient_ID) %>%
    summarise(Accessions = unique(Accessions))
  df2 <- df2 %>%
    group_by(Patient_ID) %>%
    summarise(Accessions = unique(Accessions))

  # union of unique acc --> x achse
  uniq_acc <- unique(c(df1$Accessions, df2$Accessions))

  # neues df
  waterfall_df <- data.frame(uniq_acc)
  # ratio of peptide occurences für df1 und df2
  waterfall_df$set1 <- sapply(uniq_acc, function(acc) {
    mean(ifelse(is.na(tuple::matchAll(acc, df1$Accessions)), 0,
      (length(tuple::matchAll(acc, df1$Accessions)) / length(unique(df1$Patient_ID))) * 100
    ))
  })

  waterfall_df$set2 <- sapply(uniq_acc, function(acc) {
    mean(ifelse(is.na(tuple::matchAll(acc, df2$Accessions)), 0,
      (length(tuple::matchAll(acc, df2$Accessions)) / length(unique(df2$Patient_ID))) * 100
    ))
  })

  waterfall_df$Total <- sapply(uniq_acc, function(acc) {
    (length(tuple::matchAll(acc, c(df1$Accessions, df2$Accessions))) /
      length(c(unique(df1$Patient_ID), unique(df2$Patient_ID)))) * 100
  })


  # filtern des df
  exclusive_1 <- waterfall_df[which(waterfall_df$set2 == 0), ]
  exclusive_2 <- waterfall_df[which(waterfall_df$set1 == 0), ]
  # ordnen
  exclusive_1 <- exclusive_1[order(exclusive_1$set1, decreasing = T), ]
  exclusive_2 <- exclusive_2[order(exclusive_2$set2), ]

  other <- waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0), ]
  other <- other[order(other$set1, other$set2), ]
  # create new waterfall df
  waterfall_df <- rbind(exclusive_1, other, exclusive_2)
  return(waterfall_df)
}


## TITLE: plot Waterfall plot---------------------------------------------------
## DESCRIPTION:
##
## PARAMETERS: -
##
## OUTPUT:
##
plot_waterfall <- function(waterfall_df, name1, name2, title, with_seq) {
  if ( with_seq == TRUE) {
    w <- waterfall_df %>% mutate(uniq_seq = factor(uniq_seq, levels = unique(uniq_seq)))
    x <- w$uniq_seq
    x_lab_text <- "Sequences"
    
  } else (with_seq == FALSE)
  
  
  w <-  waterfall_df %>% mutate(uniq_acc = factor(uniq_acc, levels = unique(uniq_acc)))
  x <- w$uniq_acc
  x_lab_text <- "Source Proteins"
 
  
  # plot waterfall plot 
  ggplot(w, aes(x = x)) +
    geom_bar(aes(y = set1, fill = name1), stat = "identity", width = 1) +
    geom_bar(aes(y = -set2, fill = name2), stat = "identity", width = 1) +
    scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 10)) +
    scale_fill_manual(values = c("salmon1", "steelblue3")) +
    ylab("Frequency of positive ligandomes (%)") +
    xlab(x_lab_text) +
    labs(title = title) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}


plot_venn_waterfall <- function(listarea1,listarea2,title) {
  
  VennDiagram::draw.pairwise.venn(area1 = length(unique(listarea1)),area2 = length(unique(listarea2)),
                                  cross.area = length(Reduce(dplyr::intersect, list(unique(listarea1), unique(listarea2)))),
                                  #category = category,
                                  fill = c("salmon1", "steelblue3"),
                                  cat.pos = c(0, 0), cat.dist = c(.05, .05),
                                  cex = 1.5,
                                  cat.cex = 1.5)
}



## TITLE: rewrite HLA types ----------------------------------------------------
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
rewrite_HLA_types <- function(uniqe_HLA_types, as_string) {
  new_HLA_types <- c()
  for (i in uniqe_HLA_types) {
    str2 <- strsplit(i, "*", fixed = TRUE)
    new_HLA_types <- append(
      new_HLA_types,
      paste0("HLA-", str2[[1]][1], str2[[1]][2])
    )
  }
  if (as_string == TRUE) {
    HLA_str <- ""
    for (i in unique(new_HLA_types)) {
      HLA_str <- paste0(HLA_str, i, ",")
    }
    return(HLA_str)
  } else {
    return(unique(new_HLA_types))
  }
}


## TITLE: create df from MHCquant results---------------------------------------
## DESCRIPTION: function to read in ligandomics data from csv files
##
## PARAMETERS: -
##
## OUTPUT:
##
create_data_frame <- function(filelist) {
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
    file_list <- strsplit(file, split = "/")[[1]]
    patient_id <- strsplit(grep("ZH", file_list, value = TRUE), split = "_")[[1]][1]
    t_region <- strsplit(grep("ZH", file_list, value = TRUE), split = "_")[[1]][2]
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
  df$Tumor_region <- str_replace_all(df$Tumor_region, "\\.csv", "")
  df$Sequence <- toupper(df$Sequence) # sequence all upper case
  return(df)
}

## TITLE: plot custom venn -----------------------------------------------------
## DESCRIPTION: function to plot venn diagrams
##
## PARAMETERS: -
##
## OUTPUT:
##
plot_custom_venn <- function(set, title) {
  venn_data <- process_data(Venn(set))
  venn_plot <- ggVennDiagram(set,
    label_alpha = 0,
    label = "count",
    label_size = 4,
    edge_size = 0
  ) +
    scale_fill_gradient(low = "papayawhip", high = "paleturquoise4") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  plot(venn_plot)
  return(venn_data)
}


# getProteinAcc <- function(list){
#   Pacc = c()
#   acc_only = c()
#
#   for( i in list){
#     str = strsplit(i, ";")
#     Pacc = append(Pacc, str[[1]])
#   }
#   for (i in Pacc){
#     str = strsplit(i, "|", fixed = TRUE)
#     acc_only = append(acc_only, str[[1]][2])
#   }
#   return(acc_only)
# }
## TITLE: ----------------------------------------------------------------------
## DESCRIPTION:
##
## PARAMETERS: -
##
## OUTPUT:
##
getProteinAcc_uniqemappers <- function(list) {
  acc_only <- c()
  for (i in list) {
    str <- strsplit(i, "|", fixed = TRUE)
    acc_only <- append(acc_only, str[[1]][2])
  }
  return(acc_only)
}
