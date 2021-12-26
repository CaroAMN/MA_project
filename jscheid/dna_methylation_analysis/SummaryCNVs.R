#-----------------------------------------
#     Purpose
#-----------------------------------------
# This script aims to access all .igv files containing the microarray
# intensity data as the ratio methylated / (methylated + unmethylated)
# to observe copy-number-variations (CNVs) on a chromosomic scale
# These ratios are summed up for each chromosom position and divided
# by the number of patients * 100 to obtain the rate of alterations [%]
# Source: https://link.springer.com/article/10.1007%2Fs00401-018-1879-y
# Figure 2

library(ggplot2)
library(ggHoriPlot)
library(dplyr)
library(ggthemes)
library(tidyr)
library(stringr)
library(pheatmap)
library("RColorBrewer")

setwd("~/Documents/Uni/Master/Thesis/Methylation")

#-----------------------------------------
#     Input
#-----------------------------------------
# Access all .igv files
files <- list.files(pattern = "\\.igv$", recursive = T)

# Create dataframe containing the base position and the chromosome
# information
df <- read.csv(files[1], sep = '\t')

# Iterate over each file and append the signal column of each sample 
# to the dataframe
for (i in 2:length(files)) {
  f <- read.csv(files[i], sep = '\t')
  col <- unlist(strsplit(files[i], split = '/'))[4]
  df <- cbind(df, f[,5])
  colnames(df)[i+4] <- col
}

# Throw out values for females on chromY
sex <- read.csv("sample_sex.tsv", sep = '\t')
females <- sex$sampleID[which(sex$sex == "f")]
df[which(df$Chromosome == "chrY"),
   unlist(lapply(females, function(f){ grep(f, colnames(df))}))] <- NA

# Create df for primary and recurrent GBM
prim <- cbind(df[,1:4], df[,grep("prim", colnames(df))])
rec <- cbind(df[,1:4], df[,grep("rec", colnames(df))])

# Take median over all samples
prim$summary <- apply(prim[,5:ncol(prim)], MARGIN = 1, FUN = median, na.rm=T)
rec$summary <- apply(rec[,5:ncol(rec)], MARGIN = 1, FUN = median, na.rm=T)


# Change the ticks in such a way that they are in the middle of each chromosome
ticks <- function(v) {
  new_v = c()
  for (i in 1:length(v)) {
    if (i == 1) {
      new_v = c(new_v, v[i] / 2)
    } else {
      new_v = c(new_v, v[i-1] + ((v[i] - v[i-1]) / 2) )
    }
  }
  return(new_v)
}

summary_CNV = function(df, marker){
  # Determine label params for scatterplot 
  labels <- as.data.frame(df %>% group_by(Chromosome) %>% summarise(n=n(), median=median(summary)))
  labels <- labels[match(unique(df$Chromosome), labels$Chromosome), ]
  labels$Chromosome <- gsub("chr","" ,labels$Chromosome)
  labels$n <- cumsum(labels$n)
  labels$ticks <- ticks(labels$n)
  
  # Calculate transparency values: Points close to 0 are less visible than points far away from 0
  alpha <- (abs(df$summary) - min(abs(df$summary))) / (max(abs(df$summary)) - min(abs(df$summary)))
  alpha <- sqrt(alpha)
  red <- brewer.pal('RdBu', n = 4)[1]
  blue <- brewer.pal('RdBu', n = 4)[4]
  grey <- brewer.pal('RdGy', n = 4)[4]
  
  # Save median per chromosome in a df
  chr_med <- aggregate(df$summary, list(df$Chromosome), median) %>%
                mutate(n = aggregate(df$summary, list(df$Chromosome), length)[,2])
  chr_med <- chr_med[match(unique(df$Chromosome), chr_med$Group.1), ]
  chr_med <- rep.int(chr_med$x, chr_med$n)
  
  # Get marker indices for marker
  marker$index <- unlist(sapply(marker$start, function(s){
    for (i in 1:nrow(df)) {
      marker_chr <- marker$chr[which(marker$start == s)]
      if (marker_chr == df$Chromosome[i] && between(s, df$Start[i], df$End[i])) {
        return(i)
      }
    }
  }))
  
  # Get marker values
  marker$value <- df$summary[marker$index]
  
  # Scatterplot of primary data
  p <- ggplot() +
    geom_point(aes(x=seq(nrow(df)), y = df$summary), col = ifelse(df$summary > 0, blue, red),
               size=0.5, alpha = alpha, na.rm = T) +
    geom_point(aes(x=seq(nrow(df)), y = chr_med), col = grey, size = 0.5, alpha = 0.03) +
    scale_x_continuous(name="Chromosomes", breaks=labels$ticks, labels=labels$Chromosome, expand=c(0,0)) +
    geom_vline(xintercept = labels$n, alpha = 0.4) +
    ylab(expression("log"[2]*" copy-number ratio")) + 
    ylim(c(-0.4,0.4)) +
    theme_classic() +
    geom_point(aes(x = marker$index, y = marker$value), col = "black", size=0.8) +
    geom_text(aes(x = marker$index, y = marker$value, label=marker$name), size=3, hjust=0,
              vjust=0, fontface='bold', angle = 90)
  
  return(p)
}

# Include common markers
marker <- read.csv("common_markers.txt", sep = '\t')

summary_CNV(prim, marker)
summary_CNV(rec, marker) 


###########################################################################################
# Test

# Observation: Lots of noise in there
# -> Bin data and average out noise to observe true changes in chromosome gain or loss

# Bin every 10 signals together by taking the mean
prim_binned <- as.data.frame(apply(prim[,5:(ncol(prim)-1)], MARGIN = 2, function(c){
  suppressWarnings(return(colMeans(matrix(c, 10))))
}))

# Add the metadata
prim_binned <- cbind(prim[seq(1, nrow(prim), 10), 1:4], prim_binned)
# Reset index
rownames(prim_binned) <- NULL
prim_binned$summary <- apply(prim_binned[, 5:ncol(prim_binned)], MARGIN = 1, FUN = median, na.rm=T)

# Define params for new x lables
labels_binned <-as.data.frame(prim %>% group_by(Chromosome) %>% summarise(n=n(), median=median(summary)))
labels_binned <- labels_binned[match(unique(prim_binned$Chromosome), labels_binned$Chromosome), ]
labels_binned$Chromosome <- gsub("chr","" ,labels_binned$Chromosome)
labels_binned$n <- cumsum(labels_binned$n) / 10
labels_binned$ticks <- ticks(labels_binned$n)

ggplot(data=prim_binned, aes(x=seq(nrow(prim_binned)), y= summary)) +
  geom_point() + theme_classic() +
  scale_x_continuous(name="name", breaks=labels_binned$n, labels=labels_binned$Chromosome)


# Observation: Noise could be reduced. Chromosome loss and gain can still be observed
# But binning is not needed here
###########################################################################################

#-----------------------------------------
#     Horizonplot
#-----------------------------------------

# Alter the df to "long format"
df_copy <- df
df_copy$ID <- as.double(rownames(df_copy))
df_long <- df_copy %>% gather(key = "SampleID", value = "Signal",-ID, -Chromosome, -Start, -End, -Feature)
rm(df_copy)
# Define the cutpoints of the horizonplots
cutpoints = seq(-0.5, 0.45, length.out = 6)

# TODO: Rearrange horizon plot 
# df_long %>%
#   mutate(ordr_int = strtoi(gsub("GBM(.+)_.*","\\1", df_long$SampleID))) %>%
#   mutate(ordr_str = df_long$SampleID) %>%
#   arrange(ordr_int, ordr_str)) %>%
#   mutate(SampleID = factor(ordr_str, levels = unique(ordr_str))

h <- df_long %>% ggplot() +
  geom_horizon(aes(x=ID,
                   y = Signal,
                   fill=..Cutpoints..),
               origin = 0, horizonscale = cutpoints) +
  scale_fill_hcl(palette = 'RdBu', reverse = F) +
  facet_grid(SampleID~.) +
  theme_few() +
  theme(
    panel.spacing.y=unit(0, "lines"),
    strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank()) +
  guides(fill=guide_legend(title=expression(atop("log"[2]*" copy-","number ratio")),
                           title.theme = element_text(size = 11))) +
  scale_x_continuous(name="Chromosomes",
                     breaks=labels$ticks,
                     labels=labels$Chromosome,
                     expand=c(0,0)) +
  geom_vline(xintercept = labels$n, alpha = 0.5)
h
# -> Legend labels are a little bit buggy. Adjusted manually

#-----------------------------------------
#     Sample distances
#-----------------------------------------
SampleDist <- dist(as.matrix(t(na.omit(df[,5:ncol(df)]))))

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(as.matrix(SampleDist),
         clustering_distance_rows=SampleDist,
         clustering_distance_cols=SampleDist,
         col=colors)


heatmap(cor(na.omit(t(df[,5:ncol(df)]))))
heatmap(cor(na.omit(prim[,5:(ncol(prim)-1)])))
heatmap(cor(na.omit(rec[,5:(ncol(rec)-1)])))

dd <- dist(na.omit(t(prim[,5:(ncol(prim)-1)])), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc)
str(hc)
dd <- dist(t(df_dat), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc)

df_dat<- df[,5:(ncol(df))]
df_dat <- df_dat[complete.cases(df_dat),]

# Add metadata
# primary ordering gender
# secondary ordering samplenames


