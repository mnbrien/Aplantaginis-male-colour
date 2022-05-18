#### TWISST ####
# Brien, Orteu et al. 2022

library(tidyverse)
library(data.table)
library(ape)

# Load twisst plotting functions. available from https://github.com/simonhmartin/twisst
source("plot_twisst.R")

## 1) Read in files

# Read in weights file
tgmoth.weights <- read.table("wwmapped_scaff419.weights.csv", header=T)
tgmoth.weights.na.number <- tgmoth.weights[!grepl("NaN", tgmoth.weights$topo1),]
# Change weights so they add up to 1 (i.e. divide by sum)
tgmoth.weights <- tgmoth.weights / apply(tgmoth.weights, 1, sum)

# Read in windows file
tgmoth.windows <- read.table("wwmapped_scaff419.data.tsv", header=T)

## 2) Make additive windows

# Load fasta fai 
ref.scaff<-read.table("tgmoth.WW.fa.fai")
names(ref.scaff)[1:2]<-c("rscaff" , "end")

# Reorder by scaffold size
ref.scaff <- ref.scaff[order(-(ref.scaff$end)),]

# Make additive window start, mid and end
ref.scaff$add <- c(0 , cumsum(ref.scaff$end[-1]))

tgmoth.windows$add.mid <- tgmoth.windows$mid + ref.scaff[match(tgmoth.windows$scaffold , ref.scaff$rscaff) , "add"]
tgmoth.windows$add.start <- tgmoth.windows$start + ref.scaff[match(tgmoth.windows$scaffold , ref.scaff$rscaff) , "add"]
tgmoth.windows$add.end <- tgmoth.windows$end + ref.scaff[match(tgmoth.windows$scaffold , ref.scaff$rscaff) , "add"]

## 3) Filter out small scaffolds

# Filter out small scaffolds
ref.scaff.filt <- subset(ref.scaff , ref.scaff$end > 1000000)

# Subset scaffolds which correspond to filtered scaffolds
tgmoth.windows.filt <- subset(tgmoth.windows, scaffold  %in% ref.scaff.filt$rscaff)

# Add scaffold names column
tgmoth.weights$scaffold <- tgmoth.windows$scaffold
# Subset scaffolds corresponding to filtered scaffolds
tgmoth.weights.filt <- subset(tgmoth.weights, scaffold  %in% ref.scaff.filt$rscaff)
# Remove scaffold column
tgmoth.weights.filt$scaffold <- NULL

## Remove NA windows (when used positional windows) in weights and windows file
tgmoth.weights.na <- tgmoth.weights.filt[!grepl("NaN", tgmoth.weights.filt$topo1),]

tgmoth.windows.na <- tgmoth.windows.filt
tgmoth.windows.na$nas <- tgmoth.weights.filt$topo1
tgmoth.windows.na <- tgmoth.windows.na[!grepl("NaN", tgmoth.windows.na$nas),]
tgmoth.windows.na$nas <- NULL

## 4) Reorder tables in window size 

tgmoth.weights.sort <- tgmoth.weights.na[order((tgmoth.windows.na$add.mid)),]
tgmoth.windows.sort <- tgmoth.windows.na[order((tgmoth.windows.na$add.mid)),]

## 5) Plot

## Plot raw data in "stepped" style, with polygons stacked - full scaffold
par(mfrow = c(1,1), mar = c(4,4,1,1))

plot.weights(weights_dataframe=tgmoth.weights.sort, positions=tgmoth.windows.sort[,c("start","end")],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE, xaxt="n", xlab="Scaffold 419",xlim=c(0,7300000))
axis(side=1, at = seq (0, 7500000, 500000))

## zoomed in plot
plot.weights(weights_dataframe=tgmoth.weights.sort, positions=tgmoth.windows.sort[,c("start","end")],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE, xaxt="n", xlab="Scaffold 419",xlim=c(7058000,7160000))
axis(side=1, at = seq (7058000, 7200000, 25000))

##############################################################################################################################

