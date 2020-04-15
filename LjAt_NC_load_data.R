
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# load libraries

library("ggplot2")
library("scales")
library("grid")
library("vegan")

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")

# directories

results.dir <- "/project_folder/LjAt_NC/results/"
data.dir <- "/project_folder/LjAt_NC/data/"
figures.dir <- "/project_folder/figures/"

# files

design.file <- paste(data.dir, "design.txt", sep="")
otu_table.file <- paste(results.dir, "ASV_table_rarefied_1000.txt", sep="")
taxonomy.file <- paste(results.dir, "ASV_taxonomy_rdp.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, check.names=F)

# re-order data matrices

idx <- design$Original.SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$Original.SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

colnames(otu_table) <- design$SampleID

# remove chloroplast reads

chloro <- taxonomy$ASV[taxonomy$Phylum=="Cyanobacteria/Chloroplast"]
otu_table <- otu_table[!rownames(otu_table) %in% chloro, ]

# normalize otu table

design$depth <- colSums(otu_table)
otu_table <- apply(otu_table, 2, function(x) x/sum(x))

