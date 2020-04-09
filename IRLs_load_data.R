
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# load libraries

library("ggplot2")
library("scales")
library("grid")
library("vegan")
library("Biostrings")

# load plotting functions

source("plotting_parameters.R")
source("plotting_functions.R")

# directories

results.dir <- "/project_folder/LjAt_IRL/results/"
data.dir <- "/project_folder/LjAt_IRL/data/"
figures.dir <- "/project_folder/figures/"

# files

design.file <- paste(data.dir, "design.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table.txt", sep="")
rep_seqs.file <- paste(results.dir, "rep_seqs.fasta", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
rep_seqs <- readDNAStringSet(rep_seqs.file)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# normalize otu tables

otu_table_norm <- apply(otu_table, 2, function(x) x/sum(x))

# split OTU tables

idx <- design$group %in% c("At_NC_root")
otu_table_At_NC <- otu_table[, idx]
otu_table_norm_At_NC <- otu_table_norm[, idx]

idx <- design$group %in% c("AtIRL", "At_root_CP")
otu_table_AtIRL <- otu_table[, idx]
otu_table_norm_AtIRL <- otu_table_norm[, idx]

idx <- design$group %in% c("Lj_NC_root")
otu_table_Lj_NC <- otu_table[, idx]
otu_table_norm_Lj_NC <- otu_table_norm[, idx]

idx <- design$group %in% c("LjIRL")
otu_table_LjIRL <- otu_table[, idx]
otu_table_norm_LjIRL <- otu_table_norm[, idx]

