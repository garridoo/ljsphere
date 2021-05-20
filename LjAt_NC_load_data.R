
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
source("cpcoa.func.R")
source("paths.R")

# files

design.file <- paste(data.dir, "LjAt_NC_design.txt", sep="")
otu_table.file <- paste(results.dir, "LjAt_NC_ASV_table.txt", sep="")
taxonomy.file <- paste(results.dir, "LjAt_NC_taxonomy.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, check.names=F)

# treat the order Burkholderiales (formerly Proteobacteriales) as Proteobacteria (Class)

taxonomy$Class <- as.character(taxonomy$Class)
idx <- taxonomy$Order=="Burkholderiales"
taxonomy$Class[idx] <- "Betaproteobacteria"

# rarefy ASV table

otu_table_raref <- rrarefy(otu_table, sample=1000)

# re-order data matrices

idx <- design$Original.SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$Original.SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]
otu_table_raref <- otu_table_raref[, idx]

colnames(otu_table) <- design$SampleID

# remove chloroplast reads

chloro <- taxonomy$ASV[taxonomy$Phylum=="Cyanobacteria/Chloroplast"]
otu_table <- otu_table[!rownames(otu_table) %in% chloro, ]

# normalize otu table

design$depth <- colSums(otu_table)
otu_table <- apply(otu_table, 2, function(x) x/sum(x))

