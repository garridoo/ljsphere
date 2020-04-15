
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("gridExtra")
library("vegan")

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")

# process run LjAt_004 (invasion experiment)

run <- "LjAt_004"

# directories

results.dir <- paste("/project_folder/LjAt_SC/", run, "/", sep="")
data.dir <- "/project_folder/LjAt_SC/data/"
figures.dir <- "/project_folder/figures/"

# files

design.file <- paste(data.dir, run, "_design.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table.txt", sep="")
otu_table_unfiltered.file <- paste(results.dir, "otu_table_unfiltered.txt", sep="")
taxonomy.file <- paste(data.dir, run, "_taxonomy.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
otu_table_unfiltered <- read.table(otu_table_unfiltered.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, check.names=F)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- match(design$SampleID, colnames(otu_table_unfiltered))
otu_table_unfiltered <- otu_table_unfiltered[, idx]
   
# normalize otu tables

design$depth <- colSums(otu_table)
design$depth_unfiltered <- colSums(otu_table_unfiltered)

otu_table <- apply(otu_table, 2, function(x) x/sum(x))
otu_table_unfiltered <- apply(otu_table_unfiltered, 2, function(x) x/sum(x))

idx <- rownames(otu_table_unfiltered) %in% c("good", "noisy", "chimera", "other")
design <- cbind(design, t(otu_table_unfiltered[idx, ]))
design$exact <- 1-rowSums(t(otu_table_unfiltered[!idx, ]))

# aggregate RAs of strains from each host

lj_strains <- rownames(otu_table)[grepl("Lj", rownames(otu_table))]
design$lj_strains_ra <- colSums(otu_table[rownames(otu_table) %in% lj_strains, ])

at_strains <- rownames(otu_table)[!grepl("Lj", rownames(otu_table))]
design$at_strains_ra <- colSums(otu_table[rownames(otu_table) %in% at_strains, ])

# subset samples of interest
 
design$treatment <- paste(design$treatment1, design$treatment2)

idx <- design$genotype %in% c("col0", "gifu", "none") &
       design$compartment %in% c("root") &
       # design$compartment %in% c("rhizosphere") &
       TRUE

design <- design[idx, ]
otu_table <- otu_table[, idx]
otu_table_unfiltered <- otu_table_unfiltered[, idx]

