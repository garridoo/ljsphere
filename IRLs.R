
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load data

source("IRLs_load_data.R")

# thresholding input OTU tables

otu_table_norm_At_NC <- apply(otu_table_At_NC, 2, function(x) x/sum(x))
otu_table_norm_Lj_NC <- apply(otu_table_Lj_NC, 2, function(x) x/sum(x))

# thresholding IRL OTU tables

min_reads <- 10

idx <- rowSums(otu_table_AtIRL) > min_reads
otu_table_AtIRL <- otu_table_AtIRL[idx, ]
otu_table_norm_AtIRL <- otu_table_norm_AtIRL[idx, ]

min_reads <- 10

idx <- rowSums(otu_table_LjIRL) > min_reads
otu_table_LjIRL <- otu_table_LjIRL[idx, ]
otu_table_norm_LjIRL <- otu_table_norm_LjIRL[idx, ]

# remove wells without a insufficient number of reads

min_reads <- 100

idx <- colSums(otu_table_AtIRL) > min_reads
otu_table_AtIRL <- otu_table_AtIRL[, idx]
otu_table_norm_AtIRL <- otu_table_norm_AtIRL[, idx]

min_reads <- 100

idx <- colSums(otu_table_LjIRL) > min_reads
otu_table_LjIRL <- otu_table_LjIRL[, idx]
otu_table_norm_LjIRL <- otu_table_norm_LjIRL[, idx]

# calculate recovery rates for culture-independent samples

source("IRLs_recovery_rates.R")

