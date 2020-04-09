
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load data

source("AtLj_NC_load_data.R")

# plot PCoA of Bray-Curtis

source("AtLj_NC_diversity.R")

# rank-abundance plots

source("AtLj_NC_rank.R")

