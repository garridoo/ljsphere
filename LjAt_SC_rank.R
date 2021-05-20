
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")
source("cpcoa.func.R")

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("vegan")

# process independently a given experiment (run)

run <- "LjAt_001"

# directories

results.dir <- paste("/netscratch/dep_psl/grp_rgo/garridoo/", run, "/", sep="")
data.dir <- "/biodata/dep_psl/grp_rgo/ljsphere/competition_experiments/data/"
figures.dir <- "/biodata/dep_psl/grp_rgo/ljsphere/figures/"

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

idx <- design$depth_unfiltered >= 1000
design <- design[idx, ]
otu_table <- otu_table[, idx]
otu_table_unfiltered <- otu_table_unfiltered[, idx]

otu_table <- apply(otu_table, 2, function(x) x/sum(x))
otu_table_unfiltered <- apply(otu_table_unfiltered, 2, function(x) x/sum(x))

idx <- rownames(otu_table_unfiltered) %in% c("good", "noisy", "chimera", "other")
design <- cbind(design, t(otu_table_unfiltered[idx, ]))
design$exact <- 1-rowSums(t(otu_table_unfiltered[!idx, ]))

# subset samples of interest
 
idx <- design$genotype %in% c("col0", "gifu", "nfr5", "soil", "deps", "cyp79b2b3", "Atbbc", "Atfls2", "Ljfls2", "MN47", "Lc", "wood") & 
       design$genotype %in% c("col0", "gifu", "soil") & 
       design$compartment %in% c("root", "droplet", "soil") &
       # design$treatment %in% c("fullSC") &
       !design$SampleID %in% c("w12.3") &
       TRUE

design <- design[idx, ]
otu_table <- otu_table[, idx]
otu_table_unfiltered <- otu_table_unfiltered[, idx]

### beta diversity

colors <- data.frame(group=c("col0",    "gifu",    "nfr5",    "soil",    "deps",    "cyp79b2b3", "Atbbc",    "Atfls2",  "Ljfls2", "MN47", "Lc", "wood"),
                     color=c("#f8756b", "#00b8e3", "#1d5966", "#654321", "#d48079", "#845753",   "#845753",  "#d48079", "#8aaeb6", al_color, lc_color, "brown"))

shapes <- data.frame(group=c("root", "soil", "rhizosphere", "deadroot", "wood"),
                     shape=c(root_shape, soil_shape, rhizosphere_shape, 16, 16))

# plot rank abundance boxplots per genotype

df <- melt(otu_table)
colnames(df) <- c("strain", "SampleID", "RA")
df$genotype <- design$genotype[match(df$SampleID, design$SampleID)]
df$host <- taxonomy$host[match(df$strain, taxonomy$strain)]
df$family <- taxonomy$family[match(df$strain, taxonomy$strain)]

taxonomy$p.val <- NULL
taxonomy$fc <- 1
taxonomy$mean_ra_lj <- taxonomy$mean_ra_at <- taxonomy$mean_ra_soil <- 0

df$RA[df$RA < 0.001] <- 0.001

for (strain in unique(df$strain)) {

     idx <- df$strain==strain

     a <- df$RA[idx & df$genotype=="gifu"]
     b <- df$RA[idx & df$genotype=="col0"]
     c <- df$RA[idx & df$genotype=="soil"]
   
    if (sum(a) > 0 & sum(b) > 0) {
        p.val <- wilcox.test(a, b)$p.val
    
        idx <- taxonomy$strain==strain

        taxonomy$p.val[idx] <- p.val
        taxonomy$fc[idx] <- mean(a) / mean(b)

        taxonomy$mean_ra_lj[idx] <- mean(a)
        taxonomy$mean_ra_at[idx] <- mean(b)
        taxonomy$mean_ra_soil[idx] <- mean(c)

    }
    
}

taxonomy$p.adj <- p.adjust(taxonomy$p.val, method="fdr")

# taxonomy <- taxonomy[!is.na(taxonomy$p.val), ]
taxonomy$differential <- taxonomy$p.adj < 0.05

df$fc <- df$RA / taxonomy$mean_ra_at[match(df$strain, taxonomy$strain)]
df$differential <- taxonomy$differential[match(df$strain, taxonomy$strain)]
df$differential[is.na(df$differential)] <- FALSE

# p1 <- ggplot(df, aes(x=SampleID, y=strain, fill=log(RA), color=differential)) +
#              geom_tile() +
#              # geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=.2*length(unique(df$genotype)), fill="transparent") +
#              # geom_jitter(position=position_jitterdodge(jitter.width=.1*length(unique(df$genotype)),
#              #                                           dodge.width=.2*length(unique(df$genotype))), size=.3, alpha=0.7) +
#              # # scale_y_continuous(position="right", labels=percent, limits=c(0, 1), breaks=seq(0, 1, by=.1)) +
#              # scale_colour_manual(values=as.character(colors$color)) +
#              # labs(x="", y="Relative Abundance") +
#              theme(axis.text.x=element_text(size=7.5)) +
#              theme(axis.text.y=element_text(size=7.5)) +
#              theme(axis.title=element_text(size=9)) +
#              theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
#              main_theme +
#              theme(legend.position="right")
# 
# ggsave(paste(figures.dir, run, "_fc.pdf", sep=""), p1, width=8, height=4)
# 
