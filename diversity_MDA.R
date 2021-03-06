
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("vegan")

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")
source("cpcoa.func.R")

# directories

results.dir <- "/project_folder/LjAt_SC/results/MDA10/"
data.dir <- "/project_folder/LjAt_SC/data/"
figures.dir <- "/project_folder/figures/"

# files

design.file <- paste(data.dir, "MDA10_design.txt", sep="")
otu_table.file <- paste(results.dir, "rarefied_otu_table.txt", sep="")
otu_table_unfiltered.file <- paste(results.dir, "otu_table_unfiltered.txt", sep="")
taxonomy.file <- paste(data.dir, "MDA10_taxonomy.txt", sep="")

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
    
design$depth <- colSums(otu_table)
design$depth_unfiltered <- colSums(otu_table_unfiltered)

# normalize otu tables

otu_table <- apply(otu_table, 2, function(x) x/sum(x))
otu_table_unfiltered <- apply(otu_table_unfiltered, 2, function(x) x/sum(x))

idx <- rownames(otu_table_unfiltered) %in% c("good", "noisy", "chimera", "other")
design <- cbind(design, t(otu_table_unfiltered[idx, ]))
design$exact <- 1-rowSums(t(otu_table_unfiltered[!idx, ]))

# subset samples of interest
 
idx <- design$exudate %in% c("Col0", "Gifu", "axGifu", "nfr5", "M9") &
       design$SampleID != "LA.Gifu.18" &
       TRUE

design <- design[idx, ]
otu_table <- otu_table[, idx]
otu_table_unfiltered <- otu_table_unfiltered[, idx]

### beta diversity

colors <- data.frame(group=c("Col0", "axGifu", "Gifu", "nfr5", "M9", "input"),
                     color=c(at_color, lj_color, lj_color, lj_mutant_color, soil_color, c_black))

shapes <- data.frame(group=c("exp0095", "exp0096", "exp0097"),
                     shape=c(root_shape, soil_shape, rhizosphere_shape))

# PCoA Bray-Curtis

bray_curtis <- vegdist(t(otu_table), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

colors_pcoa <- colors[colors$group %in% points$exudate, ]
points$exudate <- factor(points$exudate, levels=colors_pcoa$group)
shapes_pcoa <- shapes[shapes$group %in% points$experiment, ]
points$experiment <- factor(points$experiment, levels=shapes_pcoa$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=exudate)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     scale_colour_manual(values=as.character(colors_pcoa$color)) +
     scale_shape_manual(values=shapes_pcoa$shape) +
     geom_text(aes(label=SampleID)) +
     theme(axis.text.x=element_text(size=12)) +
     theme(axis.text.y=element_text(size=12)) +
     theme(axis.title=element_text(size=13)) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
     main_theme +
     theme(legend.position="top")

ggsave(paste(figures.dir, "MDA10_PCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### CPCoA Bray-Curtis of samples conditioned by technical factors

d <- design

sqrt_transform <- T

bray_curtis <- vegdist(t(otu_table), method="bray")

capscale.gen <- capscale(bray_curtis ~ exudate_species + Condition(experiment), data=d, add=F, sqrt.dist=sqrt_transform)

# ANOVA-like permutation analysis

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)
                                                    
# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

points_cpcoa <- capscale.gen$CCA$wa[, 1:2]
colnames(points_cpcoa) <- c("x", "y")
points_cpcoa <- cbind(points_cpcoa, d[match(rownames(points_cpcoa), d$SampleID), ])

# plot CPCo 1 and 2

colors_cpcoa <- colors[colors$group %in% points_cpcoa$exudate, ]
points_cpcoa$exudate <- factor(points_cpcoa$exudate, levels=colors_cpcoa$group)
shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$experiment, ]
points_cpcoa$experiment <- factor(points_cpcoa$experiment, levels=shapes_cpcoa$group)

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=exudate)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     theme(axis.text.x=element_text(size=12)) +
     theme(axis.text.y=element_text(size=12)) +
     theme(axis.title=element_text(size=13)) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "MDA10_CPCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

