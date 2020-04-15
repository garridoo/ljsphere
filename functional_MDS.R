
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list = ls())

# load libraries

library(utils, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(MASS, quietly=T, warn.conflicts=F)
library(gridExtra, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)
library(vegan, quietly=T, warn.conflicts=F)

options(warn=-1)

# plotting stuff

source("plotting_functions.R")
source("plotting_parameters.R")

# paths to files

data.dir <- "/project_folder/LjAt_WGS/results/"
figures.dir <- "/project_folder/figures/"
mapping.file <- paste(data.dir, "/../mapping.txt", sep="")
annotation.dir <- paste(data.dir, "/annotations/", sep="")

# load data

message("loading matrix of functional profiles...")

mapping <- read.table(mapping.file, sep="\t", header=T, colClasses="character")

idx <- mapping$quality=="clean" & ((mapping$collection=="AtSPHERE" & mapping$compartment=="Root") | mapping$collection=="LjSPHERE")
mapping <- mapping[idx, ]

### functional profiles

message("generating matrix of functional profiles...")

ko.all <- data.frame(genome=NULL, ko=NULL) 
sizes.all <- data.frame(genome=NULL, size=NULL) 

i <- 1

for (g in mapping$ID) {
 
    i <- i + 1
   
    ko <- read.table(paste(annotation.dir, g, ".ko", sep=""),
                     fill=T, header=F, sep="\t",
                     col.names=c("peg", "ko", "gene", "description"),
                     colClasses="character",
                     quote="")[, 2]
    ko.genome <- data.frame(genome=g, ko=ko)
    ko.all <- rbind(ko.all, ko.genome)

    size.genome <- data.frame(genome=g, size=dim(ko.genome)[1])
    sizes.all <- rbind(sizes.all, size.genome)

}

ko.table <- table(ko.all)
ko.table <- t(ko.table[, -1])

sizes.all$perc_annotated <- colSums(ko.table) / sizes.all$size

func <- (ko.table > 0) * 1

# calculate pairwise functional distances

message("calculating pairwise functional distances...")

d <- 1 - cor(func)
diag(d) <- 0

### PCoA of functional distances

message("calculating functional PCoA...")

k <- 2

pcoa <- cmdscale(d, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig

# points <- tsne(d, k=2, perplexity=200)
# rownames(points) <- rownames(d)

points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, mapping[match(rownames(points), mapping$ID), ])

p1 <- ggplot(points, aes(x=x, y=y, color=host, shape=phylum)) +
      geom_point(alpha=.7, size=2) +
      scale_shape_manual(values=c(16, 17, 15, 3)) +
      labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
      main_theme +
      theme(legend.position="right")

ggsave(file=paste(figures.dir, "functional_MDS.pdf", sep=""), p1, height=6, width=8)

### functional overlap

labels <- colnames(func)
idx <- grepl("Lj", labels)
lj <- rownames(func)[rowSums(func[, idx])!=0]
at <- rownames(func)[rowSums(func[, !idx])!=0]
all <- rownames(func)[rowSums(func)!=0]

lj_specific <- sum(!lj %in% at)
at_specific <- sum(!at %in% lj)
shared <- sum(lj %in% at)
total <- lj_specific + at_specific + shared

nperm <- 1000
perm <- rep(NA, nperm)

perm <- data.frame(shared_p=NULL, lj_specific_p=NULL, at_specific_p=NULL)

for (i in 1:nperm) {

    perm_labels <- sample(colnames(func), dim(func)[2], replace=T)
    
    idx <- grepl("Lj", perm_labels)
    lj <- rownames(func)[rowSums(func[, idx])!=0]
    at <- rownames(func)[rowSums(func[, !idx])!=0]
   
    shared_Fp <- sum(lj %in% at)
 
    lj_specific_p <- sum(!lj %in% at)
    at_specific_p <- sum(!at %in% lj)
    shared_p <- sum(lj %in% at)
   
    perm <- rbind(perm, data.frame(shared_p, lj_specific_p, at_specific_p))

}

p <- (sum(perm$lj_specific_p > shared) + 1) / (nperm + 1)

print(paste("P = ", p, "; F = ", shared, sep=""))

### PERMANOVA of functional distances

adonis(d ~ mapping$host)
adonis(d ~ mapping$family)

