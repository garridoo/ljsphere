
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

library(ggplot2, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)
library(Biostrings, quietly=T, warn.conflicts=F)
library(gridExtra, quietly=T, warn.conflicts=F)
library(MASS, quietly=T, warn.conflicts=F)
library(diptest, quietly=T, warn.conflicts=F)

kmer_chunks <- function(seqs, width, k){
   
    x <- unlist(seqs)
    kmers <- matrix(ncol=length(x)/width-1, nrow=k^k)
    start <- 1
    end <- 0
    
    for (i in 1:(length(x)/width-1)) {

        start <- end + 1
        end <- end + width

        if (end < length(x) & start < end)
            kmers[, i] <- oligonucleotideFrequency(x[start:end], width=k)

    }

    return(kmers)

}

args <- commandArgs(TRUE)

if(length(args) < 2) {

        cat("usage: assembly_stats.R assembly_file genome_id output_dir\n")

} else {
    
    # parse arguments

    assembly.file <- args[1]
    genome.id <- args[2]
    output.dir <- args[3]
    
    min.length <- 0
    width <- 10000

    # generate a vector of contig lengths

    seqs <- readDNAStringSet(assembly.file)
    seqs <- seqs[width(seqs) >= min.length]
    n.contigs <- sum(width(seqs) >= min.length)
    names(seqs) <- gsub("^", "contig_", 1:length(seqs))
    
    # calculate assembly statistics from contig lengths
    
    lengths.table <- sort(width(seqs))
    lengths <- data.frame(ctg_number=1:length(seqs), acc=cumsum(lengths.table))
    tot.length <- lengths[dim(lengths)[1], 2]
    n50.idx <- which(lengths[, 2] >= tot.length * .50)[1]
    n50 <- lengths.table[n50.idx]
    n90.idx <- which(lengths[, 2] >= tot.length * .10)[1]
    n90 <- lengths.table[n90.idx]

    # GC content spectrum

    window.size <- 200
    gc <- rowSums(letterFrequencyInSlidingView(unlist(seqs), 
                                               window.size, c("G","C"))) / window.size

    # kmer spectrum MDS

    if (sum(width(seqs)) < width)
        width <- sum(width(seqs)) / 10
    kmers <- kmer_chunks(seqs, width, k=4)
    d <- dist(t(kmers))
    d[d==0] <- 0.0001
    points <- isoMDS(d, k=2)$points
    points <- as.data.frame(points)
    colnames(points) <- c("x", "y")

    # Hartigan's dip test of multimodality
    
    alpha <- 0.03
    dtest_kmer <- dip(isoMDS(d, k=1)$points)
    gc_coarse <- rowSums(letterFrequencyInSlidingView(unlist(seqs),
                                                      width, c("G","C"))) / width
    dtest_gc <- dip(gc_coarse[seq(1, length(gc_coarse), window.size)])

    ### plotting

    main_theme <- theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(color="black"),
                        axis.ticks=element_line(color="black"),
                        axis.text=element_text(colour="black", size=10),
                        legend.position="top")

    p1 <- ggplot(lengths, aes(x=ctg_number, y=acc)) +
                 geom_point(cex=2, color="blue", alpha=.7) +
                 geom_line(color="blue") + 
                 geom_vline(xintercept=n50.idx, linetype=2, color="grey") +
                 annotate("text", x=n50.idx, y=1, hjust=0,
                          label=paste("N50 (", n50, " bp)", sep=""),
                          size=3, angle=90) +
                 geom_vline(xintercept=n90.idx, linetype=2, color="grey") +
                 annotate("text", x=n90.idx, y=1, hjust=0,
                          label=paste("N90 (", n90, " bp)", sep=""),
                          size=3, angle=90) +
                 labs(x="number of contigs (sorted by increasing length)",
                      y="accumulated total contig length") +
                 ggtitle(genome.id) +
                 main_theme 
   
    df <- as.data.frame(gc)

    p2 <- ggplot(data=df, aes(x=gc)) + 
                 geom_histogram(binwidth=.01, size=.5,
                                colour="red", fill="transparent") + 
                 xlim(0, 1) + 
                 labs(x="GC content (%)", y="") + 
                 theme(axis.ticks.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.line.y=element_blank()) + 
                 main_theme
   
    p3 <- ggplot(points, aes(x=x, y=y)) +
          geom_point(alpha=.6) +
          labs(x="tetranucleotide composition (MDS)", y=" ") + 
          main_theme

    pbottom <- arrangeGrob(p3, p2, nrow=1, ncol=2, heights=4, widths=c(4, 4))
    p <- arrangeGrob(p1, pbottom, ncol=1, heights=c(5, 4), widths=8)

    ggsave(filename=paste(output.dir, "/", genome.id, ".pdf", sep=""),
           plot=p, width=8, height=9)
    
    # output statistics in text file

    file <-  paste(output.dir, "/", genome.id, ".txt", sep="")
    sink(file, append=F)
    cat(paste(genome.id, length(seqs), tot.length, n50, n90,
              format(dtest_kmer, digits=2), format(dtest_gc, digits=2), sep="\t"))
    cat("\n")

}

