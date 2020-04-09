
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# process independently a given experiment (run)

# for invasion experiment (exp. F, run LjAt_004), only the competition
# treatment (design$treatment1=="SC32") should be consdiered (subset
# samples of interest in LjAt_SC_load_data.R script)

run <- "LjAt_006"

source("LjAt_SC_load_data.R")

idx <- design$genotype %in% c("col0", "gifu") &
       TRUE

design <- design[idx, ]
otu_table <- otu_table[, idx]
otu_table_unfiltered <- otu_table_unfiltered[, idx]


df <- melt(otu_table)
colnames(df) <- c("strain", "SampleID", "RA")
df$host <- design$host[match(df$SampleID, design$SampleID)]
df$family <- taxonomy$family[match(df$strain, taxonomy$strain)]
df <- df[df$host %in% c("Lj", "At"), ]

df$RA[df$RA < 0.001] <- 0.001

for (strain in unique(df$strain)) {

    idx <- taxonomy$strain==strain
    host <- as.character(taxonomy$host[idx])
    family <- as.character(taxonomy$family[idx])

    ra_host <- df$RA[df$strain==strain & df$host==host]
    ra_nonhost <- df$RA[df$strain==strain & df$host!=host]
    p <- wilcox.test(ra_host, ra_nonhost, alternative="greater")$p.value

    mean_ra_host <- mean(ra_host)
    mean_ra_nonhost <- mean(ra_nonhost)
    taxonomy$mean_ra_host[idx] <- mean_ra_host
    taxonomy$mean_ra_nonhost[idx] <- mean_ra_nonhost
    taxonomy$p[idx] <- p
    taxonomy$hpi_ratio[idx] <- mean_ra_host/mean_ra_nonhost

}

taxonomy$p_adjust <- p.adjust(taxonomy$p, method="fdr")
taxonomy$sig <- taxonomy$p_adjust < 0.05

print(max(taxonomy$hpi_ratio, na.rm=T))

write.table(taxonomy, paste(results.dir, "hpi_", run, ".txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

lim <- 9
taxonomy$hpi_ratio[taxonomy$hpi_ratio <= 1] <- 0
taxonomy$hpi_ratio[taxonomy$hpi_ratio >= lim] <- lim

p1 <- ggplot(taxonomy, aes(x=host, y=family)) +
             geom_point(aes(size=log2(mean_ra_host), fill=hpi_ratio, color=sig), shape=21) +
             labs(x="", y="") +
             scale_fill_gradientn(colors=c("darkgrey", "yellow", "red"), values=rescale(c(0, 1, lim)), limits=c(0, lim)) +
             scale_color_manual(values=c("transparent", "black")) +
             theme(axis.text.x=element_text(size=7.5)) +
             theme(axis.text.y=element_text(size=7.5)) +
             theme(axis.title=element_text(size=9)) +
             theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
             main_theme +
             theme(legend.position="left")
 
ggsave(paste(figures.dir, run, "_hpi_dotplot.pdf", sep=""), p1, width=10, height=5)
    
