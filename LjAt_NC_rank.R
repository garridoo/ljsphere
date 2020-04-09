
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# subset samples of interest

idx <- design$Study %in% c("zgadzaj_2018", "masayoshi_2016") &
       design$Soil.Batch %in% c("CAS", "CAS10", "CAS11b") &
       design$Compartment %in% c("root") &
       T
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

idx <- match(rownames(otu_table_subset), taxonomy$ASV) 
taxonomy_subset <- taxonomy[idx, ]


colors <- data.frame(group=c("Soil", "At", "Lj (Nod+)", "Lj (Nod-)"),
                     color=c(soil_color, at_color, lj_color, lj_mutant_color))

shapes <- data.frame(group=c("root", "soil"),
                     shape=c(root_shape, soil_shape))

taxonomy$taxon <- as.character(taxonomy$Class)
taxonomy$taxon[is.na(taxonomy$taxon)] <- "Other"

taxon_table <- aggregate(otu_table_subset, by=list(taxonomy_subset$taxon), FUN=sum)

df <- melt(taxon_table)
colnames(df) <- c("taxon", "SampleID", "ra")
df$host <- design_subset$Host.Species[match(df$SampleID, design_subset$SampleID)]
df <- df[df$taxon!="Other", ]

df_stats <- aggregate(df$ra, by=list(df$taxon, df$host), FUN=mean)
colnames(df_stats) <- c("taxon", "host", "mean_ra")
df_stats$sd_ra <- aggregate(df$ra, by=list(df$taxon, df$host), FUN=sd)$x

abundance_threshold <- 0.001

df_stats_at <- df_stats[df_stats$host=="arabidopsis_thaliana", ]
df_stats_at <- df_stats_at[order(df_stats_at$mean_ra, decreasing=T), ]
top_taxa_at <- df_stats_at$taxon[df_stats_at$mean_ra > abundance_threshold]

df_stats_lj <- df_stats[df_stats$host=="lotus_japonicus", ]
df_stats_lj <- df_stats_lj[order(df_stats_lj$mean_ra, decreasing=T), ]
top_taxa_lj <- df_stats_lj$taxon[df_stats_lj$mean_ra > abundance_threshold]

top_taxa <- union(top_taxa_at, top_taxa_lj)

taxa <- unique(df$taxon)
wt <- data.frame(taxon=taxa, p=NA)

for (taxon in taxa) {

    taxon_ra_at <- df$ra[df$taxon==taxon & df$host=="arabidopsis_thaliana"]
    taxon_ra_lj <- df$ra[df$taxon==taxon & df$host=="lotus_japonicus"]
    wt$p[wt$taxon==taxon] <- wilcox.test(taxon_ra_at, taxon_ra_lj)$p.value

}

wt$p.adjusted <- p.adjust(wt$p)
wt <- wt[!is.nan(wt$p), ]
wt <- wt[wt$p.adjusted < 0.001, ]
sig_taxa <- as.character(wt$taxon[wt$taxon %in% top_taxa])

print(sig_taxa)

# rank abundance plot

df_stats <- df_stats[df_stats$taxon %in% top_taxa, ]

idx <- order(aggregate(df_stats$mean_ra, by=list(df_stats$taxon), FUN=mean)$x, decreasing=T)
df_stats$taxon <- factor(df_stats$taxon, levels=df_stats$taxon[idx])

p1 <- ggplot(df_stats, aes(x=taxon, y=mean_ra, fill=host)) +
             geom_bar(stat="identity", position=position_dodge(), color="black") +
             theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1)) +
             theme(axis.text.y=element_text(size=7.5)) +
             theme(axis.title=element_text(size=9)) +
             theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
             main_theme
 
    
ggsave(paste(figures.dir, "AtLj_NC_barplots_.pdf", sep=""), p1, width=8, height=5)
 
