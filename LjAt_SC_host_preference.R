
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# process independently a given experiment (run)

run <- "LjAt_010"

source("LjAt_SC_load_data.R")

at_strains_ra <- colSums(otu_table[rownames(otu_table) %in% at_strains, ])
lj_strains_ra <- colSums(otu_table[rownames(otu_table) %in% lj_strains, ])
otu_table <- rbind(otu_table, at_strains_ra, lj_strains_ra)

df <- melt(otu_table)
colnames(df) <- c("strain", "SampleID", "RA")
df$genotype <- design$genotype[match(df$SampleID, design$SampleID)]

df <- df[!df$genotype=="input", ]

idx <- match(df$SampleID, design$SampleID)
df$treatment <- design$genotype[idx]

df$family <- taxonomy$family[match(df$strain, taxonomy$strain)]
df$family <- as.character(df$family)
df$family[is.na(df$family)] <- "All strains"

df$strain <- as.character(df$strain)
df$strain[df$strain=="at_strains_ra"] <- "At-SPHERE strains"
df$strain[df$strain=="lj_strains_ra"] <- "Lj-SPHERE strains"

df$family_strain <- paste(df$family, " (", df$strain, ")", sep="")

df$RA_log <- log10(df$RA * 1000 + 1)

colors <- data.frame(group=c("col0",    "gifu",    "nfr5",    "soil",    "deps",    "cyp79b2b3", "Atbbc",    "Atfls2",  "Ljfls2", "MN47", "Lc", "wood"),
                    color=c("#f8756b", "#00b8e3", "#1d5966", "#654321", "#d48079", "#845753",   "#845753",  "#d48079", "#8aaeb6", al_color, lc_color, "brown"))

fmt_dcimals <- function(decimals=0) {
    function(x) format(x, nsmall=decimals, scientific=F)
}

family <- "All strains"

df_family <- df[which(df$family==family), ]

if (family!="All strains") df_family$strain <- factor(df_family$strain, levels=rev(sort(unique(df_family$strain))))

colors <- colors[colors$group %in% df_family$treatment, ]
df_family$treatment <- factor(df_family$treatment, levels=colors$group)
df_family <- unique(df_family)

df_family_at <- df_family[df_family$strain %in% c("At-SPHERE strains", at_strains), ]
df_family_lj <- df_family[df_family$strain %in% c("Lj-SPHERE strains", lj_strains), ] 

if (dim(df_family_at)[1] > 0 & dim(df_family_lj)[1] > 0) {

    # Kruskal-Wallis test of group differences

    pval <- kruskal.test(RA ~ genotype, data=df_family_at)$p.value

    lim_padding <- 0.2

    idx <- df_family_at$treatment %in% c("col0", "gifu", "Atbbc", "Atfls2", "Ljfls2")
    lim <- mean(df_family_at$RA[idx])
    h_lim <- min(1, lim+lim_padding)
    l_lim <- max(0, lim-lim_padding)
    # h_lim <- 1
    # l_lim <- 0

    p1 <- ggplot(df_family_at, aes(x=strain, y=RA, color=treatment, fill="transparent")) +
                geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=.2*length(unique(df_family$treatment)), fill="transparent") +
                geom_jitter(position=position_jitterdodge(jitter.width=.1*length(unique(df_family$treatment)),
                                                          dodge.width=.2*length(unique(df_family$treatment))), size=.3, alpha=0.7) +
                scale_y_continuous(position="left", labels=percent, limits=c(l_lim, h_lim), breaks=seq(0, 1, by=.1)) +
                scale_colour_manual(values=as.character(colors$color)) +
                labs(x="", y="Relative Abundance") +
                theme(axis.text.x=element_text(size=7.5)) +
                theme(axis.text.y=element_text(size=7.5)) +
                theme(axis.title=element_text(size=9)) +
                theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
                main_theme +
                theme(legend.position="none")

    idx <- df_family_lj$treatment %in% c("col0", "gifu", "Atbbc", "Atfls2", "Ljfls2")
    lim <- mean(df_family_lj$RA[idx])
    h_lim <- min(1, lim+lim_padding)
    l_lim <- max(0, lim-lim_padding)
    # h_lim <- 1
    # l_lim <- 0
  
    p2 <- ggplot(df_family_lj, aes(x=strain, y=RA, color=treatment, fill="transparent")) +
                geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=.2*length(unique(df_family$treatment)), fill="transparent") +
                geom_jitter(position=position_jitterdodge(jitter.width=.1*length(unique(df_family$treatment)),
                                                          dodge.width=.2*length(unique(df_family$treatment))), size=.3, alpha=0.7) +
                scale_y_continuous(position="right", labels=percent, limits=c(l_lim, h_lim), breaks=seq(0, 1, by=0.1)) +
                scale_colour_manual(values=as.character(colors$color)) +
                labs(x="", y="Relative Abundance") +
                theme(axis.text.x=element_text(size=7.5)) +
                theme(axis.text.y=element_text(size=7.5)) +
                theme(axis.title=element_text(size=9)) +
                theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
                main_theme +
                theme(legend.position="none")
    
    gA <- ggplotGrob(p1)
    gB <- ggplotGrob(p2)
    maxHeight = grid::unit.pmax(gA$heights, gB$heights)
    gA$heights <- as.list(maxHeight)
    gB$heights <- as.list(maxHeight)
    pg1 <- grid.arrange(gA, gB, ncol=2, top=paste(run, "; P=", format(pval, digits=2), sep=""))
    
    ggsave(paste(figures.dir, run, "_host_preference.pdf", sep=""), pg1, width=4, height=4)

}

