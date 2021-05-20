
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# subset samples of interest

idx <- design$Study %in% c("wippel_2020") &
      design$Soil.Batch %in% c("CAS", "CAS10", "CAS11b", "CAS14") &
      design$Description %in% c("Day36") &
      design$Compartment %in% c("root", "soil", "rhizosphere") &
      design$Host.Species %in% c("arabidopsis_thaliana", "lotus_japonicus", "host_soil") &
      design$depth >= 1000 &
      T  
   
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]
otu_table_raref_subset <- otu_table_raref[, idx]

# add a new column to design with nodulation phenotype

nod_plus <- c("gifu", "ram1")
nod_minus <- c("nfr5", "symrk", "ccamk")

design_subset$nodulation[design_subset$Host.Species=="arabidopsis_thaliana"] <- "At"
design_subset$nodulation[design_subset$Host.Genotype %in% nod_plus] <- "Lj (Nod+)"
design_subset$nodulation[design_subset$Host.Genotype %in% nod_minus] <- "Lj (Nod-)"

### beta diversity

colors <- data.frame(group=c("host_soil", "arabidopsis_thaliana", "lotus_japonicus"),
                     color=c(soil_color, at_color, lj_color))

shapes <- data.frame(group=c("root", "soil", "rhizosphere", "phycosphere"),
                     shape=c(root_shape, soil_shape, rhizosphere_shape, 11))

# PCoA Bray-Curtis

bray_curtis <- vegdist(t(otu_table_subset), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_subset[match(rownames(points), design_subset$SampleID), ])

colors <- colors[colors$group %in% points$Host.Species, ]
points$Host.Species <- factor(points$Host.Species, levels=colors$group)

shapes <- shapes[shapes$group %in% points$Compartment, ]
points$Compartment <- factor(points$Compartment, levels=shapes$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Host.Species, shape=Compartment)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle("PCoA of Bray-Curtis distances") +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "AtLj_NC_PCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### shannon index

index <- design_subset
index$shannon <- diversity(otu_table_raref_subset, index="shannon", MARGIN=2)
index$Host.Species <- factor(index$Host.Species, levels=colors$group)

# reorder boxplots

index$Host.Species <- factor(index$Host.Species, levels=colors$group)

p <- ggplot(index, aes(x=Host.Species, y=shannon, color=Host.Species, shape=Compartment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(1.5), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="Host.Species", y="shannon index") +
            ggtitle("Shannon diversity") +
            main_theme

ggsave(paste(figures.dir, "LjAt_NC_shannon.pdf", sep=""), p, width=shannon_width, height=shannon_height)

