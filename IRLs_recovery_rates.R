
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# function to calculate recovery rates for abundant OTUs

recovery <- function(otu_table_norm_IRL, otu_table_norm_NC, host_IRL, host_NC, top=100, threshold=0.1) {

    # get important OTUs
    
    mean_ra <- apply(otu_table_norm_NC, 1, mean)
    
    idx <- mean_ra * 100 > threshold
    abundant_otus <- rownames(otu_table_norm_NC)[idx]
    abundant_otus <- abundant_otus[!is.na(abundant_otus)]
    
    idx <- sort(mean_ra, decreasing=T, index.return=T)$ix
    mean_ra <- mean_ra[idx]
    top_otus <- names(mean_ra)[1:top]
    top_otus <- top_otus[!is.na(top_otus)]

    # print recovery rates
    
    rc_top_abundant <- sum(abundant_otus %in% rownames(otu_table_norm_IRL)) * 100 / length(abundant_otus)
    rc_top <- sum(top_otus %in% rownames(otu_table_norm_IRL)) * 100 / length(top_otus)
    
    cat(paste("Abundant OTUs (>", threshold, "% avg. RA; ", length(abundant_otus), " OTUs): ",
              format(rc_top_abundant, digits=4), "%", "\n", sep=""))
    cat(paste("Top ", top, " OTUs by average RA: ",
              format(rc_top, digits=4), "%", "\n\n", sep=""))

    # write list of important OTUs

    write.table(abundant_otus, paste(results.dir, host_NC, "_abundant_otus.txt", sep=""),
                col.names=F, row.names=F, quote=F)
    write.table(top_otus, paste(results.dir, host_NC, "_top_otus.txt", sep=""),
                col.names=F, row.names=F, quote=F)

    # write representative sequences of important OTUs

    writeXStringSet(rep_seqs[names(rep_seqs) %in% abundant_otus],
                    paste(results.dir, host_NC, "_abundant_otus.fasta", sep=""))
    writeXStringSet(rep_seqs[names(rep_seqs) %in% top_otus],
                    paste(results.dir, host_NC, "_top_otus.fasta", sep=""))
    writeXStringSet(rep_seqs[names(rep_seqs) %in% rownames(otu_table_norm_IRL)],
                    paste(results.dir, host_IRL, "_IRL_otus.fasta", sep=""))

    mean_ra <- c(0, mean_ra)
    df <- data.frame(NC_OTU=names(mean_ra))
    df$recovered <- names(mean_ra) %in% rownames(otu_table_norm_IRL)
    df$RA <- mean_ra
    df$cumsum <- cumsum(df$RA * df$recovered)
    df$index <- 1:length(mean_ra)

    df_top <- df[1:top, ]
    df_top$NC_OTU <- factor(df_top$NC_OTU, levels=df_top$NC_OTU)
   
    # plot barchart of recovery for the top OTUs
    
    p <- ggplot(df_top, aes(x=NC_OTU, y=RA, fill=recovered)) +
         geom_bar(stat="identity", color="black", size=size_rec_barplot) +
         scale_fill_manual(values=c("transparent", "black")) +
         scale_y_continuous(labels=percent) +
         labs(x="Species found in natural communities", y="Relative abundance") +
         main_theme +
         theme(axis.ticks.x=element_blank(),
               axis.text.x=element_blank())
    
    ggsave(paste(figures.dir, "recovery_barplot_", host_IRL, host_NC, ".pdf", sep=""), p,
           width=width_rec_barplot, height=height_rec_barplot)

    # plot aggregated recovery rate

    if(host_IRL=="At") color <- at_color
    if(host_IRL=="Lj") color <- lj_color
    
    p <- ggplot(df_top, aes(x=index, y=cumsum)) +
         geom_line(stat="identity", color=color, size=size_cumsum) +
         scale_y_continuous(labels=percent, position="right") +
         geom_hline(yintercept=max(df$cumsum), color=color, size=size_cumsum) +
         labs(x="Species found in natural communities", y="Relative abundance") +
         main_theme +
         theme(axis.ticks.x=element_blank(),
               axis.text.x=element_blank())
    
    ggsave(paste(figures.dir, "recovery_cumsum_", host_IRL, host_NC, ".pdf", sep=""), p,
           width=width_rec_barplot, height=height_rec_barplot)

}

### full ICL input AL7 (both media)

cat("\nRecovery rates for AtIRL respect to At roots:\n\n")

recovery(otu_table_norm_AtIRL, otu_table_norm_At_NC, "At", "At")

cat("\nRecovery rates for AtIRL respect to Lj roots:\n\n")

recovery(otu_table_norm_AtIRL, otu_table_norm_Lj_NC, "At", "Lj")

cat("\nRecovery rates for LjIRL respect to Lj roots:\n\n")

recovery(otu_table_norm_LjIRL, otu_table_norm_Lj_NC, "Lj", "Lj")

cat("\nRecovery rates for LjIRL respect to At roots:\n\n")

recovery(otu_table_norm_LjIRL, otu_table_norm_At_NC, "Lj", "At")

