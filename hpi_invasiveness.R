
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load data from the invasion experiment (exp. F, run LjAt_004)

run <- "LjAt_004"

source("LjAt_SC_load_data.R")

full <- read.table(paste(results.dir, "hpi_004_fullSC.txt", sep=""), header=T, sep="\t")
atlj <- read.table(paste(results.dir, "hpi_004_At_Lj.txt", sep=""), header=T, sep="\t")
ljat <- read.table(paste(results.dir, "hpi_004_Lj_At.txt", sep=""), header=T, sep="\t")

idx <- full$host=="At"
at_strains <- data.frame(full[idx, c(1, 2, 3, 5, 8)], invasion=atlj$hpi_ratio[idx])

idx <- full$host=="Lj"
lj_strains <- data.frame(full[idx, c(1, 2, 3, 5, 8)], invasion=ljat$hpi_ratio[idx])

at_strains <- at_strains[rowSums(is.na(at_strains))==0, ]
lj_strains <- lj_strains[rowSums(is.na(lj_strains))==0, ]

df <- rbind(at_strains, lj_strains)
df <- df[df$mean_ra_host>0.001, ]

# correlation test (subset root or rhizosphere samples in LjAt_SC_load_data.R script)

test <- cor.test(df$hpi_ratio, df$invasion)

p1 <- ggplot(df, aes(x=invasion, y=hpi_ratio, color="black", fill=host)) +
             geom_point(aes(size=log2(mean_ra_host)), shape=21, color="black", alpha=0.8)+
             geom_smooth(method=lm, formula=y~x, color="black", size=0.5, fill=NA) +
             scale_y_log10(breaks=c(0, 0.5, 1, 5, 10), limits=c(NA, 10)) +
             scale_x_log10(breaks=c(0, 0.5, 1, 5, 10), limits=c(NA, 5)) +
             annotation_logticks(short=unit(1.5,"mm"), mid=unit(2.5,"mm"), long=unit(3.5,"mm")) +
             labs(x="Invasiveness", y="Host preference") +
             main_theme +
             theme(axis.text.x=element_text(size=7.5)) +
             theme(axis.text.y=element_text(size=7.5)) +
             theme(axis.title=element_text(size=9)) +
             theme(plot.margin=unit(c(5, 5, 5, 5), "mm")) +
             ggtitle(paste("r=", format(test$estimate, digits=2), "; P=", format(test$p.value, digits=2), sep="")) +
             theme(legend.position="none", plot.title=element_text(size=12))
 
ggsave(paste(figures.dir, "hpi_invasiveness.pdf", sep=""), p1, width=6, height=5)
 
