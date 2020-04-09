
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

pcoa_width <- 8
pcoa_height <- 6
pcoa_size <- 2
pcoa_alpha <- 0.7

shannon_width <- 4
shannon_height <- pcoa_height
shannon_alpha <- 0.7

boxplot_size <- 1
boxplot_width <- 0.75
boxplot_jitter_size <- 1

width_rec_barplot <- 5
height_rec_barplot <- 3
size_rec_barplot <- 0.35

size_cumsum <- 0.75

at_color <- "#f8766c"
lj_color <- "#00bfc4"
lj_mutant_color <- "#8aaeb6"
soil_color <- "#654321"

root_shape <- 19
rhizosphere_shape <- 3
soil_shape <- 18

# ggplot2 theme

main_theme <- theme(panel.background=element_blank(),
              panel.grid=element_blank(),
              axis.line=element_line(color="black", size=1),
              axis.ticks=element_line(color="black", size=1),
              legend.background=element_blank(),
              legend.key=element_blank(),
              text=element_text(size=18, color="black"),
              legend.position="none")

