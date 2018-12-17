library(tidyverse)
library(cowplot)
pinetree_counts <- read.table("test_protein_counts.tsv", sep = "\t", header = TRUE)
pinetree_tbl <- as.tibble(pinetree_counts)
pinetree_plot <- pinetree_tbl %>% filter(species %in% c("proteinP","__ribosome")) 
gene_plot <- ggplot(data = pinetree_plot, aes(x=time, y=protein,group=species) + geom_line(aes(color=species)) + ggtitle("Protein Production Rate and Free Ribosomes in Half-Slow Transcript", subtitle = waiver()))
gene_plot
save_plot("half_slow_example.png", gene_plot, base_aspect_ratio = 1.3)


#ribo_slopes <- pinetree_tbl %>% filter(species %in% c("__ribosome"))
#model <- ribo_slopes[c("time")]
#free_ribosomes <- ribo_slopes["protein"]
#time <- ribo_slopes["time"]
#linear_model <- lm(time ~ free_ribosomes)