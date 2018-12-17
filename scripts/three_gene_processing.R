library(tidyverse)
library(cowplot)
pinetree_counts <- read.table("three_genes_counts.tsv", sep = "\t", header = TRUE)
pinetree_tbl <- as.tibble(pinetree_counts)
pinetree_plot <- pinetree_tbl %>% filter(species %in% c("proteinX","rnapol","__ribosome")) 
gene_plot <- ggplot(data = pinetree_plot, aes(x=time,xlim = 100, y=protein,group=species)) + geom_line(aes(color=species))
save_plot("three_genes_example.png", gene_plot, base_aspect_ratio = 1.3)


iris %>% group_by(Species) %>%
  summarize(Sepal.Length.mean = mean(Sepal.Length))
