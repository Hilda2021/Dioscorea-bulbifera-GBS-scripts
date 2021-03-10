#!/usr/bin/env Rscript

library(ggplot2)
library(ggtree)
library(tidyverse)

pdf(NULL)

# Read command line arguments
## Arg1 = Path to treefile from iqtree
## Arg2 = Output image path
args = commandArgs(trailingOnly=TRUE)

# Load experiment design file
design <- read_tsv(
  '~/Interface/projects/hilda-popgen/data/vcf/design3.tsv',
  col_names = F
) %>% rename(sample = X1, population = X2)

# Plot iqtree output as dendrogram
nwk <- treeio::read.newick(args[1])

population.labels <- nwk$tip.label %>%
  as_tibble() %>%
  rename(sample = value) %>%
  inner_join(design)

group.info <- as.list(population.labels$sample)
names(group.info) <- population.labels$population
nwk.grouped <- groupOTU(nwk, group.info)

p <- ggtree(nwk.grouped, aes(color=group), layout = 'rectangular', branch.length = 'none')
p$data <- p$data %>% full_join(design, by = c('label' = 'sample'))
p + geom_tiplab(show.legend = FALSE) + geom_tippoint(aes(color=population), size=2)

# Save the dendrogram
ggsave(
  args[2],
  dpi = 150,
  units = 'cm',
  width = 30,
  height = 20,
  device = 'png'
)