#!/usr/bin/env Rscript
library(ggplot2)
library(ggtree)
library(tidyverse)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)

pdf(NULL)

# Functions
get_eigenvector_variance <- function(int, pca.obj){
  eig = pca.obj$eig[int]
  eig.sum = sum(pca.obj$eig)
  eig.var = (eig / eig.sum)*100
  return(eig.var)
}

# Load command line args
## Arg1 = Path to input vcf. Compression allowed.
## Arg2 = Output prefix
args = commandArgs(trailingOnly=TRUE)

design <- read_tsv(
  '~/Interface/projects/hilda-popgen/data/vcf/design3.tsv',
  col_names = F
) %>% rename(sample = X1, population = X2)

# Perform phylogenetic analysis with poppr (https://grunwaldlab.github.io/poppr/)
## Followed tutorial here -> https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
vcf <- read.vcfR(
  args[1]
)

vcf.gl <- vcfR2genlight(vcf)
ploidy(vcf.gl) <- 2
pop(vcf.gl) <- design$population

tree <- aboot(vcf.gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = F)

# Perform pca
pca <- glPca(vcf.gl, nf = 3, center = T, scale = T)

# Plot pca
pca.scores <- as.data.frame(pca$scores)
pca.scores$pop <- pop(vcf.gl)
pc1.eig <- paste(
  as.character(round(get_eigenvector_variance(1, pca), 1)), '%', sep = ''
)
pc2.eig <- paste(
  as.character(round(get_eigenvector_variance(2, pca), 1)), '%', sep = ''
)

set.seed(9)
cols <- brewer.pal(n = nPop(vcf.gl), name = "Dark2")
ggplot(pca.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols)  +
  geom_hline(yintercept = 0)  +
  geom_vline(xintercept = 0)  +
  labs(
    x = paste('PC1', pc1.eig),
    y = paste('PC2', pc2.eig)
  ) +
  theme_bw()

# Save pca as png
ggsave(
    paste(args[2], 'pca.png', sep = '.'),
    dpi = 150,
    units = 'cm',
    width = 30,
    height = 20,
    device = 'png'
)

# Plot poppr output as dendrogram
population.labels <- tree$tip.label %>%
  as_tibble() %>%
  rename(sample = value) %>%
  inner_join(design)

group.info <- as.list(population.labels$sample)
names(group.info) <- population.labels$population
tree.grouped <- groupOTU(tree, group.info)

p <- ggtree(tree.grouped, aes(color=group), layout = 'rectangular')
p$data <- p$data %>% full_join(design, by = c('label' = 'sample'))
p + geom_tippoint(aes(color=population), size=2) + geom_tiplab(show.legend = FALSE)

# Save dendrogram as png
ggsave(
    paste(args[2], 'poppr-dend.png', sep = '.'),
    dpi = 150,
    units = 'cm',
    width = 30,
    height = 20,
    device = 'png'
)