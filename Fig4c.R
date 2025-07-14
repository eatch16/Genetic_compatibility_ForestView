#----------------- Figure 4c --------------------------------------------------
# Load packages
library(data.table)
library(dplyr)
library(vegan)

# Plotting settings
fig_c_theme <- theme(
  axis.text.x = element_text(size = rel(1.3)),
  axis.text.y = element_text(size = rel(1.3)),
  axis.title = element_text(size = rel(1.3)),
  legend.title = element_text(size = rel(1.2)),
  legend.text = element_text(size = rel(1)),
  plot.title = element_text(size = rel(3), face = "bold"),
  plot.title.position = "plot"
)

phylum_palette <- c(
  "Actinomycetota"='#009E73',
  "Bacillota" ="#56B4E9",
  "Bacteroidota"='#E69F00',
  "Campylobacterota"="#CC79A7",
  "Pseudomonadota"="#D55E00",
  "Other"="grey"
)

# Analysis
kmers <- fread("D:/data_diver/Genetic_compatibility/data_for_figures/representative_kmer_distributions.txt") %>% as.data.frame()
taxonomy <- data.frame(fread("D:/data_diver/Genetic_compatibility/data_for_figures/assembly_full_lineage.txt", header=FALSE)) %>%
  select(c(3:5)) %>%
  na.omit() %>%
  unique() %>%
  `colnames<-`(c("phylum", "class", "order")) %>%
  `rownames<-`(.$order)

phylum <- taxonomy[kmers[, 1], 1]
rel_phyl <- c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Pseudomonadota")
phylum <- ifelse(phylum %in% rel_phyl, phylum, "Other")

data <- kmers[, -1]
data <- data[phylum != "Other", ]
phylum <- phylum[phylum != "Other"]

kmers_mds <- metaMDS(comm = data, distance = "euclidean", trace = TRUE, autotransform = FALSE)

# Complile results
MDS_xy <- data.frame(kmers_mds$points)
MDS_xy$phylum <- phylum

# Plot and save Fig 4c
fig_4c <- ggplot(MDS_xy, aes(x = MDS1, y = MDS2, color = phylum)) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = phylum_palette,
    name = "Phylum"
  ) +
  theme_minimal() +
  ggtitle("c") +
  fig_c_theme

pdf("fig_4c.pdf", height = 5, width=15)
plot(fig_4c)
dev.off()
