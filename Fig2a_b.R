#----------------- Figure 2a,b ------------------------------------------------
# Load packages
library(readxl)
library(ggplot2)

# Plotting settings
fig_ab_theme <- theme(
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = rel(1.3)),
  axis.title = element_text(size = rel(1.3)),
  legend.title = element_text(size = rel(1.2)),
  legend.text = element_text(size = rel(1)),
  plot.title = element_text(size = rel(3), face = "bold")
)

mechanism_palette <- c(
  "AAC" = '#a50026', 
  "APH" = "#d73027", 
  "Class A, C, D" = '#fdae61',
  "Class B" = "#fee090", 
  "Erm" = "#ffffff", 
  "Mph" = "#e0f3f8",
  "Tet efflux" = "#74add1", 
  "Tet enzyme" = "#4575b4", 
  "Tet RPG" = "#313695",
  "Qnr" = "#000000"
)

# Load data
df <- read_excel("D:/data_diver/Genetic_compatibility/data_for_figures/summary_transferred_args.xlsx")
df <- as.data.frame(df)

# Analysis
aggregated_data <- data.frame(matrix(nrow=10, ncol=4))
colnames(aggregated_data) <- c("mechanism", "antibiotic_class", "predicted_ARGs", "HGT_Transfers")
aggregated_data$mechanism <- c("AAC", "APH", "Class A, C, D", "Class B", "Qnr", "Erm", "Mph", "Tet efflux", "Tet enzyme", "Tet RPG")
aggregated_data$antibiotic_class <- c(rep("Aminoglycoside", 2), rep("Beta-lactam", 2), "Fluoroquinolone", rep("Macrolide", 2), rep("Tetracycline", 3))

for (mech in unique(aggregated_data$mechanism)) {
  if (mech == "Class A, C, D") {
    aggregated_data[aggregated_data$mechanism == "Class A, C, D", 3:4] <- colSums(df[df$gene_class == "Class A" | df$gene_class == "Class C" | df$gene_class == "Class D1" | df$gene_class == "Class D2", 3:4])
  }
  else {
    aggregated_data[aggregated_data$mechanism == mech, 3:4] <- colSums(df[grep(mech, df$gene_class), 3:4])
  }
}

mech_levels <- c("AAC", "APH", "Class A, C, D", "Class B", "Erm", "Mph", "Tet efflux", "Tet enzyme", "Tet RPG", "Qnr")

# Plot figures
fig_2a <- ggplot(data = aggregated_data, aes(x = factor(mechanism, levels = mech_levels), 
                                             y = predicted_ARGs, fill = factor(mechanism, levels = mech_levels))) +
  geom_bar(position="dodge", stat="identity", color="black") +
  ggtitle("a") +
  xlab("") +
  ylab("Predicted ARGs") +
  theme_minimal() +
  scale_fill_manual(values=mechanism_palette, name = "Resistance mechanism") +
  fig_ab_theme

fig_2b <- ggplot(data = aggregated_data, aes(x = factor(mechanism, levels = mech_levels), 
                                             y = HGT_Transfers, fill = factor(mechanism, levels = mech_levels))) +
  geom_bar(position="dodge", stat="identity", color="black") +
  ggtitle("b") +
  xlab("") +
  ylab("Distantly related host pairs") +
  theme_minimal() +
  scale_fill_manual(values=mechanism_palette, name = "Resistance mechanism") +
  fig_ab_theme

