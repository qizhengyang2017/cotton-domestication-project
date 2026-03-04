# Generate boxplots comparing median bias between wild and domesticated cotton groups across developmental timepoints

library(tidyverse)
library(ggpubr)

expr_data <- read_tsv("Normalized_bias_expression.tsv") %>%
  mutate(
    Group = ifelse(grepl("^NB", Sample_ID), "Wild", "Domesticated"),
    Group = factor(Group, levels = c("Wild", "Domesticated")),
    Time_num = as.numeric(gsub("DPA", "", Timepoint)),
    Timepoint = factor(Timepoint, levels = c("4DPA", "8DPA", "12DPA", "16DPA", "20DPA"))
  )

# at least 10% of samples have At or Dt > 1:
high_expr_genes <- expr_data %>%
  group_by(Gene_Pair) %>%
  summarize(high_expr_frac = mean(FPKM_At > 1 | FPKM_Dt > 1)) %>%
  filter(high_expr_frac >= 0.1) %>%
  pull(Gene_Pair)
expr_data <- expr_data %>% filter(Gene_Pair %in% high_expr_genes)

write_tsv(expr_data, 'filter_bias_10percent_fpkm1.txt')


# Step 1: calculate median bias for each gene pair by timepoint and group
gene_bias_median_tp <- expr_data %>%
  group_by(Gene_Pair, Group, Timepoint) %>%
  summarise(median_bias = median(abs(Bias), na.rm = TRUE), .groups = "drop")


gene_bias_median_tp$Group <- recode(gene_bias_median_tp$Group, "Domesticated" = "Cultivated", "Wild" = "Semi-wild")


# Step 2: plot boxplots with significance annotations, faceted by timepoint
options(repr.plot.width = 7, repr.plot.height = 4)
ggplot(gene_bias_median_tp, aes(x = Group, y = median_bias, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
  # geom_jitter(width = 0.2, size = 0.7, alpha = 0.4) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = "top",
    size = 4
  ) +
  facet_wrap(~Timepoint, nrow = 1) +
  scale_fill_manual(values = c("Semi-wild" = "#084594", "Cultivated" = "#b10026")) +
  labs(
    x = "Group",
    y = "Median bias per gene pair"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1)
  )
ggsave("bias_new1.pdf", width = 7, height = 4)
