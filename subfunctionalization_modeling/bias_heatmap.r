library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# read TF annotation file
tf_annot <- read.delim("~/semi_wild/new_analysis1124/PROJECT_DATA/raw/Ghir_AT_annot_70199_TF.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tf_annot <- tf_annot[, c("Cotton", "TF")]
colnames(tf_annot) <- c("Gene", "TF")

# read bias result file
bias_df <- read.delim("bias_model_results_quad_model2_sig.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# split Gene_Pair into Gene_A and Gene_D
library(tidyr)
bias_df <- separate(bias_df, col = Gene_Pair, into = c("Gene_A", "Gene_D"), sep = "-", remove = FALSE)

bias_df <- merge(bias_df, tf_annot, by.x = "Gene_A", by.y = "Gene", all.x = TRUE)
colnames(bias_df)[ncol(bias_df)] <- "TF_A"

bias_df <- merge(bias_df, tf_annot, by.x = "Gene_D", by.y = "Gene", all.x = TRUE)
colnames(bias_df)[ncol(bias_df)] <- "TF_D"

# write.table(bias_df, file = "bias_model_results_with_TF.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# get significantly subfunctionalized genes
sig_genes <- bias_df %>%
  filter(Trend_Category == "Group-specific trend") %>%
  pull(Gene_Pair)

bias_long <- expr_data %>%
  filter(Gene_Pair %in% sig_genes) %>%
  mutate(Sample_Label = paste0(Timepoint, "_", Sample_ID)) %>%
  select(Gene_Pair, Sample_Label, Bias, Group, Timepoint)

bias_matrix <- bias_long %>%
  select(Gene_Pair, Sample_Label, Bias) %>%
  pivot_wider(names_from = Sample_Label, values_from = Bias) %>%
  column_to_rownames("Gene_Pair") %>%
  as.matrix()

sample_anno <- bias_long %>%
  distinct(Sample_Label, Group, Timepoint) %>%
  arrange(factor(Group, levels = c("Wild", "Domesticated")), factor(Timepoint, levels = c("4DPA", "8DPA", "12DPA", "16DPA", "20DPA")))

bias_matrix <- bias_matrix[, sample_anno$Sample_Label]

sample_anno <- sample_anno %>%
  mutate(
    Group = factor(Group, levels = c("Wild", "Domesticated")),
    Timepoint = factor(Timepoint, levels = c("4DPA", "8DPA", "12DPA", "16DPA", "20DPA"))
  )

time_levels <- sort(unique(as.character(sample_anno$Timepoint)))

# set timepoint colors
time_colors <- setNames(
  RColorBrewer::brewer.pal(length(time_levels), "Set3"),
  time_levels
)

time_colors

# top annotation: Group + Timepoint
top_anno <- HeatmapAnnotation(
  Group = sample_anno$Group,
  Timepoint = sample_anno$Timepoint,
  col = list(
    Group = c(Domesticated = "#d95f02", Wild = "#1b9e77"),
    Timepoint = time_colors
  ),
  annotation_name_side = "left"
)


### annotate genes
highlight_genes <- tf_gene_pair |>
  filter(TF_A == "bHLH" | TF_A == "MYB") |>
  pull(Gene_Pair)
highlight_genes_symbol <- tf_gene_pair |>
  filter(TF_A == "bHLH" | TF_A == "MYB") |>
  pull(TF_A)
# create combined labels
combined_labels <- paste0(highlight_genes, " (", highlight_genes_symbol, ")")

# construct named vector: names are Gene_Pair, values are combined labels
gene_label_map <- setNames(combined_labels, highlight_genes)

# construct row_labels: only show those matching bias_matrix row names
row_labels <- ifelse(
  rownames(bias_matrix) %in% names(gene_label_map),
  gene_label_map[rownames(bias_matrix)],
  ""
)
pdf("bias_heatmap5.pdf", width = 6, height = 7)

Heatmap(
  bias_matrix,
  name = "Bias",
  top_annotation = top_anno,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_labels = row_labels,
  row_names_gp = gpar(
    fontsize = 6,
    col = ifelse(row_labels != "", "black", NA)
  ),
  column_split = sample_anno$Group
)

dev.off()
