# Purpose: Compare bias values between domesticated and wild groups across different timepoints,
#          perform t-tests, apply FDR correction, and filter for significant differences

library(tidyverse)

df <- read_tsv("filter_bias_10percent_fpkm1.txt")

results <- list()
gene_pairs <- unique(df$Gene_Pair)
timepoints <- unique(df$Timepoint)

# iterate over all gene-pair/timepoint combinations
for (gene in gene_pairs) {
  for (tp in timepoints) {
    subset_df <- df %>%
      filter(Gene_Pair == gene, Timepoint == tp)

    groups_present <- unique(subset_df$Group)

    # ensure both groups are present
    if (all(c("Domesticated", "Wild") %in% groups_present)) {
      dom_bias <- subset_df %>% filter(Group == "Domesticated") %>% pull(Bias)
      wild_bias <- subset_df %>% filter(Group == "Wild") %>% pull(Bias)

      # ensure each group has at least two samples
      if (length(dom_bias) >= 2 && length(wild_bias) >= 2) {
        ttest <- t.test(dom_bias, wild_bias)

        results[[length(results) + 1]] <- tibble(
          Gene_Pair = gene,
          Timepoint = tp,
          p_value = ttest$p.value,
          mean_bias_domesticated = mean(dom_bias),
          mean_bias_wild = mean(wild_bias),
          bias_diff = mean(wild_bias) - mean(dom_bias)  # 保留方向信息
        )
      }
    }
  }
}

# combine all results
result_df <- bind_rows(results)

# FDR correction (Benjamini-Hochberg)
result_df <- result_df %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

write_tsv(result_df, "bias_change_gene_pairs_FDR_all.tsv")

# filter significant results
significant_results <- result_df %>%
  filter(p_adj < 0.05, abs(bias_diff) > 0.2)  # 注意加 abs 保留方向但设阈值时不受影响

output_file <- "significant_bias_change_gene_pairs_FDR_new_1.tsv"
write_tsv(significant_results, output_file)
