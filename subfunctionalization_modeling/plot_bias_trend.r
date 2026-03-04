# Performs statistical tests for linear/quadratic trends and group interactions
# Generates visualization of top gene pairs showing significant bias changes

library(patchwork)
library(ggpubr)
library(dplyr)
library(ggplot2)


results <- read_tsv("bias_model_results_quad_model2.tsv")
results_1 <- results %>%
  mutate(
    # multiple testing correction
    FDR_linear = p.adjust(p.value_Time_num, method = "BH"),
    FDR_quadratic = p.adjust(p.value_Time_num2, method = "BH"),
    FDR_interaction_linear = p.adjust(`p.value_Time_num:GroupWild`, method = "BH"),
    FDR_interaction_quadratic = p.adjust(`p.value_GroupWild:Time_num2`, method = "BH"),

    # significance flag based on FDR and effect size thresholds
    LinearTrend = FDR_linear < 0.05 & abs(estimate_Time_num) > 0.1,
    QuadraticTrend = FDR_quadratic < 0.05 & abs(estimate_Time_num2) > 0.01,
    InteractionTrend = (FDR_interaction_linear < 0.05 & abs(`estimate_Time_num:GroupWild`) > 0.1) |
      (FDR_interaction_quadratic < 0.05 & abs(`estimate_GroupWild:Time_num2`) > 0.01),
    Trend_Category = case_when(
      (LinearTrend | QuadraticTrend) & InteractionTrend ~ "Group-specific trend",
      (LinearTrend | QuadraticTrend) & !InteractionTrend ~ "Shared trend",
      !LinearTrend & !QuadraticTrend & !InteractionTrend ~ "No trend",
      TRUE ~ "Other"
    )
  )
num_sig_slope <- sum(results_1$Trend_Category %in% c("Group-specific trend", "Shared trend"), na.rm = TRUE)

num_sig_interaction <- sum(results_1$Trend_Category == "Group-specific trend", na.rm = TRUE)

cat("Analyzed", nrow(results), "gene pairs.\n")
cat(num_sig_slope, "gene pairs showed significant time-dependent trends (subfunctionalization).\n")
cat(num_sig_interaction, "gene pairs showed differences in bias trends between wild and cultivated (significant interaction effects).\n")



plot_bias_trend <- function(gene_name, expr_data, results) {
  df <- expr_data %>% filter(.data$Gene_Pair == gene_name)

  effect_label <- results %>%
    filter(.data$Gene_Pair == gene_name) %>%
    pull(.data$`estimate_Time_num:GroupWild`) %>%
    round(4)

  p_value <- results %>%
    filter(.data$Gene_Pair == gene_name) %>%
    pull(.data$`p.value_Time_num:GroupWild`) %>%
    formatC(format = "e", digits = 2)

  ggplot(df, aes(x = .data$Time_num, y = .data$Bias, color = .data$Group)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, size = 1.2) +
    labs(
      title = gene_name,
      subtitle = paste0("ΔSlope (Dom - Wild): ", effect_label, ", P = ", p_value),
      x = "Days Post Anthesis (DPA)",
      y = "Expression Bias"
    ) +
    theme_classic() +
    theme(plot.subtitle = element_text(size = 10, face = "italic"))
}

options(repr.plot.width = 12, repr.plot.height = 7)
top_genes <- results_output %>%
  arrange(FDR_interaction_linear) %>%
  slice(1:6) %>%
  pull(Gene_Pair)

library(purrr)
wrap_plots(map(top_genes, ~plot_bias_trend(.x, expr_data, results)), nrow = 2)