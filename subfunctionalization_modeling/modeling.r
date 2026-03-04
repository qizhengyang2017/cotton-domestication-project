# This script performs subfunctionalization analysis on cotton genomes.

library(lme4)
library(lmerTest)
library(broom.mixed)

expr_data <- read_tsv("Normalization_bias_5_timepoint.tsv")

expr_data <- expr_data %>%
  mutate(
    Group = ifelse(grepl("^NB", Sample_ID), "Wild", "Domesticated"),
    Group = factor(Group, levels = c("Domesticated", "Wild")),
    Time_num = as.numeric(gsub("DPA", "", Timepoint)),
    Time_num2 = Time_num^2
  )

# Mixed-effects linear model analysis (linear + quadratic terms + group interaction)
fit_model_interaction_quad <- function(df) {
  model <- tryCatch(
    {
      lme4::lmer(Bias ~ Time_num * Group + Time_num2 * Group + (1 | Sample_ID), data = df)
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(model)) {
    return(tibble::tibble())
  }

  coefs <- broom.mixed::tidy(model, effects = "fixed") |>
    dplyr::select("term", "estimate", "std.error", "p.value") |>
    tidyr::pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "p.value"))

  coefs
}

# Batch modeling
results <- expr_data |>
  group_by(Gene_Pair) |>
  group_modify(~ fit_model_interaction_quad(.x)) |>
  ungroup()

write_tsv(results, "bias_model_results_quad_model2.tsv")
