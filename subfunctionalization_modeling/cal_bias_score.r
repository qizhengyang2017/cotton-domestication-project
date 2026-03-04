# Purpose: Calculate Homologous Gene Expression Bias Ratio
# The bias ratio is computed as: (FPKM_At - FPKM_Dt)/(FPKM_At + FPKM_Dt)

library(tidyverse)

# set file paths and timepoint labels
timepoints <- c("4DPA", "8DPA", "12DPA", "16DPA", "20DPA")
files <- paste0("effective_exp_gene/cultivar_semi_merge_", timepoints, "_filter.txt")
names(files) <- timepoints

# read homologous gene pairs
homologs <- read_tsv("/data/cotton/jqyou/hg_task/allm_hom_20189_name.txt", col_names = c("Gene_At", "Gene_Dt"))

# read data for all timepoints into a list
expression_list <- map(files, ~ read_tsv(.x))

# add timepoint label to each expression dataset
long_data_list <- map2(expression_list, names(expression_list), function(df, tp) {
  df_long <- df %>%
    pivot_longer(-GeneID, names_to = "Sample_ID", values_to = "FPKM") %>%
    mutate(Timepoint = tp)
  df_long
})

# combine all timepoint long-format data
all_expr <- bind_rows(long_data_list)

# extract At and Dt expression data separately
expr_At <- all_expr %>%
  filter(GeneID %in% homologs$Gene_At) %>%
  rename(Gene_At = GeneID, FPKM_At = FPKM)

expr_Dt <- all_expr %>%
  filter(GeneID %in% homologs$Gene_Dt) %>%
  rename(Gene_Dt = GeneID, FPKM_Dt = FPKM)

# merge expression data
merged <- homologs %>%
  inner_join(expr_At, by = "Gene_At") %>%
  inner_join(expr_Dt, by = c("Gene_Dt", "Sample_ID", "Timepoint")) %>%
  mutate(
    Gene_Pair = paste(Gene_At, Gene_Dt, sep = "-"),
    Bias = case_when(
      FPKM_At == 0 & FPKM_Dt == 0 ~ 0,
      TRUE ~ (FPKM_At - FPKM_Dt) / (FPKM_At + FPKM_Dt)
    )
  ) %>%
  select(Gene_Pair, Timepoint, Sample_ID, FPKM_At, FPKM_Dt, Bias)


write_tsv(merged, "Normalized_bias_expression.tsv")
