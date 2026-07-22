library(coloc)
library(tidyverse)
library(arrow)

# setwd('~/semi_wild/new_analysis1124/bias_coloc')

args <- commandArgs(trailingOnly = TRUE)
time <- args[1]
chr <- args[2]

eqtl <- read_parquet(paste0("~/semi_wild/new_analysis1124/bias-eQTL/parquet_file/",time,".cis_qtl_pairs.",chr,".parquet"))
# gwas <- read_tsv('~/semi_wild/GWAS/GWAS_15M/GWAS_ourlab_phenotype_filter/FL_out')
gwas <- read_tsv('/path/to/gwas/trait_blup_376_filter_ranLineRep/BLUP_UHML_out')
#gwas <- read_tsv('~/semi_wild/GWAS/GWAS_15M/GWAS_zms/FL_out')
#eqtl <- read_parquet("~/semi_wild/new_analysis1124/bias-eQTL/4DPA.cis_qtl_pairs.1.parquet")
#gwas <- read_tsv('~/semi_wild/GWAS/GWAS_15M/GWAS_ourlab_phenotype_filter/FL_out')
eqtl_selected <- eqtl |>  
  select(phenotype_id, variant_id, pval_nominal, slope, slope_se) |> 
  rename(beta = slope, varbeta = slope_se) |> 
  mutate(varbeta = varbeta^2)  # 将 slope_se 转换为效应值的方差
gwas_selected <- gwas |>  
  select(SNP, Pvalue, SNPWeight, SNPWeightSE) |> 
  rename(variant_id = SNP, 
         pval_nominal = Pvalue, 
         beta = SNPWeight, 
         varbeta = SNPWeightSE) |> 
  mutate(varbeta = varbeta^2)  # 将 SNPWeightSE 转换为效应值的方差
input <- merge(eqtl_selected, gwas_selected, by="variant_id", all=FALSE, suffixes=c("_eqtl","_gwas"))

out.GENE <- c()
out.COLOC.PP0 <- c()
out.COLOC.PP1 <- c()
out.COLOC.PP2 <- c()
out.COLOC.PP3 <- c()
out.COLOC.PP4 <- c()
out.CAUSAL.SNP <- c()
out.SNP.PP.H4 <- c()

# 按 phenotype_id 分组
gene_groups <- split(input, input$phenotype_id)

# 循环遍历每个基因
for (gene in names(gene_groups)) {
  input_gene <- gene_groups[[gene]] |>
      drop_na(beta_eqtl, varbeta_eqtl, beta_gwas, varbeta_gwas)
  
  clc <- coloc.abf(
    dataset1 = list(
      snp = input_gene$variant_id,
      beta = input_gene$beta_gwas, 
      varbeta = input_gene$varbeta_gwas, 
      type = "quant",  
      N = 400,
      sdY=1  
    ), 
    dataset2 = list(
      snp = input_gene$variant_id,
      beta = input_gene$beta_eqtl, 
      varbeta = input_gene$varbeta_eqtl, 
      type = "quant",  
      N = 400,
      sdY=1
    ) 
  )
  PP4 <- round(clc$summary[6],3)
  m <- which.max(clc$results$SNP.PP.H4)
  causal_snp <- clc$results$snp[m]
  snp_H4 <- clc$results$SNP.PP.H4[m]

  out.GENE <- c(out.GENE, gene)
  out.COLOC.PP0 <- c(out.COLOC.PP0, round(clc$summary[2],3))
  out.COLOC.PP1 <- c(out.COLOC.PP1, round(clc$summary[3],3))
  out.COLOC.PP2 <- c(out.COLOC.PP2, round(clc$summary[4],3))
  out.COLOC.PP3 <- c(out.COLOC.PP3, round(clc$summary[5],3))
  out.COLOC.PP4 <- c(out.COLOC.PP4, round(clc$summary[6],3))
  out.CAUSAL.SNP <- c(out.CAUSAL.SNP, causal_snp)
  out.SNP.PP.H4 <- c(out.SNP.PP.H4, snp_H4)

}

out.tbl = data.frame(out.GENE, out.COLOC.PP0, out.COLOC.PP1, out.COLOC.PP2,
                     out.COLOC.PP3, out.COLOC.PP4, out.CAUSAL.SNP, out.SNP.PP.H4)

write.table(out.tbl,quote=F , row.names=F , sep='\t' , file=paste0(time,'_',chr,'_coloc_result.txt'))

out.tbl.confident <- out.tbl |> filter(out.COLOC.PP4 > 0.8)
write.table(out.tbl.confident ,quote=F , row.names=F , sep='\t' , file=paste0(time,'_',chr,'_coloc_result_confident.txt'))
