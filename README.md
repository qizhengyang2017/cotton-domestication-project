## Overview

Starting from population whole-genome sequencing and developmental-stage fiber RNA-seq, the pipeline quantifies homoeologous expression bias (HEB) and its developmental dynamics, maps cis-eQTL and bias-eQTL to find regulatory variants underlying HEB, integrates ATAC-seq / DAP-seq to assess chromatin context, colocalizes bias-eQTL with fiber-quality GWAS signals, and reconstructs co-expression networks to link regulatory divergence to fiber traits.



## Modules

### `SNP_calling/` — Variant calling
Population-scale variant calling from whole-genome sequencing.
- `sra.sh` — download raw reads from SRA and convert to FASTQ.
- `fastp.sh` — read quality trimming with fastp.
- `snp_calling_sentieon.sh` — BWA-MEM alignment and per-sample variant calling with Sentieon.
- `joint_calling_filter.sh` — joint genotyping across samples.
- `filter.sh` — SNP selection and hard filtering (GATK).

### `GWAS/` — Association mapping
- `run_fastlmmc.sh` — mixed-model GWAS with FaST-LMM using PCA covariates, followed by Manhattan/QQ plotting (`qqman`).

### `eQTL_analysis/` — eQTL mapping
- `run.sh` — cis-eQTL, cis-independent, and trans mapping with [tensorQTL](https://github.com/broadinstitute/tensorqtl), using normalized expression BED files and genotype (PLINK) input. The same framework is applied to HEB phenotypes to map bias-eQTL.

### `bias_coloc/` — Colocalization
- `coloc_batch.r` — Bayesian colocalization (`coloc`) between bias-eQTL and fiber-quality GWAS loci to identify shared causal variants.

### `TIP_identification/` — TE insertion polymorphisms
Genotyping transposable-element insertion polymorphisms (TIPs) against a pangenome graph using [vg giraffe](https://github.com/vgteam/vg). 
- `toupper.sh` — normalize the pangenome VCF.
- `vg_autoindex_new.sh` — build the giraffe index from the reference and pangenome VCF.
- `vg_snarls.sh` — compute snarls for variant calling.
- `call_cultivated_new.sh` / `call_wild_new.sh` — generate per-sample genotyping jobs for cultivated and semi-wild accessions.
- `NB190177.job.sh` — example per-sample job (map → pack → call).



### `DAPseq/` — DAP-seq
- `pipeline/Chipseq.smk_new1.py` — Snakemake pipeline for DAP-seq (alignment, peak calling with MACS2); `run.sh`, `cluster.json` for cluster execution.
- `DAPseq_analysis.r` — peak annotation and analysis

### `network_MEGENA/` — Co-expression networks
Co-expression network reconstruction with [MEGENA](https://cran.r-project.org/package=MEGENA) and cross-population module comparison.
- `MEGENA_construction.r` — build the network from an expression matrix.
- `eigengenes_calculate.r` — module eigengenes and module–phenotype relationships.
- `Hierarchical_division.r` — hierarchical module structure.
- `jaccard_similarity.sh`, `workflow_similarity.sh`, `fishertest.r` — module similarity between semi-wild and cultivated networks (Jaccard and aggregated Fisher's exact test, aFETp).

### `subfunctionalization_modeling/` — HEB analysis
Quantification and modeling of homoeologous expression bias across five fiber developmental stages (4/8/12/16/20 DPA).
- `cal_bias_score.r` — bias ratio (FPKM_At − FPKM_Dt)/(FPKM_At + FPKM_Dt) per homoeolog pair.
- `modeling.r` — mixed-model analysis of bias (`lme4`/`lmerTest`).
- `significant_bias_change_pairs.r` — t-tests with FDR correction for bias differences between wild and domesticated groups.
- `plot_bias_trend.r`, `bias_boxplot.r`, `bias_heatmap.r` — visualization of bias trends and distributions.

## Requirements

Analyses were run on a Linux HPC cluster (LSF/`bsub`). Key software:

- **Read processing & alignment:** fastp, BWA, Sentieon, samtools, GATK 4.1.9
- **Data download:** SRA Toolkit 3.0
- **eQTL:** Python 3 + PyTorch, tensorQTL, PLINK
- **GWAS:** FaST-LMM (fastlmmc)
- **Pangenome / TIP:** vg 1.53.0
- **DAP-seq:** Snakemake, Bowtie2, MACS2
- **R (≥4.0):** tidyverse, coloc, MEGENA, lme4, lmerTest, broom.mixed, ComplexHeatmap, ggpubr, patchwork, arrow

