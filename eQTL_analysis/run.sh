prefix=$1
expression_bed=phenotypes/${prefix}_Sample_fpkm_filter_qqnorm_sorted.bed.gz
covariates_file=covariant/${prefix}_covariants.txt

plink_prefix_path=TIP-new-pop537/pop537_GT_numberChr
echo ${covariates_file}
echo ${expression_bed}


python -c "import torch; print(torch.__version__); print('CUDA available: {} ({})'.format(torch.cuda.is_available(), torch.cuda.get_device_name(torch.cuda.current_device())))"

# get permutation result
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis


python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --cis_output ${prefix}.cis_qtl.txt.gz \
    --mode cis_independent

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode trans

