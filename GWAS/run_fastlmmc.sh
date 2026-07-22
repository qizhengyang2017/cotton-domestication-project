module load Fastlmm/v0.2.32

read -a PHENO_NAME < header.txt



for((i=1;i<=5;i++))
do

bsub -J ${PHENO_NAME[i-1]} -o ${PHENO_NAME[i-1]}_%J.out -e ${PHENO_NAME[i-1]}_%J.err -q normal -n 5 "
 
fastlmmc -bfile semi_161 -bfilesim semi_161 -pheno phenotype_zhongmiansuo.txt -covar plink.pca.eigenvec -mpheno ${i} -out ${PHENO_NAME[i-1]}_out

cat ${PHENO_NAME[i-1]}_out|awk 'NR>1{print \$1,\$2,\$4,\$5}' OFS='\t' > ${PHENO_NAME[i-1]}_out.manhattan.txt

sed -i '1iSNP\tCHR\tBP\tP' ${PHENO_NAME[i-1]}_out.manhattan.txt

Rscript qqman.r ${PHENO_NAME[i-1]}_out.manhattan.txt ${PHENO_NAME[i-1]}_manhattan.png ${PHENO_NAME[i-1]}_qq.png
"
done 
