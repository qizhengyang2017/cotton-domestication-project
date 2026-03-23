awk 'BEGIN{OFS="\t"} /^#/ {print; next} { $4=toupper($4); $5=toupper($5); print }' pangenome.vcf > pangenome.upper.vcf
