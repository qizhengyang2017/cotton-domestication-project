bsub -q normal -K -J joint -n 36 -o %J.joint.out -e %J.joint.err -R span[hosts=1] "
module load sentieon/202112
cat gvcf_list.txt |sentieon driver -t 36 -r /path/to/vcf_call_101/Ghirsutum_genome.fasta \
--algo GVCFtyper output-joint.vcf.gz -"

bsub -q normal -K -J filter -n 1 -o %J.filter.out -e %J.filter.err -R span[hosts=1] "
sh filter.sh"
