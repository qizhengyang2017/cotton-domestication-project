# =================================
# Workflow for processing ATAC-seq data
#  Goal: processing ATAC-seq data from multiple samples and identifying
#        accessible chromatin peaks located in gene promoter regions.
#
#  Input:
#    - Raw sequencing data (*.fq.gz)
#
#  Output:
#    - ATAC-seq peaks for each sample (*peaks.narrowPeak)
#    - Peaks overlapping gene promoter regions (inter_*.txt)
#
# =================================

sample_list=("L248_1" "L347_1" "L289_1" "L209_1" "L128_1" "L248_2" "L347_2" "L289_2" "L209_2" "L128_2")
raw_dir="ATAC_seq/rawdata"
result_dir="ATAC_seq"
genome_dir="ATAC_seq/genome"

bowtie2-build -f ${genome_dir}/Ghirsutum_genome.fasta --threads 5 ${genome_dir}/Ghirsutum_genome

## --- Step1: call peaks
for sample in "${sample_list[@]}"; do
    raw_data1=${raw_dir}/${sample}_R1
    raw_data2=${raw_dir}/${sample}_R2
    bam_file=${result_dir}/bowtie2/${sample}_merge
    clean_data1=${result_dir}/data_trim/${sample}_R1
    clean_data2=${result_dir}/data_trim/${sample}_R2
	
    trim_galore -q 20 --phred33 --stringency 3 --length 20 --paired --gzip --trim-n ${raw_data1}.fq.gz ${raw_data2}.fq.gz -o ${result_dir}/data_trim 
    bowtie2 -p 2 -x ${genome_dir}/Ghir -1 ${clean_data1}_val_1.fq.gz -2 ${clean_data2}_val_2.fq.gz | samtools sort -@ 2 -o ${bam_file}.sorted.bam > ${bam_file}.log 2> ${bam_file}.err 
    samtools rmdup ${bam_file}.sorted.bam ${bam_file}.sorted.rmdup.bam 

    macs2 callpeak -t ${bam_file}.sorted.rmdup.bam -f BAM -g 2348137562 -n ${sample} -B -q 0.01 --outdir ${result_dir}/peaks/${sample}
done

## --- Step2: merge replication
sample_list=("L248" "L347" "L289" "L209" "L128")
peak_dir="ATAC_seq/peaks"
idr_dir="ATAC_seq/idr_result"
promoter_file="ATAC_seq/genome/TM1_promoter.bed"

for sample in "${sample_list[@]}"; do
    peak1="${peak_dir}/${sample}_1/${sample}_1_peaks.narrowPeak"
    peak2="${peak_dir}/${sample}_2/${sample}_2_peaks.narrowPeak"
    out_file="${idr_dir}/${sample}_IDR.txt"
    clean_file="${idr_dir}/clean_${sample}_IDR.txt"
    idr --samples ${peak1} ${peak2} --input-file-type narrowPeak --rank signal.value --output-file ${out_file} --plot
    awk '$5 > 540' ${out_file} > ${clean_file}

    inter_file="${idr_dir}/inter_${sample}.txt"
    bedtools intersect -a ${promoter_file} -b ${clean_file} -c > ${inter_file}
done

