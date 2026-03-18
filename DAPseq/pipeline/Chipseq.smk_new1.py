# _*_ coding: UTF-8 _*_
########################################################################
# Qi Zhengyang oringnal 2023-12-21
# DAP-seq pipeline
######################################################################## 
# run on legend
# MAKE SURE WORK ON CONDA ENV PY3.6


SAMPLE = ["MY2","TC2","TM_input"]
BT2_IDX = "genome/Ghirsutum_genome"
species="tm1"

rule all:
	input:
		bamCorver=expand( 
			'3.Coverage/{sample}.bt2.{species}.q15.rmdup.sort.bw',
			species=['tm1'],sample=SAMPLE
		),
		peakFile=expand(
			"4.macs2_result/{sample}_{species}_peaks.narrowPeak",
			 species=['tm1'],sample=["MY2","TC2"]
		)
rule bowtie2_mapping:
	input:
		"clean_data/{sample}_1.clean.fq.gz",
		"clean_data/{sample}_2.clean.fq.gz"
	output:
		"2.bowtie2_results/{sample}.bt2.{species}.sam",
		"2.bowtie2_results/{sample}.bt2.{species}.bam",
		"2.bowtie2_results/{sample}.bt2.{species}.sort.rmdup.bam",
		"2.bowtie2_results/{sample}.bt2.{species}.matrix",
		"2.bowtie2_results/{sample}.bt2.{species}.q15.rmdup.sort.bam"
	log:
		"2.bowtie2_results/{sample}.bt2.{species}.log",
		"2.bowtie2_results/{sample}.bt2.{species}.rmdup.log"
	threads: 20
	shell:
		"""
		module load Bowtie2/2.5.2
		module load picard/2.23.9
		module load sambamba/0.8.2
		module load SAMtools/1.9

		bowtie2 -p {threads} \
		    -x {BT2_IDX} \
		    -1 {input[0]} -2 {input[1]} \
		    -S {output[0]} > {log[0]} 2>&1 && \
		samtools sort -@ {threads} -o {output[1]} {output[0]} && \
		java -Xms12g -Xmx12g -XX:ParallelGCThreads={threads} \
		    -jar ${{EBROOTPICARD}}/picard.jar MarkDuplicates \
		    I={output[1]} O={output[2]} M={output[3]} \
		    ASO=coordinate REMOVE_DUPLICATES=true 2>{log[1]}
		samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 15 {output[2]} | \
		    samtools sort -@ {threads} -o {output[4]} -
		samtools index {output[4]}
		"""


rule bamcoverage:
	input:
		"2.bowtie2_results/{sample}.bt2.{species}.q15.rmdup.sort.bam"
	output:
		"3.Coverage/{sample}.bt2.{species}.q15.rmdup.sort.bw"
	threads: 4
	shell:
		"""
		module load deepTools/3.5.0
		bamCoverage -b {input} -o {output} -of bigwig --binSize 1 -p 10 --normalizeUsing RPKM
		"""

rule call_peak:
	input:
		treat="2.bowtie2_results/{sample}.bt2.{species}.q15.rmdup.sort.bam",
		control="2.bowtie2_results/TM_input.bt2.{species}.q15.rmdup.sort.bam"
	output:
		"4.macs2_result/{sample}_{species}_peaks.narrowPeak"
	params:
		head_outfile="{sample}_{species}",
		outdir="4.macs2_result"
	log:
		"4.macs2_result/{sample}_{species}.macs2.log"
	threads: 4
	shell:
		"""
		module load MACS2/2.2.7.1
		macs2 callpeak -t {input.treat} -c {input.control} -f BAM -g 23e8 \
		--outdir {params.outdir} -n {params.head_outfile} \
		-B -q 0.05 > {log} 2>&1
		"""

