module load fastp/0.23.0

mkdir report clean_data

basedir=/path/to/du
cat fastp_sample.txt|while read i
do
report=$basedir/report/$i.html
json=$basedir/report/$i.json

fq1_raw=$basedir/fastq_data/${i}_1.fastq.gz
fq2_raw=$basedir/fastq_data/${i}_2.fastq.gz

fq1=$basedir/clean_data/${i}_clean.1.fq.gz
fq2=$basedir/clean_data/${i}_clean.2.fq.gz

# thread
nt=5
# queue
bsub -J ${i} -o %J.${i}.out -e %J.${i}.err -n $nt -q normal "
fastp -w $nt -i $fq1_raw -I $fq2_raw -o $fq1 -O $fq2 -h $report -j $json && \
rm $fq1_raw $fq2_raw
"
done
