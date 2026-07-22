
module load  sratoolkit/3.0.0

# download

mkdir sra_data
cat remain_466_srr.txt|while read -a NAME
do
bsub -q normal -J ${NAME[1]} -o ${NAME[1]}.%J.out -e ${NAME[1]}.%J.err "
prefetch --max-size 100G ${NAME[0]} -o sra_data/${NAME[1]}.sra
"
done



# sra to fastq

nt=5
mkdir fastq_data
for i in `ls sra_data`
do
NAME=${i/.sra/}
bsub -q normal -J $NAME -o $NAME.%J.out -e $NAME.%J.err -n $nt "
fasterq-dump --threads $nt --split-3 sra_data/$NAME.sra -O fastq_data && \
pigz -p $nt fastq_data/${NAME}_1.fastq && \
pigz -p $nt fastq_data/${NAME}_2.fastq && \
rm sra_data/$NAME.sra
"
done
