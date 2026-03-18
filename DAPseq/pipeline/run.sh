module load snakemake/5.4.0
mkdir -p logs/cluster
snakemake  --cluster-config cluster.json --cluster  \
'bsub -n {threads} -q {cluster.queue} -R {cluster.resources} -o {cluster.output} -e {cluster.error}'  -j 100  -s Chipseq.smk_new1.py
