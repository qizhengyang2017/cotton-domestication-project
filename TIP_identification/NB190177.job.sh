#!/bin/bash
set -euo pipefail

module load vg/1.53.0
tmp_dist=$(mktemp -u dist_XXXXX)

cp cotton-TE.dist "$tmp_dist"

vg giraffe -p -t 20 -Z cotton-TE.giraffe.gbz -d "$tmp_dist" -m cotton-TE.min -f "data/semi-wild/NB190177_paired_R1.fq.gz" -f "data/semi-wild/NB190177_paired_R2.fq.gz" > "NB190177.gam"
rm "$tmp_dist"
vg pack -x cotton-TE.giraffe.gbz -g "NB190177.gam" -o "NB190177.pack" -Q 5 -t 20
vg call cotton-TE.giraffe.gbz -r cotton-TE.snarls -k "NB190177.pack" -s "NB190177" -z -t 20 -a > "NB190177.vcf"
rm "NB190177.gam" "NB190177.pack"
