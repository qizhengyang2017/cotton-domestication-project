while IFS= read -r i; do
    fq1="data/semi-wild/${i}_paired_R1.fq.gz"
    fq2="data/semi-wild/${i}_paired_R2.fq.gz"
    echo "Generating job script for sample: $i"
    nt=20
    mkdir -p sample_scripts
    script_name="sample_scripts/${i}.job.sh"

    cat <<EOF > "$script_name"
#!/bin/bash
set -euo pipefail

module load vg/1.53.0
tmp_dist=\$(mktemp -u dist_XXXXX)

cp cotton-TE.dist "\$tmp_dist"

vg giraffe -p -t $nt -Z cotton-TE.giraffe.gbz -d "\$tmp_dist" -m cotton-TE.min -f "$fq1" -f "$fq2" > "${i}.gam"
rm "\$tmp_dist"
vg pack -x cotton-TE.giraffe.gbz -g "${i}.gam" -o "${i}.pack" -Q 5 -t $nt
vg call cotton-TE.giraffe.gbz -r cotton-TE.snarls -k "${i}.pack" -s "$i" -z -t $nt -a > "${i}.vcf"
rm "${i}.gam" "${i}.pack"
EOF

    chmod +x "$script_name"

    echo "Saved job script: $script_name"
    cat "$script_name"
    echo "----------------------------------"

    bsub -q normal -J "$i" -n "$nt" -o "${i}_%J.out" -e "${i}_%J.err" -R "span[hosts=1]" bash "$script_name"

done < wild_redo.txt
