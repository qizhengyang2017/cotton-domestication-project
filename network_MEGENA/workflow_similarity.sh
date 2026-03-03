# =====================================================================
#  Workflow of network similarity between two populations
#  Goal: calculating Jaccard and aFETp between modules of semi-wild cotton and cultivated cotton.
#
#  Input:
#    - cultivar_*_geneID.txt : Gene IDs of each module in cultivated cotton
#    - semi_*_geneID.txt : Gene IDs of each module in semi-wild cotton
#
#  Output:
#    - Jaccard_summary.csv # Jaccard between modules
#    - summary_aFETp_table.txt # aFETp between modules
#
#  Requirements:
#    - shell script: jaccard_similarity.sh
#    - R script: fishertest.r
# =====================================================================

input_folder1="cultivar_geneID/"
input_folder2="semi_geneID/"

output_file="Jaccard_summary.csv"

echo "cultivar_module,semi_module,Jaccard_similarity" > "$output_file"

for file1 in "$input_folder1"/*; do
    filename1=$(basename "$file1")
    for file2 in "$input_folder2"/*; do
        filename2=$(basename "$file2")
        filename1_processed1="${filename1%.*}"
        filename1_processed2=$(echo "$filename1_processed1" | cut -d'_' -f2-)
        filename2_processed1="${filename2%.*}"
        filename2_processed2=$(echo "$filename2_processed1" | cut -d'_' -f2-)
		
        ## --- Step 1: Jaccard calculation
        similarity=$(sh jaccard_similarity.sh "$file1" "$file2")
        echo "$filename1_processed2,$filename2_processed2,$similarity" >> "$output_file"
		
        ## --- Step 2: FET matrix 
        file_both="Both_${filename1_processed2}_${filename2_processed2}.txt"
        cat "$file1" "$file2" | sort | uniq -d > ${file_both}
        num_both=`wc -l ${file_both} | cut -d" " -f1`
	    A=`wc -l ${file1} | cut -d" " -f1`
        Module_A=`expr $A - ${num_both}`
        B=`wc -l ${file2} | cut -d" " -f1`
        Module_B=`expr $B - ${num_both}` 
        Module_other=`expr ${all_num} - ${Module_A} - ${Module_B} - ${num_both}`
	    FET_file="./fisher_matrix/${filename1_processed2}_${filename2_processed2}.txt"
        echo "${num_both}" > "$FET_file"
	    echo "${Module_A}" >> "$FET_file"
        echo "${Module_B}" >> "$FET_file"
	    echo "${Module_other}" >> "$FET_file"
    done
done

## --- Step 3: aFETp calculation
Rscript fishertest.r
