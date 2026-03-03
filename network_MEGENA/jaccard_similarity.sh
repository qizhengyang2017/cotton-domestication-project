# =====================================================================
#  Jaccard Calculation
#  Goal: calculating Jaccard between different modules.
# =====================================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file1 file2"
    exit 1
fi

set1=($(cat "$1"))
set2=($(cat "$2"))

tmpfile1=$(mktemp)
tmpfile2=$(mktemp)

printf "%s\n" "${set1[@]}" | sort > "$tmpfile1"
printf "%s\n" "${set2[@]}" | sort > "$tmpfile2"

intersection=($(comm -12 "$tmpfile1" "$tmpfile2"))
intersection_size=${#intersection[@]}

union=($(cat "$tmpfile1" "$tmpfile2" | sort -u))
union_size=${#union[@]}

rm "$tmpfile1" "$tmpfile2"

if [ "$union_size" -eq 0 ]; then
    similarity=0
else
    similarity=$(awk "BEGIN { print $intersection_size / $union_size }")
fi

echo "$similarity"