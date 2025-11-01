#!/bin/bash

mapping_file="srr_sample_mapping.txt"
reads_dir="reads_per_gene"

# Loop through each line in the mapping file
while IFS=',' read -r srr treatment rep; do
    # Trim whitespace
    srr=$(echo "$srr" | xargs)
    treatment=$(echo "$treatment" | xargs)
    rep=$(echo "$rep" | xargs)

    old_file="${reads_dir}/${srr}_ReadsPerGene.out.tab"
    new_file="${reads_dir}/${treatment}_rep_${rep}_ReadsPerGene.out.tab"

    if [[ -f "$old_file" ]]; then
        echo "Renaming $old_file → $new_file"
        mv "$old_file" "$new_file"
    else
        echo "⚠️  Skipping: $old_file not found"
    fi
done < "$mapping_file"
