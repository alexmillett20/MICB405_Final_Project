#!/bin/bash

# Input CSV file
input="SraRunTable_PRJNA1240347.csv"
# Output file
output="sample_list.csv"

# Skip header, extract Run (1st col), culture_replicate (13th col), and treatment (last col)
awk -F',' 'NR>1 {print $1 ", " $NF ", " $13}' "$input" > "$output"

echo "✅ Extracted sample info saved to $output"

