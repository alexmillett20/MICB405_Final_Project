#!/bin/bash
# Author: Brendan Ng
# MICB405 Final Project Extract Sample Info From SraRunTable Script
# Date created: October 31, 2025
# Last updated: October 31, 2025

# This script extracts sample information from the SraRunTable CSV file.
# It retrieves the Run accession, treatment, and culture replicate number.

# ChatGPT was used to help generate this script.

# Input CSV file
input="SraRunTable_PRJNA1240347.csv"
# Output file
output="sample_list.csv"

# Skip header, extract Run (1st col), culture_replicate (13th col), and treatment (last col)
awk -F',' 'NR>1 {print $1 ", " $NF ", " $13}' "$input" > "$output"

echo "✅ Extracted sample info saved to $output"

