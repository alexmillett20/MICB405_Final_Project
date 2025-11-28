#!/usr/bin/env python3
# Author: Brendan Ng
# MICB405 Final Project GeneID to GO Mapping Script
# Date created: November 18, 2025
# Last updated: November 18, 2025

# This script maps Entrez GeneIDs from a GTF file to GO terms using the gene2go database.
# The GENE2GO_FILE and GTF_FILE were downloaded from NCBI and stored on the MICB 405 server.

import re
from collections import defaultdict

# === INPUT FILES ===
GTF_FILE = "../../data/GCF_000001635.27_GRCm39_genomic.gtf"
GENE2GO_FILE = "../../data/gene2go"

# === OUTPUT FILE ===
OUTPUT_FILE = "../gene_GO_mapping.tsv"


print("=== Step 1: Extracting Entrez GeneIDs from GTF ===")

gene_ids = set()
pattern = re.compile(r'GeneID:(\d+)')

with open(GTF_FILE) as f:
    for line in f:
        match = pattern.search(line)
        if match:
            gene_ids.add(match.group(1))

print(f"Found {len(gene_ids)} Entrez GeneIDs in GTF")


print("=== Step 2: Reading gene2go and mapping to GO terms ===")

gene_to_go = defaultdict(set)

with open(GENE2GO_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        
        tax_id, gid, go_id = parts[0], parts[1], parts[2]

        # Mouse = tax_id 10090
        if tax_id == "10090" and gid in gene_ids:
            gene_to_go[gid].add(go_id)

print(f"Mapped GO terms to {len(gene_to_go)} genes")


print("=== Step 3: Writing topGO mapping file ===")

with open(OUTPUT_FILE, "w") as out:
    for gene, go_terms in gene_to_go.items():
        go_list = ",".join(sorted(go_terms))
        out.write(f"{gene}\t{go_list}\n")

print(f"Done! Saved topGO file: {OUTPUT_FILE}")

