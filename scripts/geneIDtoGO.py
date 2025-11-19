import pandas as pd
import re

# 1. Parse GTF and extract Entrez GeneIDs
gtf_file = "../../data/GCF_000001635.27_GRCm39_genomic.gtf"

gene_ids = set()
pattern = re.compile(r'GeneID:(\d+)')

with open(gtf_file) as f:
    for line in f:
        if "GeneID:" in line:
            match = pattern.search(line)
            if match:
                gene_ids.add(match.group(1))

print(f"Found {len(gene_ids)} unique Entrez IDs in GTF.")

# 2. Load gene2go
gene2go = pd.read_csv(
    "../../data/gene2go",
    sep="\t",
    comment="#",
    names=["tax_id","GeneID","GO_ID","Evidence","Qualifier","GO_term","PubMed","Category"]
)

# 3. Filter for mouse (TaxID 10090)
gene2go_mouse = gene2go[gene2go["tax_id"] == 10090]

# 4. Keep only genes present in your GTF
gene_go_map = gene2go_mouse[gene2go_mouse["GeneID"].astype(str).isin(gene_ids)]

# 5. Save output
gene_go_map[["GeneID","GO_ID"]].drop_duplicates().to_csv(
    "../gene_GO_mapping.csv", index=False
)

print("Saved gene_GO_mapping.csv")

