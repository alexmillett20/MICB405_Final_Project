# MICB405_Final_Project

## Overview

This repository contains the documentation for the MICB 405 final project. This repository contains
the scripts used in R to generate the differential expression results, and the README.md will contain the steps
for data processing.

### **Dataset**:

#### Step 1: download fastq files

To download the fastq files for the dataset, the following instructions were used:

```bash
esearch -db sra -query PRJNA1240347 \
| efetch -format runinfo \
| cut -d',' -f1 | tail -n +2 > SRR_all_samples.list

```

Then, the samples we are testing were manually selected and added into `SRR_subset.list`, and used to upload the fastq files using the command below.

```bash
parallel --jobs 4 'fastq-dump --split-files --origfmt --gzip {}' :::: SRR_subset.list

```

The metadata for the dataset was added using `scp` and is labelled as `SraRunTable_PRJNA1240347.csv`. This metadata contains all samples in case we want to check other samples for the analysis as well.

#### Step 2: download reference genome and annotation files

To download the annotation file for the reference, the following instruction was used:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
--2025-10-30 20:00:15--  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz

```

#### Step 3: indexing the reference genome using STAR

Parameter for --sjdbOverhang was selected as SRA run table indicated an average spot length of 76, and (76 \* 2) - 1 = 151.

```bash
STAR --runMode genomeGenerate --genomeDir STARIndex --genomeFastaFiles GCF_000001635.27_GRCm39_genomic.fna  --sjdbGTFfile GCF_000001635.27_GRCm39_genomic.gtf --sjdbOverhang 151 --runThreadN 8

```

#### Step 4: aligning the reference genome using STAR

```
STAR --
```
