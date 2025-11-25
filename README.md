# MICB405_Final_Project

## Overview

This repository contains the documentation for the MICB 405 final project. This repository contains
the scripts used in R to generate the differential expression results, and the README.md will contain the steps
for data processing.

### Repository overview
The repository is organised as follows:
```
DEA_outputs/
aligned/
fastqc-reports/
old_reads_per_gene/ ############ (can remove?)
plots/
reads_per_gene/
scripts/
.gitignore # git-specific file. Not needed for analysis
.lintr
MICB405_Final_Project.Rproj # Generated when making R scripts in a project folder. Not needed for analysis
README.md
gene_GO_mapping.tsv
srr_sample_mapping.txt
srr_sample_name.txt ############### (can remove?)
```

**DEA_outputs/**: This folder contains the output dds objects that we generated for each of the comparisons (ctrl vs IL13, ctrl vs IL4, IL13 vs IL4). There are also intermediate dds objects that we generated prior to outlier removal. Anything labeled with 'final' is what we used for analysis, as these excluded outliers.
**aligned/**: This folder contains all the output files from STAR alignment (besides the aligned_sorted.bam files as these were too big to upload onto github). The ReadsPerGene.out.tab files were moved to `reads_per_gene/` for DEA.
**fastqc-reports/**: This folder contains all the fastqc.html reports for all samples.
**plots/**: This folder contains all the plots generated from intermediate steps and for the final report.
**reads_per_gene/**: This folder contains all the ReadsPerGene.out.tab files that were renamed according to the original SRR accessions, in `srr_sample_mapping.txt`. 
**scripts/**: This folder contains all the scripts that were used for DEA, GOEA, and plot generation
**gene_GO_mapping.tsv**: .tsv file with the GO annotations labelled for _Mus musculus_ that was used for GOEA. Script used to generate is available in `scripts/geneIDtoGO.py`
**srr_sample_mapping.txt**: this file contains the SRR accession IDs and the corresponding sample condition that was used for renaming the ReadsPerGene.out.tab files



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

#### Step 3: index the reference genome using STAR

Parameter for --sjdbOverhang was selected as SRA run table indicated an average spot length of 76, and (76 \* 2) - 1 = 151.
This was an incorrect interpretation of the read length, as the GEO files (and FASTQC report) showed each read was 38 bp, totalling 76 bp per sample. Though, this does **not** change any downstream analyses except for faster processing time, so we did not change the --sjdbOverhang prior to alignment. 

```bash
STAR --runMode genomeGenerate --genomeDir STARIndex --genomeFastaFiles GCF_000001635.27_GRCm39_genomic.fna  --sjdbGTFfile GCF_000001635.27_GRCm39_genomic.gtf --sjdbOverhang 151 --runThreadN 8

```

#### Step 4: align the reference genome using STAR

All default parameters were used according to the lecture.

```
STAR --genomeDir /work/STARIndex --readFilesIn /work/data/raw-data/*_1.fastq.gz /work/data/raw-data/*_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --runThreadN 8
```

#### Step 5: perform DEA using read counts

Since we specified `--quantMode GeneCounts`, output files labelled `*_ReadsPerGene.out.tab` were downloaded locally and uploaded to this github repository. Since samples were labelled with their access numbers, we manually relabelled the samples based on the metadata file from the SRA. These were double checked by team members to ensure there was no mistake in relabelling. 
