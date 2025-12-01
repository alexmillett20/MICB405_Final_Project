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
plots/
reads_per_gene/
scripts/
.gitignore # git-specific file. Not needed for analysis
.lintr
MICB405_Final_Project.Rproj # Generated when making R scripts in a project folder. Not needed for analysis
README.md
gene_GO_mapping.tsv
srr_sample_mapping.txt
```

* **DEA_outputs/**: This folder contains the output dds objects that we generated for each of the comparisons (ctrl vs IL13, ctrl vs IL4, IL13 vs IL4). There are also intermediate dds objects that we generated prior to outlier removal. Anything labeled with 'final' is what we used for analysis, as these excluded outliers.
* **aligned/**: This folder contains all the output files from STAR alignment (besides the aligned*sorted.bam files as these were too big to upload onto github). The ReadsPerGene.out.tab files were moved to `reads_per_gene/` for DEA.
* **fastqc-reports/**: This folder contains all the fastqc.html reports for all samples.
* **plots/**: This folder contains all the plots generated from intermediate steps and for the final report. Additional folders are found within this folder that contain the plots, organised by plot type (e.g. heatmap, go plots etc)
* **reads_per_gene/**: This folder contains all the ReadsPerGene.out.tab files that were renamed according to the original SRR accessions, in `srr_sample_mapping.txt`.
* **scripts/**: This folder contains all the scripts that were used for DEA, GOEA, and plot generation
* **gene_GO_mapping.tsv**: .tsv file with the GO annotations labelled for \_Mus musculus* that was used for GOEA. Script used to generate is available in `scripts/geneIDtoGO.py`. This was not used for figures in the final report, but was tested prior to using the annFUN.org function. 
* **srr_sample_mapping.txt**: this file contains the SRR accession IDs and the corresponding sample condition that was used for renaming the ReadsPerGene.out.tab files

### **Dataset**:

#### Step 1: download fastq files

To download the fastq files for the dataset, the following instructions were used:

```bash
esearch -db sra -query PRJNA1240347 \
| efetch -format runinfo \
| cut -d',' -f1 | tail -n +2 > SRR_all_samples.list

```

Then, the samples we are testing were manually selected and added into `SRR_subset.list`, and used to upload the fastq files onto the course server using the command below. `SRR_subset.list` contains all samples from IL-13, IL-4, and control conditions only. 

```bash
parallel --jobs 4 'fastq-dump --split-files --origfmt --gzip {}' :::: SRR_subset.list

```

The metadata for the dataset was added using `scp` and is labelled as `SraRunTable_PRJNA1240347.csv`, which contains the metadata for all samples in the dataset (which includes samples that we did not include for our analyses). 

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

Using `run_star_alignment.sh`, aligned each of the read pairs using STAR.

```
STAR \
        --runThreadN 8 \
        --genomeDir "/work/data/STARIndex" \
        --readFilesIn "/work/data/raw-data/Read_1.fastq.gz" "/work/data/raw-data/Read_2.fastq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "/work/data/aligned/Read_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode GeneCounts
```

#### Step 5: perform DEA using read counts

Since we specified `--quantMode GeneCounts`, output files labelled `*_ReadsPerGene.out.tab` were downloaded locally and uploaded to this github repository. Since samples were labelled with their access numbers, we manually relabelled the samples based on the metadata file from the SRA. These were double checked by team members to ensure there was no mistake in relabelling.

Under `scripts/DEA_script_DESeq2.R`, this is where the actual analysis occurred. We first performed DEA using all samples, to identify if there were any outliers. Then, we generated a PCA plot based on the rlog transformed dds object, and determined that control replicate 8 was an outlier. Next, we generated the results for the differential expression for each comparison of interest: control vs IL13, control vs IL4, and IL13 vs IL4. From this, we generated PCA plots as well as heatmaps based on the sample distances of the dds object, and identified control replicate 3 and IL13 replicate 6 as outliers. We reran the DEA using the finalised dataset without these outliers, and used these for all downstream analyses.

We kept the scripts and plots with these intermediate steps for reproducibility and transparency, and showcase the plots in our supplementary material. 

We determined differentially expressed genes as ones that satisfied a log2FC cut off of +/-1 and padj < 0.05. 

Next, we generated heatmaps based on the Z scores of the genes, volcano plots to identify the DEGs between the conditions, MA plots to visualise this with respect to normalized expression of the genes, and then performed GO enrichment analysis. All scripts for these analyses are found under `scripts/`. 
