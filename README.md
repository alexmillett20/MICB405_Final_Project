# MICB405_Final_Project

##### Overview

This repository contains the documentation for the MICB 405 final project. This repository contains
the scripts used in R to generate the differential expression results, and the README.md will contain the steps
for data processing.

**Dataset**: 

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




