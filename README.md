
# Modular Nextflow Pipeline for SRR → Assembly → Amplicons → Pairwise Matrix

This repository provides a **modular, Singularity-based Nextflow (DSL2)** pipeline for:

1. Downloading raw reads from NCBI SRA (SRR accessions)
2. Assembling reads with **Shovill**
3. Predicting in‑silico amplicons using **EMBOSS primersearch**
4. Generating a **pairwise difference matrix** from amplicon FASTA files

---

## Pipeline stages and data flow

```
SRR list
   │
   ▼
fastqs/          (download)
   │
   ▼
assemblies/      (shovill)
   │
   ▼
amplicons/       (primersearch)
   │
   ▼
pairwise matrix  (CSV)
```

---

## Directory structure

When all stages are run, outputs are written under a single `--outdir`.  

```
results/
├── fastqs/
│   ├── SRR123_1.fastq.gz
│   ├── SRR123_2.fastq.gz
│   └── ...
├── assemblies/
│   │   └── SRR123_assembled.fasta
│   └── ...
├── primersearch/
│   ├── amplicon/
│   │   ├── SRR123_assembled_extractedAmplicons.fasta
│   │   ├── SRR123_assembled_not_matched_primers.txt
│   ├── SRR123_assembled.fasta
│   ├── SRR123.json
│   ├── SRR123.ps
│   └── ...
├── matrix_rows/
│   │   └── SRR123_assembled_extractedAmplicons.txt
│   └── ...
└── matrix/
    └── pairwise_diff_matrix.csv
```

---



## Common usage patterns

### 1. Full pipeline: SRR → matrix

```bash
nextflow run main.nf \
   --srr_list <SRR list file (one SRR per row)> \
   --outdir $PWD/<output folder> (need to use absolute path, or prefix it with $PWD\
   --primers <primer file (3-column primersearch format)> \ 
    --output_csv pairwise_test_50srrs.csv 
```

### 2. Download-only (raw FASTQ files from SRA)  

```bash
nextflow run main.nf -entry download_only \
    --srr_list srr_ids.txt \
    --outdir results
```

### 3. Assembly-only (skip SRA download)   

```bash
nextflow run main.nf -entry assembly_only \
    --reads /path/to/fastqs \
    --outdir results
```

### 4. Download + assembly (da)

```bash
nextflow run main.nf -entry da \
    --srr_list srr_ids.txt \
    --outdir results
```

### 5. Primersearch-only  

```bash
nextflow run main.nf -entry primersearch_only \
    --reads /path/to/assemblies \
    --primers primers.txt \
    --outdir results \
```

### 6. Pairwise-only

```bash
nextflow run main.nf -entry pairwise_only \
    --amplicon_folder $PWD/amplicons (need to use absolute path, or prefix it with $PWD \
    --primers primers.txt \
    --output_csv pairwise.csv \
    --outdir output_folder
```

### 7. Primersearch + pairwise matrix (pp)

```bash
nextflow run main.nf -entry pp \
    --reads assemblies \
    --primers primers.txt \
    --outdir $PWD/results  (need to use absolute path, or prefix it with $PWD\
    --output_csv pairwise.csv
```




---

## Repository layout

```
.
├── main.nf
├── nextflow.config
├── modules/
│   ├── download_reads.nf
│   ├── assemble_shovill.nf
│   ├── primersearch_amplicons.nf
│   └── pairwise_matrix.nf
├── bin/
│   ├── parse_primersearch.py
│   ├── fasta_to_json.py
│   ├── pairwise_compare.py
│   └── merge_rows.py
└── README.md
```

---

## Notices

### Public Domain Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest. 

### License Standard Notice

### Privacy Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

### Contributing Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov). 
