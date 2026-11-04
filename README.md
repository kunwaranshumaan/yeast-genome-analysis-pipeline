# End-to-End NGS Genome Analysis Pipeline of *Saccharomyces cerevisiae* YJM451

## Project Overview

This project aims to perform an end-to-end analysis of whole-genome sequencing (WGS) data from *Saccharomyces cerevisiae* strain YJM451 using the S288C reference genome.

The workflow includes raw sequencing quality assessment, reference genome preparation, paired-end read alignment, alignment quality assessment, and chromosome-wise sequencing coverage visualization.

The project is being developed as a reproducible bioinformatics workflow using command-line tools and Python.

## Dataset Information

| Parameter | Details |
|---|---|
| Organism | *Saccharomyces cerevisiae* |
| Sample strain | YJM451 |
| SRA accession | SRR800775 |
| Sequencing platform | Illumina HiSeq 2000 |
| Sequencing strategy | Whole-genome sequencing |
| Library layout | Paired-end |
| Read length | 101 bp |
| Number of read pairs | 4,336,904 |

**Dataset source:** [NCBI Sequence Read Archive – SRR800775](https://www.ncbi.nlm.nih.gov/sra/SRR800775)

## Reference Genome

- Organism: *Saccharomyces cerevisiae* S288C
- Assembly accession: `GCF_000146045.2`
- Assembly name: R64
- Reference sequences: 16 nuclear chromosomes and one mitochondrial genome.
- Total reference sequence length: 12,157,105 bp.

**Reference source:** [NCBI Assembly – GCF_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/)

## Analysis Workflow

### 1. Raw Read Quality Control

Tools:
- FastQC
- MultiQC

Raw paired-end FASTQ files were assessed for:
- Per-base sequence quality
- GC content distribution
- Sequence duplication
- Adapter contamination
- Sequence length distribution
- Overrepresented sequences

### 2. Reference Genome Preparation

The S288C reference genome was downloaded from NCBI and indexed using Bowtie2.

A reference analysis script was developed in Python to calculate chromosome-wise sequence lengths.

### 3. Paired-End Read Alignment

Tool: Bowtie2

YJM451 sequencing reads were aligned against the S288C reference genome using the sensitive alignment preset.

**Overall alignment rate: 94.80%.**

### 4. Alignment Processing and Quality Assessment

Tool: SAMtools

The alignment workflow included:
- SAM-to-BAM conversion
- Coordinate sorting
- BAM indexing
- Alignment statistics
- Chromosome-wise mapping statistics
- Coverage assessment

### 5. Chromosome-wise Coverage Visualization

Python libraries:
- Pandas
- Matplotlib

Mean sequencing depth was calculated and visualized across the 16 nuclear chromosomes.

## Preliminary Results

| Parameter | Result |
|---|---:|
| Overall Bowtie2 alignment rate | 94.80% |
| Properly paired reads (SAMtools) | 97.15% |
| Average read length | 101 bp |
| Average base quality in BAM | 34.4 |
| Average insert size | 278 bp |
| Nuclear chromosome coverage | 96.35–99.96% |
| Nuclear mean sequencing depth | 53.40–88.92× |
| Mitochondrial mean sequencing depth | 1001.85× |

The mitochondrial genome showed substantially higher sequencing depth than the nuclear chromosomes. Chromosome XII also showed elevated mean depth and warrants further investigation.

## Repository Structure

```text
yeast-genome-analysis-pipeline/
│
├── docs/
│   └── quality_control_interpretation.md
│
├── results/
│   ├── figures/
│   │   └── chromosome_coverage.png
│   └── quality_control/
│       ├── raw_fastqc/
│       ├── multiqc_raw/
│       └── alignment/
│
├── scripts/
│   ├── analyze_reference.py
│   └── plot_chromosome_coverage.py
│
├── requirements.txt
├── README.md
└── .gitignore
```

Raw sequencing data, reference genome files, Bowtie2 indexes, and large alignment files are excluded from Git version control. They can be downloaded independently using the accessions and reference assembly listed above.

## Software and Tools

- SRA Toolkit
- FastQC
- MultiQC
- Bowtie2
- SAMtools
- Python
- Pandas
- Matplotlib
- Git and GitHub

Python dependencies are listed in `requirements.txt`.

## Reproducibility

Clone the repository:

```bash
git clone https://github.com/kunwaranshumaan/yeast-genome-analysis-pipeline.git
cd yeast-genome-analysis-pipeline
```

Create and activate a Python virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install Python dependencies:

```bash
python -m pip install -r requirements.txt
```

## Upcoming Work

- Further investigation of GC-content distribution.
- Investigation of elevated chromosome XII and mitochondrial coverage.
- Read trimming and post-trimming quality assessment.
- Variant calling using BCFtools.
- Variant filtering and annotation.
- Biological interpretation of genomic differences between YJM451 and S288C.

## Author

**Kunwar Anshumaan Singh**

M.Sc. Biotechnology | Bioinformatics and Molecular Biology

GitHub: [kunwaranshumaan](https://github.com/kunwaranshumaan)

---

*This repository is an ongoing bioinformatics project. Results and workflow components will be updated as additional analysis stages are completed.*
