# End-to-End Whole-Genome Sequencing Analysis Pipeline of *Saccharomyces cerevisiae* YJM451

## Project Overview

An end-to-end bioinformatics workflow for analysing paired-end whole-genome sequencing (WGS) data from *Saccharomyces cerevisiae* strain YJM451 against the S288C reference genome.

The project covers raw-read quality control, adapter trimming, read alignment, alignment quality assessment, duplicate marking, variant calling, RAW-versus-TRIMMED comparisons, and functional consequence annotation.

The workflow combines established command-line bioinformatics tools with Python-based analysis and visualization.

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
| Initial read pairs | 4,336,904 |

Dataset: [NCBI SRA – SRR800775](https://www.ncbi.nlm.nih.gov/sra/SRR800775)

## Reference Genome

| Parameter | Details |
|---|---|
| Reference strain | S288C |
| Assembly accession | GCF_000146045.2 |
| Assembly name | R64 |
| Nuclear chromosomes | 16 |
| Mitochondrial genome | 1 |
| Total reference length | 12,157,105 bp |

Reference: [NCBI Assembly – GCF_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/)

## Analysis Workflow

### 1. Raw sequencing quality control

Tools: FastQC and MultiQC.

Quality assessment included:
- Per-base sequence quality
- GC-content distribution
- Sequence duplication
- Adapter contamination
- Sequence length distribution
- Overrepresented sequences

### 2. Read trimming

Tool: fastp.

Paired-end reads were processed using adapter detection and trimming.

- Initial read pairs: 4,336,904
- Retained read pairs: 3,947,176
- Retention: 91.03%
- Removed read pairs: approximately 8.97%

Post-trimming quality was assessed using FastQC and MultiQC.

### 3. Reference preparation and alignment

Tool: Bowtie2.

Reads were aligned against the S288C reference genome using the `--very-sensitive` preset.

| Alignment metric | RAW | TRIMMED |
|---|---:|---:|
| Overall alignment rate | 94.80% | 95.91% |
| Properly paired reads | 97.15% | 97.96% |

### 4. Alignment processing and quality assessment

Tool: SAMtools.

The workflow included:
- BAM conversion and coordinate sorting
- BAM indexing and integrity checks
- Alignment statistics
- Chromosome-wise coverage assessment
- Duplicate marking using `fixmate` and `markdup`

Duplicate marking identified:
- RAW: 154,800 duplicate records (approximately 1.88%)
- TRIMMED: 120,515 duplicate records (approximately 1.59%)

Duplicate marking flags duplicates; it does not remove them.

### 5. Coverage analysis

Chromosome-wise sequencing coverage and depth were assessed using SAMtools and Python.

The nuclear chromosomes showed broad reference coverage. The mitochondrial genome had substantially higher sequencing depth than the nuclear chromosomes.

Chromosome XII also showed elevated depth and a high-depth region requiring cautious interpretation.

### 6. Variant calling

Tool: BCFtools.

SNPs and indels were called from duplicate-marked RAW and TRIMMED alignments.

The corrected mpileup workflow used:
- Mapping quality threshold: 20
- Base quality threshold: 20
- Maximum per-file depth for mpileup: 10,000
- Maximum per-file depth for indel calling: 10,000

Variant records were subsequently normalized, compared and filtered.

### 7. RAW-versus-TRIMMED variant comparison

The two call sets were compared to investigate the effect of read trimming on variant detection.

The normalized call-set comparison showed:
- Shared allele records: 61,319
- RAW-only records: 191
- TRIMMED-only records: 115
- Union: 61,625
- Jaccard similarity: 99.5034%

Jaccard similarity describes overlap between the normalized call sets. It is not a measure of variant-calling accuracy, sensitivity or specificity.

### 8. Functional consequence annotation

Tool: BCFtools `csq`.

The NCBI GFF annotation corresponding to the S288C reference assembly was used.

Separate genetic codes were applied:
- Nuclear genome: genetic code 1
- Mitochondrial genome: genetic code 3

The annotation workflow was run separately for RAW and TRIMMED datasets.

#### Annotation summary

| Dataset | Genome | Total records | Annotated | Unannotated |
|---|---|---:|---:|---:|
| RAW | Nuclear | 54,856 | 32,982 | 21,874 |
| TRIMMED | Nuclear | 54,811 | 32,964 | 21,847 |
| RAW | Mitochondrial | 505 | 86 | 419 |
| TRIMMED | Mitochondrial | 503 | 85 | 418 |

#### Major nuclear consequences

| Consequence | RAW records | TRIMMED records |
|---|---:|---:|
| Synonymous | 21,203 | 21,194 |
| Missense | 10,719 | 10,711 |
| Intron | 513 | 513 |
| Frameshift | 130 | 129 |
| In-frame insertion | 121 | 121 |
| In-frame deletion | 111 | 110 |
| Stop gained | 28 | 28 |

A variant may have multiple consequence annotations because of overlapping transcripts or features. Record counts and annotation counts are therefore reported separately in the detailed CSV output.

## Selected Investigation: Chromosome XII rDNA Region

Two candidate substitutions were investigated in the chromosome XII rDNA region:

- NC_001144.5:460103 T>C
- NC_001144.5:460152 T>G

Both positions fall within the annotated NTS1-2 non-transcribed spacer.

Read-level analyses identified differences in mapping-quality support between the candidates. The second candidate showed a more substantial reduction in alternate-allele support after applying stricter mapping-quality thresholds.

These are unresolved repeat-associated candidate variants. They should not be treated as confirmed single-copy mutations because the yeast rDNA region is repetitive and may be incompletely represented in the reference assembly.

## Software and Tools

- SRA Toolkit
- FastQC
- MultiQC
- fastp
- Bowtie2
- SAMtools
- BCFtools
- Python
- Pandas
- Matplotlib
- pysam
- Git and GitHub

Python dependencies are listed in `requirements.txt`.

## Repository Structure

The GitHub repository contains analysis scripts, documentation, and selected quality-control and comparison outputs.

```text
yeast-genome-analysis-pipeline/
├── docs/
│   └── quality_control_interpretation.md
├── results/
│   ├── figures/
│   └── quality_control/
│       ├── alignment/
│       ├── fastp/
│       ├── multiqc_raw/
│       ├── multiqc_trimmed/
│       └── raw_fastqc/
├── scripts/
│   ├── analyze_discordant_variants.py
│   ├── analyze_reference.py
│   ├── check_allele_mapping_quality.py
│   ├── check_allele_read_support.py
│   ├── check_paired_allele_agreement.py
│   ├── compare_variant_consequences.py
│   ├── investigate_discordant_depth.py
│   ├── plot_chromosome_coverage.py
│   └── plot_raw_vs_trimmed_coverage.py
├── requirements.txt
├── README.md
└── LICENSE
```

**Local-only analysis directories:** `data/`, `results/alignment/`, and `results/variants/` contain raw sequencing data, reference files, intermediate alignment files, and variant outputs. These are excluded from GitHub using `.gitignore`.

## Reproducibility

### Dataset and reference

- **Sequencing dataset:** SRR800775
- **Organism:** *Saccharomyces cerevisiae* YJM451
- **Sequencing platform:** Illumina HiSeq 2000, paired-end reads
- **Reference genome:** *S. cerevisiae* S288C, NCBI accession GCF_000146045.2_R64
- **Reference annotation:** Corresponding NCBI GFF3 annotation

Raw sequencing data, processed FASTQ files, reference files, BAM files, and VCF files are not included in this GitHub repository. Obtain the original dataset from the NCBI Sequence Read Archive and the reference genome and annotation from NCBI before reproducing the complete analysis.

### Python environment

Clone the repository:

```bash
git clone https://github.com/kunwaranshumaan/yeast-genome-analysis-pipeline.git
cd yeast-genome-analysis-pipeline
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### External software

The following command-line tools were used and must be installed separately:

| Tool | Version | Purpose |
|---|---:|---|
| Bowtie2 | 2.5.4 | Read alignment |
| SAMtools | 1.24 | BAM processing and coverage analysis |
| BCFtools | 1.24 | Variant calling and functional annotation |
| FastQC | 0.12.1 | Read quality assessment |
| MultiQC | 1.35 | Aggregated quality-control reports |
| fastp | 1.0.1 | Adapter trimming and read filtering |
| SRA Toolkit | 3.2.1 | SRA data retrieval |

Python 3.14.4 was used. Python package versions are pinned in `requirements.txt`.

### Analysis scripts

Run scripts from the repository root because input and output paths are defined relative to that directory. The scripts are analysis utilities, not a single automated workflow.

| Script | Required input or dependency |
|---|---|
| `analyze_reference.py` | Reference FASTA |
| `plot_chromosome_coverage.py` | Alignment coverage summary TXT |
| `plot_raw_vs_trimmed_coverage.py` | Raw-versus-trimmed coverage CSV |
| `analyze_discordant_variants.py` | Raw-only and trimmed-only discordant variant TSV files |
| `investigate_discordant_depth.py` | Discordant variant TSV files |
| `check_allele_mapping_quality.py` | Duplicate-marked RAW and TRIMMED BAM files; pysam |
| `check_allele_read_support.py` | Duplicate-marked RAW and TRIMMED BAM files; pysam |
| `check_paired_allele_agreement.py` | Duplicate-marked RAW and TRIMMED BAM files; pysam |
| `compare_variant_consequences.py` | Four annotated VCF files; BCFtools |

Example execution from the repository root:

```bash
python scripts/analyze_reference.py
python scripts/plot_chromosome_coverage.py
python scripts/plot_raw_vs_trimmed_coverage.py
```

The remaining analysis scripts require their corresponding intermediate BAM, TSV, or annotated VCF files to have been generated first. Their input paths are currently predefined in the scripts.

### Workflow overview

1. Retrieve paired-end sequencing reads from SRA.
2. Assess raw read quality using FastQC and MultiQC.
3. Perform adapter trimming and read filtering using fastp.
4. Assess post-trimming read quality.
5. Align raw and trimmed reads independently against the S288C reference using Bowtie2.
6. Process, sort, index, and mark duplicates in BAM files using SAMtools.
7. Calculate alignment statistics and chromosome-level coverage.
8. Call variants using BCFtools and normalize and compare the resulting variant records.
9. Annotate variants using the reference genome and NCBI GFF3 annotation.
10. Compare variant consequences and investigate selected discordant and repeat-associated candidates.

For complete reproduction, the reference files, intermediate results, directory structure, and software versions must match the analysis requirements. The GitHub repository provides selected scripts, reports, and summary outputs rather than all intermediate data.

## Limitations

- Functional consequences are computational predictions, not experimental validation.
- The variant-calling workflow used the default diploid assumption; sample ploidy was not independently established.
- Unphased genotype handling does not establish biological haplotype phase.
- Aggregate consequence counts do not identify individual variants gained or lost after trimming.
- The repetitive rDNA region requires cautious interpretation.

## Author

**Kunwar Anshumaan Singh**

M.Sc. Biotechnology | Bioinformatics and Molecular Biology

GitHub: [kunwaranshumaan](https://github.com/kunwaranshumaan)
