# Raw Read Quality Control Report

## Project
Yeast Genome Analysis Pipeline

## Dataset Information
- Organism: Saccharomyces cerevisiae
- Strain: YJM451
- SRA Accession: SRR800775
- Sequencing Platform: Illumina HiSeq 2000
- Sequencing Strategy: Whole Genome Sequencing (WGS)
- Library Layout: Paired-end
- Read Length: 101 bp

## Quality Control Tools
- FastQC v0.12.1
- MultiQC v1.35

## Sequencing Statistics

| Parameter | Forward Reads | Reverse Reads |
|---|---:|---:|
| Total reads | 4,336,904 | 4,336,904 |
| Mean GC content | 36.0% | 36.0% |
| Duplication | 17.6% | 17.3% |

## FastQC Findings

Most quality-control modules passed for both sequencing files.

The Per Sequence GC Content module failed in both files, showing a dominant distribution around 38–40% and a secondary low-GC population.

Other findings:
- Per-base sequence quality: PASS
- Per-sequence quality scores: PASS
- Per-base sequence content: PASS
- Per-base N content: PASS
- Sequence length distribution: PASS
- Sequence duplication levels: PASS
- Overrepresented sequences: PASS
- Adapter content: PASS

## Interpretation

The sequencing reads generally demonstrate good base-call quality and minimal adapter contamination.

The unusual GC-content distribution requires further investigation. The observed secondary low-GC population cannot independently establish contamination.

## Conclusion

The raw sequencing dataset is suitable for further computational investigation. Additional reference-based alignment and, if warranted, taxonomic screening will help assess the unusual GC distribution.

## Reproducibility

Raw sequencing data: NCBI SRA SRR800775

Quality-control reports: results/quality_control/

Software versions are recorded above.
