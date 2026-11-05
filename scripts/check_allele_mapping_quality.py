from pathlib import Path
import pysam

bam_files = {
    "RAW": Path("results/alignment/duplicate_marked/raw/SRR800775_raw_markdup.bam"),
    "TRIMMED": Path("results/alignment/duplicate_marked/trimmed/SRR800775_trimmed_markdup.bam")
}

candidates = [
    ("NC_001144.5", 460103, "T", "C"),
    ("NC_001144.5", 460152, "T", "G")
]

thresholds = [0, 20, 30, 40]

for dataset, bam_path in bam_files.items():

    if not bam_path.is_file():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        if not bam.has_index():
            raise FileNotFoundError(f"BAM index missing: {bam_path}.bai")

        print(f"\n========== {dataset} ==========")

        for chrom, pos, ref, alt in candidates:

            counts = {
                threshold: {ref: 0, alt: 0}
                for threshold in thresholds
            }

            for read in bam.fetch(chrom, pos - 1, pos):

                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_duplicate
                    or read.is_qcfail
                    or read.mapping_quality == 255
                ):
                    continue

                for query_pos, ref_pos in read.get_aligned_pairs():

                    if ref_pos == pos - 1 and query_pos is not None:

                        base = read.query_sequence[query_pos].upper()
                        quality = read.query_qualities[query_pos]

                        if quality < 20 or base not in (ref, alt):
                            break

                        for threshold in thresholds:
                            if read.mapping_quality >= threshold:
                                counts[threshold][base] += 1

                        break

            print(f"\n{chrom}:{pos} {ref}>{alt}")

            for threshold in thresholds:
                r = counts[threshold][ref]
                a = counts[threshold][alt]
                total = r + a
                fraction = 100 * a / total if total else 0

                print(
                    f"MAPQ >= {threshold}: "
                    f"REF={r}, ALT={a}, "
                    f"ALT fraction={fraction:.2f}%"
                )
