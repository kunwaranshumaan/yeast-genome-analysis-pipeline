import pysam

bam_files = {
    "RAW": "results/alignment/duplicate_marked/raw/SRR800775_raw_markdup.bam",
    "TRIMMED": "results/alignment/duplicate_marked/trimmed/SRR800775_trimmed_markdup.bam"
}

candidates = [
    ("NC_001144.5", 460103, "T", "C"),
    ("NC_001144.5", 460152, "T", "G")
]

thresholds = [0, 20, 30, 40]

for dataset, bam_path in bam_files.items():
    print(f"\n========== {dataset} ==========")

    with pysam.AlignmentFile(bam_path, "rb") as bam:

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

                positions = read.get_reference_positions(full_length=True)

                try:
                    index = positions.index(pos - 1)
                except ValueError:
                    continue

                base = read.query_sequence[index].upper()
                quality = read.query_qualities[index]

                if quality < 20 or base not in (ref, alt):
                    continue

                for threshold in thresholds:
                    if read.mapping_quality >= threshold:
                        counts[threshold][base] += 1

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
