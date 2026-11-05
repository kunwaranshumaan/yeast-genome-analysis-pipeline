import pysam
from collections import defaultdict

bam_files = {
    "RAW": "results/alignment/duplicate_marked/raw/SRR800775_raw_markdup.bam",
    "TRIMMED": "results/alignment/duplicate_marked/trimmed/SRR800775_trimmed_markdup.bam"
}

candidates = [
    ("NC_001144.5", 460103, "T", "C"),
    ("NC_001144.5", 460152, "T", "G")
]

for dataset, bam_path in bam_files.items():
    print(f"\n========== {dataset} ==========")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, pos, ref, alt in candidates:
            fragments = defaultdict(dict)

            for read in bam.fetch(chrom, pos - 1, pos):
                if (
                    read.is_unmapped or read.is_secondary
                    or read.is_supplementary or read.is_duplicate
                    or read.is_qcfail or read.mapping_quality < 30
                    or read.mapping_quality == 255
                ):
                    continue

                positions = read.get_reference_positions(full_length=True)

                try:
                    index = positions.index(pos - 1)
                except ValueError:
                    continue

                if read.query_qualities[index] < 20:
                    continue

                base = read.query_sequence[index].upper()

                if base in (ref, alt):
                    mate = "R1" if read.is_read1 else "R2"
                    fragments[read.query_name][mate] = base

            both = {"REF": 0, "ALT": 0, "DISCORDANT": 0}
            single = {"REF": 0, "ALT": 0}

            for observations in fragments.values():
                if "R1" in observations and "R2" in observations:
                    a, b = observations["R1"], observations["R2"]

                    if a == b == ref:
                        both["REF"] += 1
                    elif a == b == alt:
                        both["ALT"] += 1
                    else:
                        both["DISCORDANT"] += 1
                else:
                    allele = next(iter(observations.values()))
                    single["ALT" if allele == alt else "REF"] += 1

            print(f"\n{chrom}:{pos} {ref}>{alt}")
            print("Both mates cover the position:")
            print(f"  Both REF: {both['REF']}")
            print(f"  Both ALT: {both['ALT']}")
            print(f"  Discordant: {both['DISCORDANT']}")
            print("Only one mate covers the position:")
            print(f"  REF: {single['REF']}")
            print(f"  ALT: {single['ALT']}")
