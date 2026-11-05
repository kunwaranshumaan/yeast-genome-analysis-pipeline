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

            support = defaultdict(set)
            fragments = defaultdict(set)

            for read in bam.fetch(chrom, pos - 1, pos):

                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_duplicate
                    or read.is_qcfail
                    or read.mapping_quality < 30
                    or read.mapping_quality == 255
                ):
                    continue

                positions = read.get_reference_positions(
                    full_length=True
                )

                try:
                    index = positions.index(pos - 1)
                except ValueError:
                    continue

                base = read.query_sequence[index].upper()

                if read.query_qualities[index] < 20:
                    continue

                if base not in (ref, alt):
                    continue

                allele = base
                mate = "R1" if read.is_read1 else "R2"
                strand = "Reverse" if read.is_reverse else "Forward"

                support[(allele, mate, strand)].add(read.query_name)
                fragments[allele].add(read.query_name)

            print(f"\n{chrom}:{pos} {ref}>{alt}")

            for allele in (ref, alt):
                print(f"\n{allele} allele:")

                for mate in ("R1", "R2"):
                    for strand in ("Forward", "Reverse"):
                        n = len(support[(allele, mate, strand)])
                        print(f"  {mate}, {strand}: {n}")

                print(
                    f"  Unique read names supporting allele: "
                    f"{len(fragments[allele])}"
                )

            print(
                "\nUnique fragments with ALT support: "
                f"{len(fragments[alt])}"
            )

