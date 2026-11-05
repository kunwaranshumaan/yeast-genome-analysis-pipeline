import csv
import subprocess
from collections import Counter
from pathlib import Path

files = {
    ("RAW", "Nuclear"): Path(
        "results/variants/annotation/full/raw_nuclear_annotated.vcf"
    ),
    ("TRIMMED", "Nuclear"): Path(
        "results/variants/annotation/full/trimmed_nuclear_annotated.vcf"
    ),
    ("RAW", "Mitochondrial"): Path(
        "results/variants/annotation/full/raw_mito_annotated.vcf"
    ),
    ("TRIMMED", "Mitochondrial"): Path(
        "results/variants/annotation/full/trimmed_mito_annotated.vcf"
    ),
}

output_dir = Path("results/variants/annotation/comparison")
output_dir.mkdir(parents=True, exist_ok=True)

summary = []
category_rows = []

for (dataset, genome), filename in files.items():

    if not filename.is_file():
        raise FileNotFoundError(f"Input VCF not found: {filename}")

    result = subprocess.run(
        ["bcftools", "query", "-f", "%INFO/BCSQ\n", str(filename)],
        capture_output=True,
        text=True,
        check=True
    )

    total = 0
    annotated = 0
    unannotated = 0

    record_counts = Counter()
    annotation_counts = Counter()

    for line_number, line in enumerate(result.stdout.splitlines(), start=1):
        total += 1
        bcsq = line.strip()

        if not bcsq or bcsq == ".":
            unannotated += 1
            continue

        consequences = []

        for item in bcsq.split(","):
            if "|" not in item:
                raise ValueError(
                    f"Malformed BCSQ annotation in {filename}, "
                    f"query line {line_number}: {item}"
                )

            consequence = item.split("|", 1)[0]

            if not consequence:
                raise ValueError(
                    f"Missing consequence in {filename}, "
                    f"query line {line_number}"
                )

            consequences.append(consequence)

        annotated += 1

        annotation_counts.update(consequences)
        record_counts.update(set(consequences))

    summary.append({
        "dataset": dataset,
        "genome": genome,
        "total_records": total,
        "annotated_records": annotated,
        "unannotated_records": unannotated,
        "annotated_percent": (
            round(100 * annotated / total, 2) if total else 0
        )
    })

    for consequence in sorted(
        set(record_counts) | set(annotation_counts)
    ):
        category_rows.append({
            "dataset": dataset,
            "genome": genome,
            "consequence": consequence,
            "records": record_counts[consequence],
            "annotations": annotation_counts[consequence]
        })

summary_fields = [
    "dataset",
    "genome",
    "total_records",
    "annotated_records",
    "unannotated_records",
    "annotated_percent"
]

category_fields = [
    "dataset",
    "genome",
    "consequence",
    "records",
    "annotations"
]

with open(
    output_dir / "raw_vs_trimmed_annotation_summary.csv",
    "w",
    newline=""
) as f:
    writer = csv.DictWriter(f, fieldnames=summary_fields)
    writer.writeheader()
    writer.writerows(summary)

with open(
    output_dir / "raw_vs_trimmed_consequence_categories.csv",
    "w",
    newline=""
) as f:
    writer = csv.DictWriter(f, fieldnames=category_fields)
    writer.writeheader()
    writer.writerows(category_rows)

print("Comparison files generated successfully.")

print("\nDataset summary:")
for row in summary:
    print(row)

print("\nSaved files:")
print(output_dir / "raw_vs_trimmed_annotation_summary.csv")
print(output_dir / "raw_vs_trimmed_consequence_categories.csv")
