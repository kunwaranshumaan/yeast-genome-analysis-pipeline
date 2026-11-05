import subprocess
import csv
from collections import Counter
from pathlib import Path

files = {
    ("RAW", "Nuclear"): "results/variants/annotation/full/raw_nuclear_annotated.vcf",
    ("TRIMMED", "Nuclear"): "results/variants/annotation/full/trimmed_nuclear_annotated.vcf",
    ("RAW", "Mitochondrial"): "results/variants/annotation/full/raw_mito_annotated.vcf",
    ("TRIMMED", "Mitochondrial"): "results/variants/annotation/full/trimmed_mito_annotated.vcf",
}

output_dir = Path("results/variants/annotation/comparison")
output_dir.mkdir(parents=True, exist_ok=True)

summary = []

for (dataset, genome), filename in files.items():
    result = subprocess.run(
        ["bcftools", "query", "-f", "%INFO/BCSQ\n", filename],
        capture_output=True, text=True, check=True
    )

    total = annotated = unannotated = 0
    record_counts = Counter()
    annotation_counts = Counter()

    for line in result.stdout.splitlines():
        total += 1
        bcsq = line.strip()

        if not bcsq or bcsq == ".":
            unannotated += 1
            continue

        annotated += 1
        record_categories = set()

        for item in bcsq.split(","):
            consequence = item.split("|")[0]
            annotation_counts[consequence] += 1
            record_categories.add(consequence)

        for consequence in record_categories:
            record_counts[consequence] += 1

    summary.append({
        "dataset": dataset,
        "genome": genome,
        "total_records": total,
        "annotated_records": annotated,
        "unannotated_records": unannotated,
        "annotated_percent": round(100 * annotated / total, 2) if total else 0
    })

    for consequence in sorted(set(record_counts) | set(annotation_counts)):
        summary_row = {
            "dataset": dataset,
            "genome": genome,
            "consequence": consequence,
            "records": record_counts[consequence],
            "annotations": annotation_counts[consequence]
        }
        category_rows.append(summary_row) if "category_rows" in globals() else None

# Rebuild category rows in a separate pass for clarity
category_rows = []

for (dataset, genome), filename in files.items():
    result = subprocess.run(
        ["bcftools", "query", "-f", "%INFO/BCSQ\n", filename],
        capture_output=True, text=True, check=True
    )

    record_counts = Counter()
    annotation_counts = Counter()

    for line in result.stdout.splitlines():
        bcsq = line.strip()
        if not bcsq or bcsq == ".":
            continue

        record_categories = set()
        for item in bcsq.split(","):
            consequence = item.split("|")[0]
            annotation_counts[consequence] += 1
            record_categories.add(consequence)

        record_counts.update(record_categories)

    for consequence in sorted(set(record_counts) | set(annotation_counts)):
        category_rows.append({
            "dataset": dataset,
            "genome": genome,
            "consequence": consequence,
            "records": record_counts[consequence],
            "annotations": annotation_counts[consequence]
        })

with open(output_dir / "raw_vs_trimmed_annotation_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary[0].keys())
    writer.writeheader()
    writer.writerows(summary)

with open(output_dir / "raw_vs_trimmed_consequence_categories.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=category_rows[0].keys())
    writer.writeheader()
    writer.writerows(category_rows)

print("Comparison files generated successfully.")
print("\nDataset summary:")
for row in summary:
    print(row)

print("\nSaved files:")
print(output_dir / "raw_vs_trimmed_annotation_summary.csv")
print(output_dir / "raw_vs_trimmed_consequence_categories.csv")
