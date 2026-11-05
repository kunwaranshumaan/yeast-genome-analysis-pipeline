from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


input_file = Path(
    "results/quality_control/alignment/SRR800775_coverage.txt"
)
output_file = Path("results/figures/chromosome_coverage.png")

columns = [
    "chromosome", "start", "end", "numreads",
    "covbases", "coverage", "meandepth",
    "meanbaseq", "meanmapq"
]

if not input_file.is_file():
    raise FileNotFoundError(f"Coverage file not found: {input_file}")

df = pd.read_csv(
    input_file,
    sep="\t",
    comment="#",
    names=columns
)

if df.empty:
    raise ValueError("Coverage file contains no data.")

if df["chromosome"].duplicated().any():
    raise ValueError("Duplicate chromosome entries found.")

nuclear = df[df["chromosome"] != "NC_001224.1"].copy()

if len(nuclear) != 16:
    raise ValueError(
        f"Expected 16 nuclear chromosomes, found {len(nuclear)}."
    )

if not pd.to_numeric(nuclear["meandepth"], errors="coerce").notna().all():
    raise ValueError("Invalid mean-depth values detected.")

nuclear["chromosome_number"] = range(1, len(nuclear) + 1)

output_file.parent.mkdir(parents=True, exist_ok=True)

plt.figure(figsize=(12, 6))

plt.bar(
    nuclear["chromosome_number"],
    nuclear["meandepth"]
)

plt.xlabel("S288C Chromosome")
plt.ylabel("Mean sequencing depth (×)")
plt.title("Chromosome-wise Sequencing Coverage of Yeast Strain YJM451")
plt.xticks(range(1, 17))
plt.grid(axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()

print(f"Coverage plot saved to: {output_file}")
