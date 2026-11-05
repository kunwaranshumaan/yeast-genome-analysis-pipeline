from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


input_file = Path(
    "results/quality_control/alignment/raw_vs_trimmed_coverage.csv"
)

output_dir = Path("results/figures")

if not input_file.is_file():
    raise FileNotFoundError(f"Input CSV not found: {input_file}")

df = pd.read_csv(input_file)

required_columns = [
    "rname",
    "coverage_raw",
    "coverage_trimmed",
    "meandepth_raw",
    "meandepth_trimmed",
]

missing = set(required_columns) - set(df.columns)

if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

if df.empty:
    raise ValueError("Input CSV contains no data.")

if df["rname"].duplicated().any():
    raise ValueError("Duplicate chromosome identifiers found.")

expected_chromosomes = [
    "NC_001133.9", "NC_001134.8", "NC_001135.5",
    "NC_001136.10", "NC_001137.3", "NC_001138.5",
    "NC_001139.9", "NC_001140.6", "NC_001141.2",
    "NC_001142.9", "NC_001143.9", "NC_001144.5",
    "NC_001145.3", "NC_001146.8", "NC_001147.6",
    "NC_001148.4",
]

nuclear = df[df["rname"] != "NC_001224.1"].copy()

if nuclear["rname"].tolist() != expected_chromosomes:
    raise ValueError(
        "Nuclear chromosome identifiers or reference order do not match."
    )

for column in required_columns[1:]:
    nuclear[column] = pd.to_numeric(
        nuclear[column], errors="coerce"
    )

    if not nuclear[column].notna().all():
        raise ValueError(f"Invalid numeric values in {column}.")

output_dir.mkdir(parents=True, exist_ok=True)

nuclear["chromosome"] = [
    f"Chr {i}" for i in range(1, 17)
]

# Plot 1: Genome coverage
plt.figure(figsize=(13, 6))

plt.plot(
    nuclear["chromosome"],
    nuclear["coverage_raw"],
    marker="o",
    label="Raw reads"
)

plt.plot(
    nuclear["chromosome"],
    nuclear["coverage_trimmed"],
    marker="o",
    label="Trimmed reads"
)

plt.ylabel("Genome coverage (%)")
plt.xlabel("Chromosome")
plt.title("Raw vs Trimmed Read Alignment: Chromosome Coverage")
plt.xticks(rotation=45)
plt.ylim(90, 101)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(
    output_dir / "raw_vs_trimmed_chromosome_coverage.png",
    dpi=300
)
plt.close()

# Plot 2: Mean sequencing depth
plt.figure(figsize=(13, 6))

plt.plot(
    nuclear["chromosome"],
    nuclear["meandepth_raw"],
    marker="o",
    label="Raw reads"
)

plt.plot(
    nuclear["chromosome"],
    nuclear["meandepth_trimmed"],
    marker="o",
    label="Trimmed reads"
)

plt.ylabel("Mean sequencing depth (×)")
plt.xlabel("Chromosome")
plt.title("Raw vs Trimmed Read Alignment: Mean Depth")
plt.xticks(rotation=45)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(
    output_dir / "raw_vs_trimmed_chromosome_depth.png",
    dpi=300
)
plt.close()

print("Plots generated successfully:")
print(output_dir / "raw_vs_trimmed_chromosome_coverage.png")
print(output_dir / "raw_vs_trimmed_chromosome_depth.png")
