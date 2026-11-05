import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

input_file = Path(
    "results/quality_control/alignment/raw_vs_trimmed_coverage.csv"
)

output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file)

# Remove mitochondrial chromosome from nuclear comparison
nuclear = df[df["rname"] != "NC_001224.1"].copy()

# Short chromosome labels
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
