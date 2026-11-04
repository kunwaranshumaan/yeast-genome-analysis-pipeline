import pandas as pd
import matplotlib.pyplot as plt

# Input and output files
input_file = "results/quality_control/alignment/SRR800775_coverage.txt"
output_file = "results/figures/chromosome_coverage.png"

# Read SAMtools coverage output
columns = [
    "chromosome", "start", "end", "numreads",
    "covbases", "coverage", "meandepth",
    "meanbaseq", "meanmapq"
]

df = pd.read_csv(
    input_file,
    sep="\t",
    comment="#",
    names=columns
)

# Separate nuclear chromosomes from mitochondrial genome
nuclear = df[df["chromosome"] != "NC_001224.1"].copy()

# Extract chromosome number for plotting
nuclear["chromosome_number"] = range(1, len(nuclear) + 1)

# Plot mean depth
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
