from pathlib import Path

fasta_file = Path(
    "data/reference/GCF_000146045.2_R64_genomic.fna"
)

sequences = {}
current_name = None

with fasta_file.open() as file:
    for line in file:
        line = line.strip()

        if line.startswith(">"):
            current_name = line[1:]
            sequences[current_name] = 0

        elif line:
            sequences[current_name] += len(line)

print("Reference Genome Summary")
print("-" * 60)

total_length = 0

for name, length in sequences.items():
    print(f"{name}: {length:,} bp")
    total_length += length

print("-" * 60)
print(f"Total genome size: {total_length:,} bp")
print(f"Total sequences: {len(sequences)}")
