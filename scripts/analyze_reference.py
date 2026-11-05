from pathlib import Path

fasta_file = Path("data/reference/GCF_000146045.2_R64_genomic.fna")


def analyze_fasta(file_path):
    if not file_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {file_path}")

    sequences = {}
    current_name = None

    with file_path.open() as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_name = line[1:].split()[0]

                if current_name in sequences:
                    raise ValueError(f"Duplicate FASTA header: {current_name}")

                sequences[current_name] = 0

            elif current_name is None:
                raise ValueError(
                    f"Sequence found before FASTA header at line {line_number}"
                )

            else:
                sequences[current_name] += len(line)

    if not sequences:
        raise ValueError("No sequences found in the FASTA file.")

    return sequences


sequences = analyze_fasta(fasta_file)

print("Reference Genome Summary")
print("-" * 60)

total_length = 0

for name, length in sequences.items():
    print(f"{name}: {length:,} bp")
    total_length += length

print("-" * 60)
print(f"Total genome size: {total_length:,} bp")
print(f"Total sequences: {len(sequences)}")
