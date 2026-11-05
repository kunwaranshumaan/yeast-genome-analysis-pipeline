import pandas as pd
from pathlib import Path

base = Path(
    "results/variants/depth_corrected/comparison/discordance_analysis"
)

files = {
    "RAW_ONLY": base / "raw_only.tsv",
    "TRIMMED_ONLY": base / "trimmed_only.tsv"
}

columns = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AD"]


def calculate_vaf(ad):
    values = [int(x) for x in ad.split(",")]

    if len(values) != 2 or any(x < 0 for x in values):
        raise ValueError(f"Invalid AD value: {ad}")

    total = sum(values)

    if total == 0:
        return None

    return values[1] / total


for dataset, filepath in files.items():

    if not filepath.is_file():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=columns,
        dtype={"CHROM": str, "REF": str, "ALT": str, "AD": str}
    )

    if df.empty or df.isna().any().any():
        raise ValueError(f"Empty data or missing values in {dataset}.")

    for column in ["POS", "DP", "QUAL"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")

        if df[column].isna().any():
            raise ValueError(f"Invalid numeric values in {column}.")

    if (df["POS"] < 1).any() or (df["DP"] < 0).any():
        raise ValueError(f"Invalid position or depth in {dataset}.")

    df["VAF"] = df["AD"].apply(calculate_vaf)

    if df["VAF"].isna().any():
        raise ValueError(f"Zero or invalid allele depth in {dataset}.")

    print("\n" + "=" * 65)
    print(dataset)
    print("=" * 65)

    print("\nVariants by chromosome:")
    print(
        df.groupby("CHROM").size()
        .sort_values(ascending=False)
        .head(10).to_string()
    )

    print("\nTop 15 variants by read depth:")
    print(
        df.nlargest(15, "DP")[
            ["CHROM", "POS", "REF", "ALT", "DP", "AD", "VAF", "QUAL"]
        ].to_string(index=False)
    )

    print("\nVAF distribution:")
    print(df["VAF"].describe().round(3).to_string())

    print("\nVariants with DP >= 100:")
    print((df["DP"] >= 100).sum())
