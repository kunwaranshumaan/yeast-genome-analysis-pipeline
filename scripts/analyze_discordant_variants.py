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


def classify_variant(row):
    ref = row["REF"]
    alt = row["ALT"]

    if pd.isna(ref) or pd.isna(alt) or not ref or not alt:
        return "Invalid"

    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    elif len(ref) < len(alt):
        return "Insertion"
    elif len(ref) > len(alt):
        return "Deletion"

    return "Complex"


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

    if df.empty:
        raise ValueError(f"{dataset} contains no variant records.")

    if df.isna().any().any():
        raise ValueError(f"Missing values detected in {dataset}.")

    for column in ["POS", "QUAL", "DP"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")

        if df[column].isna().any():
            raise ValueError(f"Invalid numeric values in {dataset}: {column}")

    if (df["POS"] < 1).any() or (df["DP"] < 0).any():
        raise ValueError(f"Invalid position or depth values in {dataset}.")

    df["REF"] = df["REF"].str.upper()
    df["ALT"] = df["ALT"].str.upper()

    df["TYPE"] = df.apply(classify_variant, axis=1)

    if (df["TYPE"] == "Invalid").any():
        raise ValueError(f"Invalid REF/ALT alleles detected in {dataset}.")

    ad_lengths = df["AD"].apply(
        lambda x: len(x.split(","))
    )

    ad_sums = df["AD"].apply(
        lambda x: sum(int(v) for v in x.split(","))
    )

    print("\n" + "=" * 55)
    print(dataset)
    print("=" * 55)

    print("\nTotal variants:", len(df))

    print("\nVariant types:")
    print(df["TYPE"].value_counts().to_string())

    print("\nRead depth:")
    print(df["DP"].describe().round(2).to_string())

    print("\nQUAL:")
    print(df["QUAL"].describe().round(2).to_string())

    print("\nDepth categories:")
    print("DP < 5:", (df["DP"] < 5).sum())
    print("DP < 10:", (df["DP"] < 10).sum())
    print("DP < 20:", (df["DP"] < 20).sum())

    print("\nAD field structure:")
    print(ad_lengths.value_counts().sort_index().to_string())

    print("\nAD sum versus DP:")
    print(
        "Records where AD sum differs from DP:",
        (ad_sums != df["DP"]).sum()
    )

    print("\nMedian DP by variant type:")
    print(df.groupby("TYPE")["DP"].median().to_string())

    print("\nMedian QUAL by variant type:")
    print(df.groupby("TYPE")["QUAL"].median().to_string())
