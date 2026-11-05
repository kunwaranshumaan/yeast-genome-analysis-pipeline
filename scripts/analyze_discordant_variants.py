import pandas as pd
from pathlib import Path

base = Path(
    "results/variants/depth_corrected/comparison/discordance_analysis"
)

files = {
    "RAW_ONLY": base / "raw_only.tsv",
    "TRIMMED_ONLY": base / "trimmed_only.tsv"
}

for dataset, filepath in files.items():

    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AD"]
    )

    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["QUAL"] = pd.to_numeric(df["QUAL"], errors="coerce")
    df["DP"] = pd.to_numeric(df["DP"], errors="coerce")

    df["REF"] = df["REF"].str.upper()
    df["ALT"] = df["ALT"].str.upper()

    def classify_variant(row):
        ref = row["REF"]
        alt = row["ALT"]

        if len(ref) == 1 and len(alt) == 1:
            return "SNP"
        elif len(ref) < len(alt):
            return "Insertion"
        elif len(ref) > len(alt):
            return "Deletion"
        return "Complex"

    df["TYPE"] = df.apply(classify_variant, axis=1)

    ad_lengths = df["AD"].fillna("").apply(
        lambda x: len(x.split(",")) if x else 0
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
    ad_sums = df["AD"].fillna("").apply(
        lambda x: sum(int(v) for v in x.split(",") if v.isdigit())
        if x else None
    )

    comparable = ad_sums.notna() & df["DP"].notna()

    print(
        "Records where AD sum differs from DP:",
        (ad_sums[comparable] != df.loc[comparable, "DP"]).sum()
    )

    print("\nMedian DP by variant type:")
    print(df.groupby("TYPE")["DP"].median().to_string())

    print("\nMedian QUAL by variant type:")
    print(df.groupby("TYPE")["QUAL"].median().to_string())

