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

    df["POS"] = pd.to_numeric(df["POS"])
    df["DP"] = pd.to_numeric(df["DP"])
    df["QUAL"] = pd.to_numeric(df["QUAL"])

    def calculate_vaf(ad):
        try:
            values = [int(x) for x in str(ad).split(",")]
            if len(values) == 2 and sum(values) > 0:
                return values[1] / sum(values)
        except (ValueError, TypeError):
            pass
        return None

    df["VAF"] = df["AD"].apply(calculate_vaf)

    print("\n" + "=" * 65)
    print(dataset)
    print("=" * 65)

    print("\nVariants by chromosome:")
    print(df.groupby("CHROM").size().sort_values(
        ascending=False
    ).head(10).to_string())

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

