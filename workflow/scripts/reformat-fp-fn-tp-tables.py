import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

df = pd.read_csv(snakemake.input.table, sep="\t")

# Remove all normal samples
if "SAMPLE" in df.columns:
    if df["SAMPLE"].nunique(dropna=False) == 1:
        pass
    else:
        df = df[df["SAMPLE"] == str(snakemake.params.tumor_sample_name)]

# Rename columns with the expression handed to the script
df = df.rename(columns=snakemake.params.expression)

# If vaf column contains per-sample dicts (from FORMAT fields),
# extract only the tumor sample value.
if "vaf" in df.columns:
    sample = str(snakemake.params.tumor_sample_name)
    df["vaf"] = df["vaf"].apply(
        lambda v: float(v.split(f"'{sample}': ")[1].split(",")[0].strip("}"))
        if isinstance(v, str) and v.startswith("{")
        else v
    )

df.to_csv(snakemake.output.renamed_table, sep="\t", index=False)