import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import re

def parse_vaf(v):
    if isinstance(v, str) and v.startswith("{"):
        v = v.split(f"'{sample}': ")[1].split(",")[0].strip("}")
    if isinstance(v, str) and "%" in v:
        # Handle entries like "np.str_('100%')" or plain "100%"
        match = re.search(r"([\d.]+)%", v)
        if match:
            return float(match.group(1)) / 100
    return float(v) if isinstance(v, str) else v

def check_samplename_existence(df, sample_name):
    exists = True
    if "SAMPLE" in df.columns:
        if sample_name not in df["SAMPLE"].unique():
            print(f"Sample name '{sample_name}' not found in SAMPLE column.\nSample column is not filtered and removed for this callset.", file=sys.stderr)
            exists = False
    else:
        print("SAMPLE column not found in the input table.\nSample column is not filtered and removed for this callset", file=sys.stderr)
        exists = False
    return exists

df = pd.read_csv(snakemake.input.table, sep="\t")

if check_samplename_existence(df, snakemake.params.tumor_sample_name):
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

    df["vaf"] = df["vaf"].apply(parse_vaf)

# Remove SAMPLE column
if "SAMPLE" in df.columns:
    df = df.drop(columns=["SAMPLE"])

df.to_csv(snakemake.output.renamed_table, sep="\t", index=False)
