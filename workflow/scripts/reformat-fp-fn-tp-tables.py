import sys

sys.stderr = open(snakemake.log[0], "w")

import ast
import pandas as pd

df = pd.read_csv(snakemake.input.table, sep="\t")
df = df.rename(columns=snakemake.params.expression)

# If vaf column contains per-sample dicts (from FORMAT fields),
# extract only the tumor sample value.
if "vaf" in df.columns:
    sample = str(snakemake.params.tumor_sample_name)
    df["vaf"] = df["vaf"].apply(
        lambda v: ast.literal_eval(v).get(sample)
        if isinstance(v, str) and v.startswith("{")
        else v
    )

df.to_csv(snakemake.output.renamed_table, sep="\t", index=False)