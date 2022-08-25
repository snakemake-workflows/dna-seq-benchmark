import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


def load_data(f, coverage):
    d = pd.read_csv(f)
    d.insert(0, "coverage", coverage)
    return d


if snakemake.input:
    report = pd.concat(
        load_data(f, cov) for cov, f in zip(snakemake.params.coverages, snakemake.input)
    )

    report.to_csv(snakemake.output[0], sep="\t", index=False)
else:
    pd.DataFrame({}).to_csv(snakemake.output[0], sep="\t")
