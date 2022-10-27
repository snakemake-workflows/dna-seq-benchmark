import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


def get_cov_label(coverage):
    lower = snakemake.params.coverage_lower_bounds[coverage]
    bounds = [
        bound
        for bound in snakemake.params.coverage_lower_bounds.values()
        if bound > lower
    ]
    if bounds:
        upper = min(bounds)
        return f"{lower}..{upper}"
    else:
        return f"≥{lower}"


def load_data(f, coverage):
    d = pd.read_csv(f, sep="\t")
    d.insert(0, "coverage", get_cov_label(coverage))
    return d


if snakemake.input:
    report = pd.concat(
        load_data(f, cov) for cov, f in zip(snakemake.params.coverages, snakemake.input)
    )

    report.to_csv(snakemake.output[0], sep="\t", index=False)
else:
    pd.DataFrame({}).to_csv(snakemake.output[0], sep="\t")
