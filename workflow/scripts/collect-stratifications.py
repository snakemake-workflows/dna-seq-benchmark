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

    # TODO With separate files for SNVs and indels with e.g. STRELKA no predicted variants for the other type are expected
    # If later relevant, add annotation to the report
    # if (report["tp_truth"] == 0).all():
    #     raise ValueError(
    #         f"The callset {snakemake.wildcards.callset} does not predict any variant from the truth. "
    #         "This is likely a technical issue in the callset and should be checked before further evaluation."
    #     )

    report.to_csv(snakemake.output[0], sep="\t", index=False)
else:
    pd.DataFrame(
        {
            col: []
            for col in [
                "coverage",
                "precision",
                "tp_query",
                "fp",
                "recall",
                "tp_truth",
                "fn",
                "genotype_mismatch_rate",
                "F*",
            ]
        }
    ).to_csv(snakemake.output[0], sep="\t")
