import pandas as pd


def load_data(path, callset):
    d = pd.read_csv(path, sep="\t")
    d.insert(0, "callset", callset)
    return d


results = pd.concat(
    [
        load_data(f, callset)
        for f, callset in zip(snakemake.input, snakemake.params.callsets)
    ],
    axis="rows",
)
results.sort_values("callset", inplace=True)

results.to_csv(snakemake.output[0], sep="\t", index=False)
