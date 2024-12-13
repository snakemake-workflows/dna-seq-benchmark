import sys
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


def load_data(path, callset):
    d = pd.read_csv(path, sep="\t")
    d.insert(0, "callset", callset)
    return d


results = pd.concat(
    [
        load_data(f, callset)
        for f, callset in zip(snakemake.input.tables, snakemake.params.callsets)
    ],
    axis="rows",
)

def cov_key(cov_label):
    # return lower bound as integer for sorting
    if ".." in cov_label:
        return int(cov_label.split("..")[0])
    else:
        return int(cov_label[1:])



def sort_key(col):
    if col.name == "callset":
        return col
    if col.name == "coverage":
        return col.apply(cov_key)
    else:
        return col


results.sort_values(["callset", "coverage"], inplace=True, key=sort_key)
results["sort_index"] = results["coverage"].apply(cov_key)

results.to_csv(snakemake.output[0], sep="\t", index=False)
