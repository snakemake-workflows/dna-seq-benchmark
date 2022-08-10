import pandas as pd
from scipy.cluster.hierarchy import ward, leaves_list
from scipy.spatial.distance import pdist


def read_data(f, callset):
    data = pd.read_csv(
        f,
        sep="\t",
        index_col=[
            "chromosome",
            "position",
            "ref_allele",
            "alt_allele",
            "true_genotype",
        ],
    )

    data.columns = [callset]
    return data


def get_idx_sorted_by_clustering(data):
    cluster_matrix = ward(pdist(~data.isna(), metric="hamming"))
    idx = leaves_list(cluster_matrix)
    return idx


data = pd.concat(
    [
        read_data(f, callset)
        for f, callset in zip(snakemake.input, snakemake.params.callsets)
    ],
    axis="columns",
)

data = data.loc[
    data.isna().sum(axis="columns").sort_values().index,
]

if not data.empty:
    idx_rows = get_idx_sorted_by_clustering(data)
    idx_cols = get_idx_sorted_by_clustering(data.T)
    data = data.iloc[idx_rows, idx_cols]

data.to_csv(snakemake.output[0], sep="\t")
