import os
import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from scipy.cluster.hierarchy import ward, leaves_list
from scipy.spatial.distance import pdist
from sklearn.feature_selection import chi2
from statsmodels.stats.multitest import fdrcorrection
import numpy as np


def collect_chromosomes(f):
    d = pd.read_csv(f, sep="\t", dtype=str, usecols=["chromosome"]).drop_duplicates()
    return d["chromosome"].tolist()


def read_data(f, callset, chromosome=None):
    print("reading", f, "...", file=sys.stderr)
    data = pd.read_csv(
        f,
        sep="\t",
        dtype=str,
    )

    if chromosome is not None:
        data = data.loc[data["chromosome"] == chromosome]
    print(data.head(), file=sys.stderr)

    data = data.set_index(
        [
            "vaf",
            "chromosome",
            "position",
            "ref_allele",
            "alt_allele",
            "true_genotype",
        ]
    )
    data.drop("class", axis="columns", inplace=True)

    assert (
        not data.index.duplicated().any()
    ), f"bug: not expecting any duplicates in FP/FN table {f}"

    data.columns = [callset]

    if snakemake.wildcards.classification == "fn":
        data.loc[:, callset] = "FN"
    return data


def get_idx_sorted_by_clustering(data):
    cluster_matrix = ward(pdist(~data.isna(), metric="hamming"))
    idx = leaves_list(cluster_matrix)
    return idx


chromosomes = sorted(
    {chrom for f in snakemake.input.tables for chrom in collect_chromosomes(f)}
)

if not chromosomes:
    chromosomes = [None]

n_written = 0

# process data for each chromosome separately and append to the same files
for i, chromosome in enumerate(chromosomes):
    if n_written > snakemake.params.max_entries:
        break

    data = pd.concat(
        [
            read_data(f, callset, chromosome)
            for f, callset in zip(snakemake.input.tables, snakemake.params.callsets)
        ],
        axis="columns",
    )

    data = data.loc[data.isna().sum(axis="columns").sort_values().index,]

    data = data.dropna(how="all")

    if data.shape[1] > 1 and data.shape[0] > 1:
        idx_rows = get_idx_sorted_by_clustering(data)
        idx_cols = get_idx_sorted_by_clustering(data.T)
        data = data.iloc[idx_rows, idx_cols]

    label_df = pd.DataFrame(snakemake.params.labels)

    def store(data, output, label_idx=None):
        _label_df = label_df
        if label_idx is not None:
            _label_df = label_df.iloc[[label_idx]]

        if i == 0:
            mode = "w"
            header = True
             # add labels
            index_cols = data.index.names
            cols = data.columns
            data = pd.concat([_label_df, data.reset_index()]).set_index(index_cols)
            # restore column order
            data = data[cols]
        else:
            mode = "a"
            header = False
        data.to_csv(output, sep="\t", mode=mode, header=header)

    store(data, snakemake.output.main)
    n_written += len(data)

    label_df.index = snakemake.params.label_names

    os.makedirs(snakemake.output.dependency_sorting, exist_ok=True)

    # for each label, calculate mutual information and sort FP/FN entries in descending order
    for label_idx, label in enumerate(snakemake.params.label_names):
        outdata = data
        if data.shape[1] > 1 and data.shape[0] > 1:
            # oe = OrdinalEncoder()
            # le = LabelEncoder()
            # feature matrix: genotypes, transposed into samples x features, replace NA with False (match)
            # and any genotype with True (mismatch with truth).
            feature_matrix = data.reset_index(drop=True).T.copy()
            feature_matrix[~pd.isna(feature_matrix)] = True
            feature_matrix[pd.isna(feature_matrix)] = False

            # target vector: label values, converted into factors
            target_vector = label_df.loc[label]

            # ignore any NA in the target vector and correspondingly remove the rows in the feature matrix
            not_na_target_vector = target_vector[~pd.isna(target_vector)]
            # clone data
            sorted_data = data.copy(deep=True)
            # sort by label
            sorted_target_vector = target_vector.sort_values()
            sorted_data = sorted_data[sorted_target_vector.index]
            if not not_na_target_vector.empty:

                feature_matrix = feature_matrix.loc[not_na_target_vector.index]

                # perfom chi² test against each feature
                _, pvals = chi2(feature_matrix, not_na_target_vector)
                # TODO: this only sorts per chromosome, not globally
                sorted_idx = np.argsort(pvals)

                _, fdr = fdrcorrection(
                    pvals, method="negcorr"
                )  # use Benjamini/Yekutieli as variants might be both positively or negatively correlated

                # add pvalue and FDR
                sorted_data.insert(0, "FDR dependency", np.around(fdr, 3))
                sorted_data.insert(0, "p-value dependency", np.around(pvals, 3))
                sorted_data = sorted_data.iloc[sorted_idx]

            else:
                sorted_data.insert(0, "FDR dependency", pd.NA)
                sorted_data.insert(0, "p-value dependency", pd.NA)

            outdata = sorted_data

            # only keep the rather significant entries (but be a bit more permissive than 0.05)
            outdata = outdata.loc[outdata["p-value dependency"] <= 0.25]

        outpath = os.path.join(snakemake.output.dependency_sorting, f"{label}.tsv")
        store(outdata, outpath, label_idx=label_idx)
