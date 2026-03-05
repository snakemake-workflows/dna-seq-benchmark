import sys, os

# sys.path.insert(0, os.path.dirname(__file__))
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import numpy as np


def vaf_to_bin(vaf):
    """Convert a VAF value to a bin index (0-9).

    Bins correspond to VAF ranges: [0, 0.1), [0.1, 0.2), ..., [0.9, 1.0].
    """
    if isinstance(vaf, str):
        vaf = float(vaf.replace("%", "")) / 100 if "%" in vaf else float(vaf)
    return max(0, min(9, int(float(vaf) * 10) - 1))


def bin_vafs(series):
    """Bin a series of VAF values into 10 bins and return counts per bin."""
    counts = np.zeros(10, dtype=np.uint64)
    for vaf in series.dropna():
        counts[vaf_to_bin(vaf)] += 1
    return counts


class Classifications:
    def __init__(self, tp_df, tp_baseline_df, fp_df, fn_df, vaf_stratify=False):
        self.tp_query = len(tp_df)
        self.tp_truth = len(tp_baseline_df)
        self.fp = len(fp_df)
        self.fn = len(fn_df)

        self.stratify_by_vaf = vaf_stratify
        if vaf_stratify:
            self.tp_query_vaf = bin_vafs(tp_df["vaf"])
            self.tp_truth_vaf = bin_vafs(tp_baseline_df["vaf"])
            self.fp_vaf = bin_vafs(fp_df["vaf"])
            self.fn_vaf = bin_vafs(fn_df["vaf"])

    def precision(self):
        p = self.tp_query + self.fp
        if p == 0:
            return 1.0
        return float(self.tp_query) / float(p)

    def recall(self):
        t = self.tp_truth + self.fn
        if t == 0:
            return 1.0
        return float(self.tp_truth) / float(t)

    def fstar(self):
        """Calculate F* score,
        see https://link.springer.com/article/10.1007/s10994-021-05964-1.

        This is an interpretable alternative to the F-measure.
        Proportion of correct predictions among all predictions and missed variants.
        Or in other words, the probability that a variant taken from the union of
        prediction and truth is correctly predicted.
        It is a monotonic transformation of the F-measure.
        """
        a = self.tp_query + self.fn + self.fp
        if a == 0:
            return 1.0
        return float(self.tp_query) / float(a)

    def _safe_divide_vaf(self, numerator, denominator):
        """Element-wise division with 1.0 default where denominator is zero."""
        result = np.ones(10, dtype=np.float64)
        mask = denominator > 0
        result[mask] = numerator[mask] / denominator[mask]
        return result

    def precision_vaf(self):
        denom = self.tp_query_vaf.astype(np.float64) + self.fp_vaf
        return self._safe_divide_vaf(self.tp_query_vaf, denom)

    def recall_vaf(self):
        denom = self.tp_truth_vaf.astype(np.float64) + self.fn_vaf
        return self._safe_divide_vaf(self.tp_truth_vaf, denom)

    def fstar_vaf(self):
        denom = self.tp_query_vaf.astype(np.float64) + self.fn_vaf + self.fp_vaf
        return self._safe_divide_vaf(self.tp_query_vaf, denom)


def collect_results():
    tp_df = pd.read_csv(snakemake.input.tp, sep="\t")
    tp_baseline_df = pd.read_csv(snakemake.input.tp_baseline, sep="\t")
    fp_df = pd.read_csv(snakemake.input.fp, sep="\t")
    fn_df = pd.read_csv(snakemake.input.fn, sep="\t")

    vaf_status = snakemake.params.vaf_status
    has_vaf = vaf_status and all(
        "vaf" in df.columns for df in [tp_df, tp_baseline_df, fp_df, fn_df]
    )

    classifications = Classifications(
        tp_df, tp_baseline_df, fp_df, fn_df, vaf_stratify=has_vaf
    )

    d = pd.DataFrame(
        {
            "#variants_truth": [classifications.tp_truth + classifications.fn],
            "precision": [classifications.precision()],
            "tp_query": [classifications.tp_query],
            "fp": [classifications.fp],
            "recall": [classifications.recall()],
            "tp_truth": [classifications.tp_truth],
            "fn": [classifications.fn],
            "F*": [classifications.fstar()],
        }
    )

    if has_vaf:
        vafs = [0.1 * x for x in range(1, 11)]
        d_vaf = pd.DataFrame(
            {
                "vaf": vafs,
                "#variants_truth": classifications.tp_truth_vaf + classifications.fn_vaf,
                "precision": classifications.precision_vaf(),
                "tp_query": classifications.tp_query_vaf,
                "fp": classifications.fp_vaf,
                "recall": classifications.recall_vaf(),
                "tp_truth": classifications.tp_truth_vaf,
                "fn": classifications.fn_vaf,
                "F*": classifications.fstar_vaf(),
            }
        )
    else:
        d_vaf = None

    return (d, d_vaf)


assert snakemake.wildcards.vartype in ["snvs", "indels"]

results, results_vaf = collect_results()
if snakemake.wildcards.mode == "base":
    results.to_csv(snakemake.output[0], sep="\t", index=False)
else:
    if results_vaf is not None:
        results_vaf.to_csv(snakemake.output[0], sep="\t", index=False)
    else:
        # if benchmark has stratified VAF results, but this callset not,
        # then also write unstratified results
        results.to_csv(snakemake.output[0], sep="\t", index=False)

