import sys, os

# sys.path.insert(0, os.path.dirname(__file__))
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import numpy as np

class Classifications:
    def __init__(self):
        tp_df = pd.read_csv(snakemake.input.tp, sep="\t")
        tp_baseline_df = pd.read_csv(snakemake.input.tp_baseline, sep="\t")
        fp_df = pd.read_csv(snakemake.input.fp, sep="\t")
        fn_df = pd.read_csv(snakemake.input.fn, sep="\t")
        tp_query = len(tp_df["chromosome"])
        tp_truth = len(tp_baseline_df["chromosome"])
        fp = len(fp_df["chromosome"])
        fn = len(fn_df["chromosome"])
        
        self.tp_query = tp_query
        self.tp_truth = tp_truth
        self.fp = fp
        self.fn = fn

    def precision(self):
        p = self.tp_query + self.fp
        if p == 0:
            return 1.0
        p = float(self.tp_query) / float(p)
        return p

    def recall(self):
        t = self.tp_truth + self.fn
        if t == 0:
            return 1.0
        t = float(self.tp_truth) / float(t)
        return t


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
        a = float(self.tp_query) / float(a)
        return a

def collect_results():
    classifications_existence = Classifications()

    d = pd.DataFrame(
        {   "#variants_truth": [classifications_existence.tp_truth + classifications_existence.fn],
            "precision": [classifications_existence.precision()],
            "tp_query": [classifications_existence.tp_query],
            "fp": [classifications_existence.fp],
            "recall": [classifications_existence.recall()],
            "tp_truth": [classifications_existence.tp_truth],
            "fn": [classifications_existence.fn],
            "F*": [classifications_existence.fstar()],
        }
    )

    d =  d[
        [   "#variants_truth",
            "precision",
            "tp_query",
            "fp",
            "recall",
            "tp_truth",
            "fn",
            "genotype_mismatch_rate",
            "F*",
        ]
    ]

    return (d, None)

