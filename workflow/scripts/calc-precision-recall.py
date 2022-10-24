from collections import defaultdict
import sys, os

# sys.path.insert(0, os.path.dirname(__file__))
sys.stderr = open(snakemake.log[0], "w")

from abc import ABC, abstractmethod
from enum import Enum
import pandas as pd
import pysam

from common.happy_report import CompareExactGenotype, CompareExistence, Class


class Classifications:
    def __init__(self, comparator):
        self.tp_query = 0
        self.tp_truth = 0
        self.fn = 0
        self.fp = 0
        self.comparator = comparator
        self.visited = defaultdict(set)

    def register(self, record):
        for c in self.comparator.classify(record):
            if c.variant in self.visited[c.cls]:
                # skip if exactly this has been reported before, 
                # see workflow/scripts/common/happy_report.py for explanation
                continue
            if c.cls is Class.TP_truth:
                self.tp_truth += 1
            elif c.cls is Class.TP_query:
                self.tp_query += 1
            elif c.cls is Class.FN:
                self.fn += 1
            elif c.cls is Class.FP:
                self.fp += 1
            else:
                assert False, "unexpected case"

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


def collect_results(vartype):
    classifications_exact = Classifications(CompareExactGenotype(vartype))
    classifications_existence = Classifications(CompareExistence(vartype))

    for record in pysam.VariantFile(snakemake.input[0]):
        classifications_exact.register(record)
        classifications_existence.register(record)

    vartype = "indels" if vartype == "INDEL" else "snvs"

    mismatched_genotype = (
        classifications_existence.tp_query - classifications_exact.tp_query
    )
    if classifications_existence.tp_query > 0:
        mismatched_genotype_rate = (
            mismatched_genotype / classifications_existence.tp_query
        )
    else:
        mismatched_genotype_rate = 0.0

    d = pd.DataFrame(
        {
            "precision": [classifications_existence.precision()],
            "tp_query": [classifications_existence.tp_query],
            "fp": [classifications_existence.fp],
            "recall": [classifications_existence.recall()],
            "tp_truth": [classifications_existence.tp_truth],
            "fn": [classifications_existence.fn],
            "genotype_mismatch_rate": [mismatched_genotype_rate],
        }
    )
    return d[
        [
            "precision",
            "tp_query",
            "fp",
            "recall",
            "tp_truth",
            "fn",
            "genotype_mismatch_rate",
        ]
    ]


vartype = "SNP" if snakemake.wildcards.vartype == "snvs" else "INDEL"

collect_results(vartype).to_csv(snakemake.output[0], sep="\t", index=False)
