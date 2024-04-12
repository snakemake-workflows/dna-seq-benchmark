from collections import defaultdict
import sys, os

# sys.path.insert(0, os.path.dirname(__file__))
sys.stderr = open(snakemake.log[0], "w")

from abc import ABC, abstractmethod
from enum import Enum
import pandas as pd
import numpy as np
import pysam

from common.classification import CompareExactGenotype, CompareExistence, Class


class Classifications:
    def __init__(self, comparator, vaf_fields):
        if vaf_fields[0] is None or vaf_fields[1] is None:  
            self.tp_query = 0
            self.tp_truth = 0
            self.fn = 0
            self.fp = 0
            self.comparator = comparator
        else:
            self.tp_query = np.array([]) # 0-100%
            self.tp_truth = 0
            self.fn = 0
            self.fp = 0

    def increment_counter(self, counter, vaf):
        if vaf is None:
            counter += 1
        else:
            # 10 equally sized bins
            bin = int(vaf*10)
            counter[bin] += 1

    def register(self, record):
        for c in self.comparator.classify(record):
            # TODO: depending on case, fetch VAF from truth or query record (FP: from query record, field configurable by callset (e.g. FORMAT/AF, INFO/AF, ...)
            # for truth record, field configurable by benchmark preset (same syntax as above)
            # increment counters for bins, bins given to constructor as list of tuples or some numpy equivalent.
            # Default: None. If no VAF field given for either truth or callset, don't bin at all. 
            if c.cls is Class.TP_truth:
                vaf = truth.samples[record.name]["AF"] # still needs to be implemented
                self.increment_counter(self.tp_truth, vaf)
                # self.tp_truth += 1
            elif c.cls is Class.TP_query:
                vaf = truth.samples[record.name]["AF"] # still needs to be implemented
                self.increment_counter(self.tp_query, vaf)
                # self.tp_query += 1
            elif c.cls is Class.FN:
                vaf = truth.samples[record.name]["AF"] # still needs to be implemented
                self.increment_counter(self.fn, vaf)
                # self.fn += 1
            elif c.cls is Class.FP:
                vaf = record.samples[0]["AF"] # only works for FORMAT field
                self.increment_counter(self.fp, vaf)
                # self.fp += 1
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

    for record in pysam.VariantFile(snakemake.input.calls):
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

snakemake.params.vaf_fields

assert snakemake.wildcards.vartype in ["snvs", "indels"]
vartype = "SNV" if snakemake.wildcards.vartype == "snvs" else "INDEL"

collect_results(vartype).to_csv(snakemake.output[0], sep="\t", index=False)
