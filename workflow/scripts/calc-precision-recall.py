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
        self.comparator = comparator
        if vaf_fields[0] is None or vaf_fields[1] is None:
            self.stratify_by_vaf = False
            self.tp_query = 0
            self.tp_truth = 0
            self.fn = 0
            self.fp = 0
        else:
            self.stratify_by_vaf = True
            self.vaf_field_query = vaf_fields[0]["field"]
            self.vaf_field_name_query = vaf_fields[0]["name"]
            self.vaf_field_truth = vaf_fields[1]["field"]
            self.vaf_field_name_truth = vaf_fields[1]["name"]
            # arrays with 10 fields (VAF from 0% to 100%)
            self.tp_query = np.zeros(10, dtype=np.uint)
            self.tp_truth = np.zeros(10, dtype=np.uint)
            self.fn = np.zeros(10, dtype=np.uint)
            self.fp = np.zeros(10, dtype=np.uint)

    def increment_counter(self, current_record, other_record, counter, fp=False):
        if self.stratify_by_vaf:
            r = list(other_record.fetch(current_record.contig, current_record.start, current_record.stop))[0]
            if fp:
                vaf = r.info[self.vaf_field_name_query] if self.vaf_field_query == "INFO" else r.samples[0][self.vaf_field_name_query]
            else:
                vaf = r.info[self.vaf_field_name_truth] if self.vaf_field_truth == "INFO" else r.samples[0][self.vaf_field_name_truth]
            if type(vaf) == tuple:
                vaf = vaf[0]
            if type(vaf) == str:
                vaf = float(vaf.replace("%", "")) / 100
            # 10 equally sized bins
            bin = max(0, int(vaf*10) - 1)
            counter[bin] += 1
        else:
            counter += 1
        return counter

    def register(self, record, truth, query):
        for c in self.comparator.classify(record):
            # depending on case, fetch VAF from truth or query record (FP: from query record, field configurable by callset (e.g. FORMAT/AF, INFO/AF, ...)
            # for truth record, field configurable by benchmark preset (same syntax as above)
            # increment counters for bins, bins given to constructor as list of tuples or some numpy equivalent.
            # Default: None. If no VAF field given for either truth or callset, don't bin at all.
            if c.cls is Class.TP_truth:
                self.tp_truth = self.increment_counter(record, truth, self.tp_truth)
            elif c.cls is Class.TP_query:
                self.tp_query = self.increment_counter(record, truth, self.tp_query)
            elif c.cls is Class.FN:
                self.fn = self.increment_counter(record, truth, self.fn)
            elif c.cls is Class.FP:
                self.fp = self.increment_counter(record, query, self.fp, True)
            else:
                assert False, "unexpected case"

    def precision(self):
        if self.stratify_by_vaf:
            p = self.tp_query.astype(np.float32) + self.fp
            for (i, x) in enumerate(p):
                if x == 0:
                    p[i] = 1.0
                else:
                    p[i] = self.tp_query[i] / p[i]
            return p
        else:
            p = self.tp_query + self.fp
            if p == 0:
                return 1.0
            return float(self.tp_query) / float(p)

    def recall(self):
        if self.stratify_by_vaf:
            t = self.tp_truth.astype(np.float32) + self.fn
            for (i, x) in enumerate(t):
                if x == 0:
                    t[i] = 1.0
                else:
                    t[i] = self.tp_truth[i] / t[i]
            return t
        else:
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
        if self.stratify_by_vaf:
            a = self.tp_query.astype(np.float32) + self.fn + self.fp
            for (i, x) in enumerate(a):
                if a[i] == 0:
                    a[i] = 1.0
                else:
                    a[i] = self.tp_query[i] / a[i]
            return a
        else:
            a = self.tp_query + self.fn + self.fp
            if a == 0:
               return 1.0
            return float(self.tp_query) / float(a)



def collect_results(vartype):
    vaf_fields = snakemake.params.vaf_fields
    classifications_exact = Classifications(CompareExactGenotype(vartype), vaf_fields)
    classifications_existence = Classifications(CompareExistence(vartype), vaf_fields)
    truth = pysam.VariantFile(snakemake.input.truth)
    query = pysam.VariantFile(snakemake.input.query)

    for record in pysam.VariantFile(snakemake.input.calls):
        classifications_exact.register(record, truth, query)
        classifications_existence.register(record, truth, query)

    vartype = "indels" if vartype == "INDEL" else "snvs"

    # no stratification by VAF
    if vaf_fields[0] is None or vaf_fields[1] is None:
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
            {   "#variants_truth": [classifications_existence.tp_truth + classifications_existence.fn],
                "precision": [classifications_existence.precision()],
                "tp_query": [classifications_existence.tp_query],
                "fp": [classifications_existence.fp],
                "recall": [classifications_existence.recall()],
                "tp_truth": [classifications_existence.tp_truth],
                "fn": [classifications_existence.fn],
                "genotype_mismatch_rate": [mismatched_genotype_rate],
                "F*": [classifications_existence.fstar()],
            }
        )

        return d[
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

    # stratification by VAF
    else:
        vafs = [0.1*x for x in range(1,11)]

        d = pd.DataFrame(
            {
                "vaf": vafs,
                "#variants_truth": classifications_existence.tp_truth + classifications_existence.fn,
                "precision": classifications_existence.precision(),
                "tp_query": classifications_existence.tp_query,
                "fp": classifications_existence.fp,
                "recall": classifications_existence.recall(),
                "tp_truth": classifications_existence.tp_truth,
                "fn": classifications_existence.fn,
                "F*": classifications_existence.fstar(),
            }
        )

        return d[
            [
                "vaf",
                "#variants_truth",
                "precision",
                "tp_query",
                "fp",
                "recall",
                "tp_truth",
                "fn",
                "F*",
            ]
        ]


assert snakemake.wildcards.vartype in ["snvs", "indels"]
vartype = "SNV" if snakemake.wildcards.vartype == "snvs" else "INDEL"

collect_results(vartype).to_csv(snakemake.output[0], sep="\t", index=False)
