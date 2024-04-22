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
            self.vaf_field_query = vaf_fields[0][0]
            self.vaf_field_name_query = vaf_fields[0][1]
            self.vaf_field_truth = vaf_fields[1][0]
            self.vaf_field_name_truth = vaf_fields[1][1]
            # arrays with 10 fields (VAF from 0% to 100%)
            self.tp_query = np.zeros(10)
            self.tp_truth = np.zeros(10)
            self.fn = np.zeros(10)
            self.fp = np.zeros(10)

    def increment_counter(self, counter, vaf):
        if self.stratify_by_vaf:
            # 10 equally sized bins
            bin = int(vaf*10) - 1
            counter[bin] += 1
        else:
            counter += 1

    def register(self, record, truth, query):
        for c in self.comparator.classify(record):
            # depending on case, fetch VAF from truth or query record (FP: from query record, field configurable by callset (e.g. FORMAT/AF, INFO/AF, ...)
            # for truth record, field configurable by benchmark preset (same syntax as above)
            # increment counters for bins, bins given to constructor as list of tuples or some numpy equivalent.
            # Default: None. If no VAF field given for either truth or callset, don't bin at all. 
            if c.cls is Class.TP_truth:
                if self.stratify_by_vaf:
                    r = list(truth.fetch(record.contig, record.start, record.stop))[0]
                    vaf = r.info[self.vaf_field_name_truth] if self.vaf_field_truth == "INFO" else r.format[self.vaf_field_name_truth]
                else:
                    vaf = None
                self.increment_counter(self.tp_truth, vaf)
            elif c.cls is Class.TP_query:
                if self.stratify_by_vaf:
                    r = list(truth.fetch(record.contig, record.start, record.stop))[0]
                    vaf = r.info[self.vaf_field_name_truth] if self.vaf_field_truth == "INFO" else r.format[self.vaf_field_name_truth]
                else:
                    vaf = None
                self.increment_counter(self.tp_query, vaf)
            elif c.cls is Class.FN:
                if self.stratify_by_vaf:
                    r = list(truth.fetch(record.contig, record.start, record.stop))[0]
                    vaf = r.info[self.vaf_field_name_truth] if self.vaf_field_truth == "INFO" else r.format[self.vaf_field_name_truth]
                else:
                    vaf = None
                self.increment_counter(self.fn, vaf)
            elif c.cls is Class.FP:
                if self.stratify_by_vaf:
                    r = list(query.fetch(record.contig, record.start, record.stop))[0]
                    vaf = r.info[self.vaf_field_name_query][0] if self.vaf_field_query == "INFO" else r.format[self.vaf_field_name_query][0]
                else:
                    vaf = None
                self.increment_counter(self.fp, vaf)
            else:
                assert False, "unexpected case"

    def precision(self):
        if self.stratify_by_vaf:
            p = self.tp_query + self.fp
            for (i,x) in enumerate(p):
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
            t = self.tp_truth + self.fn
            for (i,x) in enumerate(t):
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

    # stratification by VAF
    else:
        vafs = [0.1*x for x in range(1,11)]

        d = pd.DataFrame(
            {
                "vaf": vafs,
                "precision": classifications_existence.precision(),
                "tp_query": classifications_existence.tp_query,
                "fp": classifications_existence.fp,
                "recall": classifications_existence.recall(),
                "tp_truth": classifications_existence.tp_truth,
                "fn": classifications_existence.fn,
            }
        )

        return d[
            [
                "vaf",
                "precision",
                "tp_query",
                "fp",
                "recall",
                "tp_truth",
                "fn",
            ]
        ]
    

assert snakemake.wildcards.vartype in ["snvs", "indels"]
vartype = "SNV" if snakemake.wildcards.vartype == "snvs" else "INDEL"

collect_results(vartype).to_csv(snakemake.output[0], sep="\t", index=False)
