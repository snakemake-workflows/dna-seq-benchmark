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
        # if vaf information is given: stratify results by vaf and add addtional counters for stratified tp, fp, fn
        if vaf_fields[0] is not None and vaf_fields[1] is not None:
            self.stratify_by_vaf = True
            self.vaf_field_query = vaf_fields[0]["field"]
            self.vaf_field_name_query = vaf_fields[0]["name"]
            self.vaf_field_truth = vaf_fields[1]["field"]
            self.vaf_field_name_truth = vaf_fields[1]["name"]
            # arrays with 10 fields (VAF from 0% to 100%)
            self.tp_query_vaf = np.zeros(10, dtype=np.uint)
            self.tp_truth_vaf = np.zeros(10, dtype=np.uint)
            self.fn_vaf = np.zeros(10, dtype=np.uint)
            self.fp_vaf = np.zeros(10, dtype=np.uint)
        else:
            self.tp_query_vaf = None
            self.tp_truth_vaf = None
            self.fn_vaf = None
            self.fp_vaf = None
            self.stratify_by_vaf = False

        # always initalize the counters without vaf stratification
        self.tp_query = 0
        self.tp_truth = 0
        self.fn = 0
        self.fp = 0

    def increment_counter(self, current_record, other_record, counter, counter_vaf, fp=False):
        if self.stratify_by_vaf:
            r = list(other_record.fetch(current_record.contig, current_record.start, current_record.stop))
            if len(r) > 0:
                r = r[0]
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
                counter_vaf[bin] += 1

        counter += 1
        return (counter, counter_vaf)

    def register(self, record, truth, query):
        for c in self.comparator.classify(record):
            # depending on case, fetch VAF from truth or query record (FP: from query record, field configurable by callset (e.g. FORMAT/AF, INFO/AF, ...)
            # for truth record, field configurable by benchmark preset (same syntax as above)
            # increment counters for bins, bins given to constructor as list of tuples or some numpy equivalent.
            # Default: None. If no VAF field given for either truth or callset, don't bin at all.
            if c.cls is Class.TP_truth:
                (self.tp_truth, self.tp_truth_vaf) = self.increment_counter(record, truth, self.tp_truth, self.tp_truth_vaf)
            elif c.cls is Class.TP_query:
                (self.tp_query, self.tp_query_vaf) = self.increment_counter(record, truth, self.tp_query, self.tp_query_vaf)
            elif c.cls is Class.FN:
                (self.fn, self.fn_vaf) = self.increment_counter(record, truth, self.fn, self.fn_vaf)
            elif c.cls is Class.FP:
                (self.fp, self.fp_vaf) = self.increment_counter(record, query, self.fp, self.fp_vaf, True)
            else:
                assert False, "unexpected case"

    def precision(self):
        if self.stratify_by_vaf:
            p_vaf = self.tp_query_vaf.astype(np.float32) + self.fp_vaf
            for (i, x) in enumerate(p_vaf):
                if x == 0:
                    p_vaf[i] = 1.0
                else:
                    p_vaf[i] = self.tp_query_vaf[i] / p_vaf[i]
        else:
            p_vaf = None

        p = self.tp_query + self.fp
        if p == 0:
            return (1.0, p_vaf)
        p = float(self.tp_query) / float(p)
        return (p, p_vaf)

    def recall(self):
        if self.stratify_by_vaf:
            t_vaf = self.tp_truth_vaf.astype(np.float32) + self.fn_vaf
            for (i, x) in enumerate(t_vaf):
                if x == 0:
                    t_vaf[i] = 1.0
                else:
                    t_vaf[i] = self.tp_truth_vaf[i] / t_vaf[i]
        else:
            t_vaf = None

        t = self.tp_truth + self.fn
        if t == 0:
            return (1.0, t_vaf)
        t = float(self.tp_truth) / float(t)
        return (t, t_vaf)


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
            a_vaf = self.tp_query_vaf.astype(np.float32) + self.fn_vaf + self.fp_vaf
            for (i, x) in enumerate(a_vaf):
                if a_vaf[i] == 0:
                    a_vaf[i] = 1.0
                else:
                    a_vaf[i] = self.tp_query_vaf[i] / a_vaf[i]
        else:
            a_vaf = None

        a = self.tp_query + self.fn + self.fp
        if a == 0:
            return (1.0, a_vaf)
        a = float(self.tp_query) / float(a)
        return (a, a_vaf)



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

    if vaf_fields[0] is not None and vaf_fields[1] is not None:

        vafs = [0.1*x for x in range(1,11)]

        d = pd.DataFrame(
            {
                "vaf": vafs,
                "#variants_truth": classifications_existence.tp_truth_vaf + classifications_existence.fn_vaf,
                "precision": classifications_existence.precision()[1],
                "tp_query": classifications_existence.tp_query_vaf,
                "fp": classifications_existence.fp_vaf,
                "recall": classifications_existence.recall()[1],
                "tp_truth": classifications_existence.tp_truth_vaf,
                "fn": classifications_existence.fn_vaf,
                "F*": classifications_existence.fstar()[1],
            }
        )

        d_vaf =  d[
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
    else:
        d_vaf = None

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
            "precision": [classifications_existence.precision()[0]],
            "tp_query": [classifications_existence.tp_query],
            "fp": [classifications_existence.fp],
            "recall": [classifications_existence.recall()[0]],
            "tp_truth": [classifications_existence.tp_truth],
            "fn": [classifications_existence.fn],
            "genotype_mismatch_rate": [mismatched_genotype_rate],
            "F*": [classifications_existence.fstar()[0]],
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

    return (d, d_vaf)



assert snakemake.wildcards.vartype in ["snvs", "indels"]
vartype = "SNV" if snakemake.wildcards.vartype == "snvs" else "INDEL"

results, results_vaf = collect_results(vartype)
if snakemake.wildcards.mode == "base":
    print("base", results)
    results.to_csv(snakemake.output[0], sep="\t", index=False)
else:
    print("vaf:", results_vaf)
    results_vaf.to_csv(snakemake.output[0], sep="\t", index=False)
