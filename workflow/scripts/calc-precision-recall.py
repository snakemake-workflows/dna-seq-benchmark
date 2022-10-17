import sys

sys.stderr = open(snakemake.log[0], "w")

from abc import ABC, abstractmethod
from enum import Enum
import pandas as pd
import pysam


class Class(Enum):
    FP = 0
    TP_query = 1
    TP_truth = 2
    FN = 3


def get_gt(record, sample):
    gt = record.samples[sample]["GT"]
    if gt[0] is None:
        return None
    else:
        return sorted(gt)


def get_decision(record, sample):
    bd = record.samples[sample]["BD"]
    if bd is None:
        return None
    else:
        return bd


def get_vartypes(record):
    return record.samples[0]["BVT"], record.samples[1]["BVT"]


def same_alt_alleles(gt1, gt2):
    def alt_alleles(gt):
        if gt is not None:
            return set(a for a in gt if a > 0)
        else:
            return set()

    return alt_alleles(gt1) == alt_alleles(gt2)


class GenotypeCompare(ABC):
    def __init__(self, vartype):
        self.vartype = vartype

    @abstractmethod
    def classify(self, record):
        ...


class CompareExactGenotype(GenotypeCompare):
    def classify(self, record):
        truth_decision, query_decision = get_decision(record, 0), get_decision(
            record, 1
        )
        truth_vartype, query_vartype = get_vartypes(record)
        if truth_decision == "TP" and truth_vartype == self.vartype:
            yield Class.TP_truth
        if query_decision == "TP" and query_vartype == self.vartype:
            yield Class.TP_query
        if truth_decision == "FN" and truth_vartype == self.vartype:
            yield Class.FN
        if query_decision == "FP" and query_vartype == self.vartype:
            yield Class.FP


class CompareExistence(CompareExactGenotype):
    def classify(self, record):
        truth_gt, query_gt = get_gt(record, 0), get_gt(record, 1)

        for classification in super().classify(record):
            if classification == Class.FP:
                # determine whether this is only FP because the genotype does not match
                if same_alt_alleles(truth_gt, query_gt):
                    yield Class.TP_query
                else:
                    yield classification
            elif classification == Class.FN:
                # determine whether this is only FN because the genotype does not match
                if same_alt_alleles(truth_gt, query_gt):
                    yield Class.TP_truth
                else:
                    yield classification
            else:
                yield classification


class Classifications:
    def __init__(self, comparator):
        self.tp_query = 0
        self.tp_truth = 0
        self.fn = 0
        self.fp = 0
        self.comparator = comparator

    def register(self, record):
        for c in self.comparator.classify(record):
            if c is Class.TP_truth:
                self.tp_truth += 1
            elif c is Class.TP_query:
                self.tp_query += 1
            elif c is Class.FN:
                self.fn += 1
            elif c is Class.FP:
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

    # if vartype == "INDEL":
    #     print(
    #         sorted(classifications_exact.visited),
    #         file=open(snakemake.output[0] + "visited-dump.py", "w"),
    #     )

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

    return pd.DataFrame(
        {
            "vartype": [vartype],
            "precision": [classifications_existence.precision()],
            "tp_query": [classifications_existence.tp_query],
            "fp": [classifications_existence.fp],
            "recall": [classifications_existence.recall()],
            "tp_truth": [classifications_existence.tp_truth],
            "fn": [classifications_existence.fn],
            "mismatched_genotype_rate": [mismatched_genotype_rate],
        }
    )


data = pd.concat([collect_results("SNP"), collect_results("INDEL")], axis="rows")

data.to_csv(snakemake.output[0], sep="\t", index=False)
