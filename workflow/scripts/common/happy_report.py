from abc import ABC, abstractmethod
from collections import namedtuple, defaultdict
from enum import Enum


class Class(Enum):
    FP = 0
    TP_query = 1
    TP_truth = 2
    FN = 3


Variant = namedtuple("Variant", "contig, pos, ref, alts")
Classification = namedtuple("Classification", "cls, variant")


def get_gt(record, sample):
    gt = record.samples[sample]["GT"]
    if gt[0] is None:
        return None
    else:
        return sorted(gt)


def get_variant(record, sample):
    gt = get_gt(record, sample)
    if gt is None:
        return None
    alts = {record.alleles[i] for i in gt}
    try:
        alts.remove(record.ref)
    except KeyError:
        # allele was not referenced by the genotype
        pass
    return Variant(record.contig, record.pos, record.ref, frozenset(alts))


def is_het(record, sample, variant):
    if len(variant.alts) > 1:
        return True
    gt = get_gt(record, sample)
    if len(set(gt)) != 1:
        return True


def get_decision(record, sample):
    bd = record.samples[sample]["BD"]
    if bd is None:
        return None
    else:
        return bd


def get_vartypes(record):
    ref_allele = record.alleles[0]
    alts_truth = alt_alleles(record, 0)
    alts_query = alt_alleles(record, 1)
    is_snv = lambda alt: len(alt) == 1 and len(ref_allele) == 1
    def vartype(alts):
        _is_snv = set(map(is_snv, alts))
        if not _is_snv:
            return None
        # we expect different types
        #assert len(_is_snv) == 1, f"{record.pos}, {_is_snv}"
        return "SNV" if all(_is_snv) else "INDEL"
    
    return vartype(alts_truth), vartype(alts_query)


def alt_alleles(record, sample):
    gt = get_gt(record, sample)
    if gt is not None:
        return set(record.alleles[a] for a in gt if a > 0)
    else:
        return set()


class GenotypeCompare(ABC):
    def __init__(self, vartype=None):
        self.vartype = vartype

    @abstractmethod
    def classify(self, record):
        ...

    def is_valid_vartype(self, vartype):
        return self.vartype is None or vartype == self.vartype


class CompareExactGenotype(GenotypeCompare):
    def classify(self, record):
        truth_decision, query_decision = get_decision(record, 0), get_decision(
            record, 1
        )
        truth_vartype, query_vartype = get_vartypes(record)

        truth_variant, query_variant = get_variant(record, 0), get_variant(record, 1)

        if truth_decision == "TP" and self.is_valid_vartype(truth_vartype):
            yield Classification(Class.TP_truth, truth_variant)
        if query_decision == "TP" and self.is_valid_vartype(query_vartype):
            yield Classification(Class.TP_query, query_variant)
        if truth_decision == "FN" and self.is_valid_vartype(truth_vartype):
            yield Classification(Class.FN, truth_variant)
        if query_decision == "FP" and self.is_valid_vartype(query_vartype):
            yield Classification(Class.FP, query_variant)


class CompareExistence(CompareExactGenotype):
    def __init__(self, vartype=None):
        super().__init__(vartype=vartype)
        self.visited = defaultdict(set)

    def classify(self, record):
        # Inferring TP, FP, and FN is complicated from the happy output when trying
        # to not count twice and ignoring genotypes. Consider the following example:
        # 12      122278216       .       G       T,GAA   .       .       BS=122278216    GT:BD:BK:QQ:BI:BVT:BLT  0/2:FN:am:.:i1_5:INDEL:het      1/0:FP:lm:0:tv:SNP:het
        # 12      122278216       .       G       GAA,GAA .       .       BS=122278216    GT:QQ:BD:BK:BI:BVT:BLT  .:.:.:.:.:NOCALL:nocall 2/1:0:FP:am:i1_5:INDEL:hetalt
        # When considering genotypes, this should be one FN and two FP. When ignoring genotypes,
        # it should be only a single FP and one TP actually. We therefore also yield some variant identifying info, so that
        # code using the comparator can decide whether a certain FP is reported twice.
        #
        # Another example:
        # 10      24984382        .       G       C,T,GA  .       .       BS=24984382     GT:BD:BK:QQ:BI:BVT:BLT  3/3:TP:gm:0:i1_5:INDEL:homalt   1/2:FP:lm:0:tv:SNP:hetalt
        # 10      24984382        .       G       GA,GAA  .       .       BS=24984382     GT:QQ:BD:BK:BI:BVT:BLT  .:.:.:.:.:NOCALL:nocall 1/2:0:FP:lm:i1_5:INDEL:hetalt
        # 10      24984382        .       G       GA,GA   .       .       BS=24984382     GT:QQ:BD:BK:BI:BVT:BLT  .:.:.:.:.:NOCALL:nocall 2/1:0:TP:gm:i1_5:INDEL:hetalt
        # 10      24984382        .       G       T       .       .       BS=24984382     GT:QQ:BD:BK:BI:BVT:BLT  .:.:.:.:.:NOCALL:nocall 1/0:0:FP:lm:tv:SNP:het
        #
        # Hence, we need to be very careful below when relabeling calls and also ensure that duplicates are
        # removed before we yield anything.

        toyield = []

        def register(classification):
            if classification.variant not in self.visited[classification.cls]:
                toyield.append(classification)
                self.visited[classification.cls].add(classification.variant)

        truth_gt, query_gt = get_gt(record, 0), get_gt(record, 1)
        truth_alts, query_alts = alt_alleles(record, 0), alt_alleles(record, 1)

        for classification in super().classify(record):
            var = classification.variant

            def modified(new_cls, alt):
                return Classification(
                    new_cls,
                    Variant(
                        var.contig,
                        var.pos,
                        var.ref,
                        frozenset([alt]),
                    ),
                )

            if classification.cls == Class.FP:

                # determine whether this is only FP because the genotype does not match
                for alt in truth_alts & query_alts:
                    # alleles match, yield as TP
                    register(modified(Class.TP_query, alt))
                for alt in query_alts - truth_alts:
                    # not seen as TP before, yield as FP
                    new_classification = modified(Class.FP, alt)
                    if new_classification.variant not in self.visited[Class.TP_truth]:
                        register(new_classification)
            elif classification.cls == Class.FN:
                # determine whether this is only FN because the genotype does not match
                for alt in truth_alts & query_alts:
                    # alleles match, yield as TP
                    register(modified(Class.TP_truth, alt))
                for alt in truth_alts - query_alts:
                    register(modified(Class.FN, alt))
            else:
                for alt in classification.variant.alts:
                    register(modified(classification.cls, alt))
        yield from toyield
