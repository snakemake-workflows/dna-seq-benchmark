from collections import defaultdict
import sys, os

# sys.path.insert(0, os.path.dirname(__file__))
sys.stderr = open(snakemake.log[0], "w")

import csv
import pysam

from common.classification import CompareExistence, Class, is_het

cmp = CompareExistence()
varfile = pysam.VariantFile(snakemake.input.calls)

with open(snakemake.output[0], "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(
        [
            "class",
            "chromosome",
            "position",
            "ref_allele",
            "alt_allele",
            "true_genotype",
            "predicted_genotype",
        ]
    )
    for record in varfile:
        for c in cmp.classify(record):
            if c.cls == Class.FP and snakemake.wildcards.classification == "fp":
                classification = "FP"
                truth_gt = "0/0"
                query_gt = "0/1" if is_het(record, 1, c.variant) else "1/1"
            elif c.cls == Class.FN and snakemake.wildcards.classification == "fn":
                classification = "FN"
                truth_gt = "0/1" if is_het(record, 0, c.variant) else "1/1"
                query_gt = ""
            else:
                continue

            for alt in c.variant.alts:
                writer.writerow(
                    [
                        classification,
                        record.contig,
                        record.pos,
                        record.ref,
                        alt,
                        truth_gt,
                        query_gt,
                    ]
                )
