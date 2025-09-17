from collections import defaultdict
import sys, os

# sys.path.insert(0, os.path.dirname(__file__))
sys.stderr = open(snakemake.log[0], "w")

import csv
import pysam

from common.classification import CompareExistence, Class, is_het

cmp = CompareExistence()
varfile = pysam.VariantFile(snakemake.input.calls)
truth = pysam.VariantFile(snakemake.input.truth)
query = pysam.VariantFile(snakemake.input.query)
vaf_fields = snakemake.params.vaf_fields
if vaf_fields[0] is not None:
    vaf_field_query = vaf_fields[0]["field"]
    vaf_field_name_query = vaf_fields[0]["name"]
if vaf_fields[1] is not None:
    vaf_field_truth = vaf_fields[1]["field"]
    vaf_field_name_truth = vaf_fields[1]["name"]

with open(snakemake.output[0], "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(
        [
            "vaf",
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
                if vaf_fields[0] is not None:
                    r = list(query.fetch(record.contig, record.start, record.stop))
                    if len(r) > 0:
                        r = r[0]
                    vaf = r.info[vaf_field_name_query] if vaf_field_query == "INFO" else r.samples[0][vaf_field_name_query]
                else:
                    #No VAF information available -> float na
                    vaf = float('nan')

            elif c.cls == Class.FN and snakemake.wildcards.classification == "fn":
                classification = "FN"
                truth_gt = "0/1" if is_het(record, 0, c.variant) else "1/1"
                query_gt = ""
                if vaf_fields[1] is not None:
                    r = list(truth.fetch(record.contig, record.start, record.stop))
                    if len(r) > 0:
                        r = r[0]
                    vaf = r.info[vaf_field_name_truth] if vaf_field_truth == "INFO" else r.samples[0][vaf_field_name_truth]
                else:
                    #No VAF information available -> float na
                    vaf = float('nan')
            else:
                continue

            if isinstance(vaf, tuple):
                vaf = vaf[0]
            if isinstance(vaf, str):
                vaf = float(vaf.replace("%", "")) / 100
            for alt in c.variant.alts:
                writer.writerow(
                    [
                        vaf,
                        classification,
                        record.contig,
                        record.pos,
                        record.ref,
                        alt,
                        truth_gt,
                        query_gt,
                    ]
                )
