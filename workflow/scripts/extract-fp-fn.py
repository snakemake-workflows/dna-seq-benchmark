from collections import defaultdict
import sys
import numpy as np

sys.stderr = open(snakemake.log[0], "w")

import csv
import pysam

from common.classification import CompareExistence, Class, is_het


def get_vaf_from_record(record, field, name):
    """Extract VAF value from a record, handling FORMAT and INFO fields."""
    try:
        if field == "INFO":
            vaf = record.info.get(name)
        else:
            sample_name = list(record.samples.keys())[0]
            vaf = record.samples[sample_name].get(name)
    except (KeyError, IndexError, AttributeError):
        return float('nan')

    if isinstance(vaf, (list, tuple, np.ndarray, np.generic)):
        if hasattr(vaf, 'item'):
            vaf = vaf.item()
        elif len(vaf) > 0:
            vaf = vaf[0]
        else:
            return float('nan')

    if isinstance(vaf, str):
        vaf = float(vaf.replace("%", "")) / 100

    if hasattr(vaf, 'item'):
        vaf = vaf.item()

    try:
        return float(vaf)
    except (ValueError, TypeError):
        return float('nan')

cmp = CompareExistence()
varfile = pysam.VariantFile(snakemake.input.calls)
truth = pysam.VariantFile(snakemake.input.truth)
query = pysam.VariantFile(snakemake.input.query)
vaf_fields = snakemake.params.vaf_fields
if vaf_fields and len(vaf_fields) > 0 and vaf_fields[0] is not None:
    vaf_field_query = vaf_fields[0]["field"]
    vaf_field_name_query = vaf_fields[0]["name"]
else:
    vaf_field_query = None
    vaf_field_name_query = None
if vaf_fields and len(vaf_fields) > 1 and vaf_fields[1] is not None:
    vaf_field_truth = vaf_fields[1]["field"]
    vaf_field_name_truth = vaf_fields[1]["name"]
else:
    vaf_field_truth = None
    vaf_field_name_truth = None


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
                if vaf_field_query is not None:
                    r = list(query.fetch(record.contig, record.start, record.stop))
                    if len(r) > 0:
                        r = r[0]
                    vaf = get_vaf_from_record(r, vaf_field_query, vaf_field_name_query)
                else:
                    # No VAF information available -> float na
                    vaf = float('nan')

            elif c.cls == Class.FN and snakemake.wildcards.classification == "fn":
                classification = "FN"
                truth_gt = "0/1" if is_het(record, 0, c.variant) else "1/1"
                query_gt = ""
                if vaf_field_truth is not None:
                    r = list(truth.fetch(record.contig, record.start, record.stop))
                    if len(r) > 0:
                        r = r[0]
                    vaf = get_vaf_from_record(r, vaf_field_truth, vaf_field_name_truth)
                else:
                    # No VAF information available -> float na
                    vaf = float('nan')
            else:
                continue

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
