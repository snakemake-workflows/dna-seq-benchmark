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


def _get_vaf(record, varfile, vaf_field, vaf_field_name):
    """Fetch VAF from a variant file."""
    if vaf_field is None:
        return float('nan')
    try:
        r = list(varfile.fetch(record.contig, record.start, record.stop))
        if len(r) > 0:
            r = r[0]
        else:
            return float('nan')
        vaf = r.info[vaf_field_name] if vaf_field == "INFO" else r.samples[0][vaf_field_name]
        return vaf
    except (KeyError, IndexError):
        return float('nan')


def _gt_str(gt):
    """Format GT tuple/list into a string, safe for partial no-calls."""
    if any(g is None for g in gt):
        return "/".join("." if g is None else str(g) for g in gt)
    return "/".join(str(g) for g in sorted(gt))


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
                    try:
                        vaf = r.info[vaf_field_name_query] if vaf_field_query == "INFO" else r.samples[0][vaf_field_name_query]
                    except (KeyError, IndexError):
                        vaf = float('nan')
                else:
                    #No VAF information available -> float na
                    vaf = float('nan')

            elif c.cls == Class.TP_query and snakemake.wildcards.classification == "tp":
                classification = "TP"
                truth_gt = _gt_str(record.samples[0].get("GT", [None, None]))
                query_gt = _gt_str(record.samples[1].get("GT", [None, None]))
                # For TP, try VAF from query first, fall back to truth
                vaf = _get_vaf(record, query, vaf_field_query, vaf_field_name_query)
                if vaf != vaf and vaf_field_truth is not None:
                    vaf = _get_vaf(record, truth, vaf_field_truth, vaf_field_name_truth)

            elif c.cls == Class.TP_truth and snakemake.wildcards.classification == "tp-baseline":
                classification = "TP-BASELINE"
                truth_gt = _gt_str(record.samples[0].get("GT", [None, None]))
                query_gt = _gt_str(record.samples[1].get("GT", [None, None]))
                # For TP-baseline, use truth VAF
                vaf = _get_vaf(record, truth, vaf_field_truth, vaf_field_name_truth)

            elif c.cls == Class.FN and snakemake.wildcards.classification == "fn":
                classification = "FN"
                truth_gt = "0/1" if is_het(record, 0, c.variant) else "1/1"
                query_gt = ""
                if vaf_fields[1] is not None:
                    r = list(truth.fetch(record.contig, record.start, record.stop))
                    if len(r) > 0:
                        r = r[0]
                    try:
                        vaf = r.info[vaf_field_name_truth] if vaf_field_truth == "INFO" else r.samples[0][vaf_field_name_truth]
                    except (KeyError, IndexError):
                        vaf = float('nan')
                else:
                    #No VAF information available -> float na
                    vaf = float('nan')
            else:
                continue

            # For TP and TP-BASELINE, VAF can come from either sample
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
