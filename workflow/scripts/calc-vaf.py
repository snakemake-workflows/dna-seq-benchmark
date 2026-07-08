import sys
import shlex

sys.stderr = open(snakemake.log[0], "w")

from cyvcf2 import VCF, Writer
import numpy as np


def calculate_vaf_from_ad(variant, samples):
    """
    Calculates the Variant Allele Frequency (VAF) for each sample and
    each alternate allele based on the AD (Allelic Depth) FORMAT field.
    """
    try:
        ad = variant.format('AD')
    except KeyError:
        return None
    if ad is None:
        return None

    n_samples = len(samples)
    n_alt_alleles = len(variant.ALT)
    vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)

    for i in range(n_samples):
        sample_ad = ad[i]
        if np.any(sample_ad < 0):
            continue
        total_depth = np.sum(sample_ad)
        if total_depth == 0:
            continue
        for j in range(n_alt_alleles):
            vaf_values[i, j] = sample_ad[j + 1] / total_depth

    return vaf_values


def calculate_vaf_from_strelka(variant, samples):
    """
    Calculate VAF from Strelka-specific FORMAT fields.
    """
    n_samples = len(samples)
    n_alt_alleles = len(variant.ALT)
    if n_alt_alleles == 0:
        return np.full((n_samples, 0), np.nan, dtype=np.float32)
    if n_alt_alleles > 1:
        return None

    vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
    ref_base = variant.REF
    alt_base = variant.ALT[0]

    if len(ref_base) == 1 and len(alt_base) == 1 and ref_base != alt_base and ref_base in "ACGT" and alt_base in "ACGT":
        ref_field = f"{ref_base}U"
        alt_field = f"{alt_base}U"
        try:
            ref_arr = variant.format(ref_field)
            alt_arr = variant.format(alt_field)
        except KeyError:
            return None
        if ref_arr is None or alt_arr is None:
            return None
        for i in range(n_samples):
            total = float(ref_arr[i, 0]) + float(alt_arr[i, 0])
            if total > 0:
                vaf_values[i, 0] = float(alt_arr[i, 0]) / total

    elif len(ref_base) != len(alt_base):
        try:
            tar_arr = variant.format("TAR")
            tir_arr = variant.format("TIR")
        except KeyError:
            return None
        if tar_arr is None or tir_arr is None:
            return None
        for i in range(n_samples):
            total = float(tar_arr[i, 0]) + float(tir_arr[i, 0])
            if total > 0:
                vaf_values[i, 0] = float(tir_arr[i, 0]) / total

    else:
        return None

    return vaf_values


def calculate_vaf_from_fields(variant, samples, field, name, den_field=None, den_name=None):
    """
    Calculates the VAF for each sample by dividing two numeric fields.
    """
    n_samples = len(samples)
    n_alt_alleles = len(variant.ALT)
    if n_alt_alleles == 0:
        return np.array([], dtype=np.float32).reshape(n_samples, 0)

    den_key = den_name
    dfield = den_field if den_field is not None else field

    if field == "INFO":
        num_val = variant.INFO.get(name)
    elif field == "FORMAT":
        try:
            num_val = variant.format(name)
        except KeyError:
            return None
    else:
        raise ValueError(f"Unsupported numerator field location: {field}. Must be 'INFO' or 'FORMAT'.")
    if num_val is None:
        return None

    if dfield == "INFO":
        den_val = variant.INFO.get(den_key)
    elif dfield == "FORMAT":
        try:
            den_val = variant.format(den_key)
        except KeyError:
            return None
    else:
        raise ValueError(f"Unsupported denominator field location: {dfield}. Must be 'INFO' or 'FORMAT'.")
    if den_val is None:
        return None

    vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)

    if field == "INFO" and isinstance(num_val, (list, np.ndarray)):
        if len(num_val) != n_alt_alleles:
            print(f"Warning: mismatch in {name} array length for variant "
                  f"{variant.CHROM}:{variant.POS}", file=sys.stderr)
            return None
        if isinstance(den_val, (list, np.ndarray)) and len(den_val) == n_alt_alleles:
            vaf = np.array(den_val, dtype=np.float32)
            for i in range(n_alt_alleles):
                if vaf[i] != 0:
                    vaf[i] = num_val[i] / vaf[i]
                else:
                    vaf[i] = np.nan
            for i in range(n_samples):
                vaf_values[i, :] = vaf
        else:
            den_scalar = float(den_val) if isinstance(den_val, (int, float, np.integer, np.floating)) else 0
            for i in range(n_samples):
                for j in range(n_alt_alleles):
                    if den_scalar > 0:
                        vaf_values[i, j] = float(num_val[j]) / den_scalar
    elif field == "INFO":
        num_scalar = float(num_val)
        den_scalar = float(den_val) if isinstance(den_val, (int, float, np.integer, np.floating)) else 0
        for i in range(n_samples):
            for j in range(n_alt_alleles):
                if den_scalar > 0:
                    vaf_values[i, j] = num_scalar / den_scalar
    else:  # FORMAT
        num_arr = np.array(num_val, dtype=np.float32)
        den_arr = np.array(den_val, dtype=np.float32)

        if num_arr.shape == (n_samples,):
            den_arr = np.atleast_1d(np.squeeze(den_arr))
            if len(den_arr) == 1:
                den_scalar = float(den_arr[0])
            else:
                den_scalar = den_arr
            for i in range(n_samples):
                d = float(den_scalar[i]) if np.ndim(den_scalar) > 0 else float(den_scalar)
                if d != 0:
                    for j in range(n_alt_alleles):
                        vaf_values[i, j] = float(num_arr[i]) / d
        elif num_arr.shape == (n_samples, n_alt_alleles):
            if den_arr.shape == (n_samples,):
                for i in range(n_samples):
                    d = float(den_arr[i])
                    if d != 0:
                        for j in range(n_alt_alleles):
                            vaf_values[i, j] = float(num_arr[i, j]) / d
            elif den_arr.shape == (n_samples, n_alt_alleles):
                for i in range(n_samples):
                    for j in range(n_alt_alleles):
                        if den_arr[i, j] != 0:
                            vaf_values[i, j] = float(num_arr[i, j]) / float(den_arr[i, j])
            else:
                den_arr = np.atleast_1d(np.squeeze(den_arr))
                if len(den_arr) == n_samples:
                    for i in range(n_samples):
                        d = float(den_arr[i])
                        if d != 0:
                            for j in range(n_alt_alleles):
                                vaf_values[i, j] = float(num_arr[i, j]) / d
                else:
                    return None

    return vaf_values


# Parse vaf_args string from snakemake.params
def parse_vaf_args(vaf_args_str):
    """Parse vaf_args string into structured params for the script directive."""
    parts = shlex.split(vaf_args_str)
    result = {"from_ad": False, "num_field": None, "num_name": None,
              "den_field": None, "den_name": None}
    i = 0
    while i < len(parts):
        if parts[i] == "--from-ad":
            result["from_ad"] = True
        elif parts[i] == "--num-field" and i + 1 < len(parts):
            result["num_field"] = parts[i + 1]
        elif parts[i] == "--num-name" and i + 1 < len(parts):
            result["num_name"] = parts[i + 1]
        elif parts[i] == "--den-field" and i + 1 < len(parts):
            result["den_field"] = parts[i + 1]
        elif parts[i] == "--den-name" and i + 1 < len(parts):
            result["den_name"] = parts[i + 1]
        i += 1
    if result["den_name"] is None and result["num_name"] is not None:
        result["den_name"] = f"{result['num_name']}_den"
    return result


# Main top-level code
vaf_args = parse_vaf_args(snakemake.params.get("vaf_args", ""))
calc_from_ad = vaf_args["from_ad"]
calc_from_fields = (vaf_args["num_field"] is not None
                    and vaf_args["num_name"] is not None
                    and vaf_args["den_field"] is not None
                    and vaf_args["den_name"] is not None)

if not calc_from_ad and not calc_from_fields:
    print("Error: specify --from-ad or provide --num-field/--num-name/--den-field/--den-name",
          file=sys.stderr)
    sys.exit(1)

try:
    vcf_reader = VCF(snakemake.input.bcf)
except Exception as e:
    print(f"Error opening input VCF file '{snakemake.input.bcf}': {e}", file=sys.stderr)
    sys.exit(1)

vcf_reader.add_format_to_header({
    'ID': 'VAF',
    'Description': 'Variant Allele Frequency',
    'Type': 'Float',
    'Number': 'A'
})

vcf_writer = Writer(snakemake.output.vcf, vcf_reader)
processed_count = 0

for variant in vcf_reader:
    vaf_array = None

    if calc_from_ad:
        vaf_array = calculate_vaf_from_strelka(variant, vcf_reader.samples)
        if vaf_array is None:
            vaf_array = calculate_vaf_from_ad(variant, vcf_reader.samples)
    elif calc_from_fields:
        vaf_array = calculate_vaf_from_fields(
            variant, vcf_reader.samples,
            vaf_args["num_field"], vaf_args["num_name"],
            den_field=vaf_args["den_field"], den_name=vaf_args["den_name"]
        )

    if vaf_array is not None:
        try:
            variant.set_format('VAF', vaf_array)
        except Exception as e:
            print(f"Error setting VAF format for variant at "
                  f"{variant.CHROM}:{variant.POS}: {e}", file=sys.stderr)

    vcf_writer.write_record(variant)
    processed_count += 1

vcf_reader.close()
vcf_writer.close()
