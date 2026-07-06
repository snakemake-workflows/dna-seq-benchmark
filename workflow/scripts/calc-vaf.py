from cyvcf2 import VCF, Writer
import sys
import argparse
import numpy as np


def calculate_vaf_from_ad(variant, samples):
    """
    Calculates the Variant Allele Frequency (VAF) for each sample and
    each alternate allele based on the AD (Allelic Depth) FORMAT field.

    Returns:
        numpy.ndarray: shape (n_samples, n_alt_alleles) with VAFs.
                       Returns None if AD field is missing.
                       Missing values or divisions by zero result in np.nan.
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
            alt_depth = sample_ad[j + 1]
            vaf_values[i, j] = alt_depth / total_depth

    return vaf_values


def calculate_vaf_from_strelka(variant, samples):
    """
    Calculate VAF from Strelka-specific FORMAT fields.

    For biallelic SNVs: uses FORMAT/{BASE}U fields (AU, TU, GU, CU)
        - ref_counts = FORMAT/{REF}U (e.g., FORMAT/AU if REF='A')
        - alt_counts = FORMAT/{ALT}U (e.g., FORMAT/TU if ALT='T')
        - VAF = alt_counts[tier0] / (alt_counts[tier0] + ref_counts[tier0])

    For biallelic indels: uses FORMAT/TAR and FORMAT/TIR
        - VAF = TAR[tier0] / (TAR[tier0] + TIR[tier0])

    Args:
        variant: cyvcf2.Variant record
        samples: list of sample names

    Returns:
        numpy.ndarray: shape (n_samples, n_alt_alleles), or None if fields are
                       missing or variant has >1 alternate allele.
    """
    n_samples = len(samples)
    n_alt_alleles = len(variant.ALT)

    if n_alt_alleles == 0:
        return np.full((n_samples, 0), np.nan, dtype=np.float32)

    vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)

    # Strelka per-nucleotide fields use Number=A, so we only handle
    # single-allele (biallelic) variants reliably.
    if n_alt_alleles > 1:
        return None

    ref_base = variant.REF
    alt_base = variant.ALT[0]

    # Biallelic SNV: REF and ALT are both single nucleotides, different bases
    if (
        len(ref_base) == 1
        and len(alt_base) == 1
        and ref_base in "ACGT"
        and alt_base in "ACGT"
        and ref_base != alt_base
    ):
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
            ref_count = float(ref_arr[i])
            alt_count = float(alt_arr[i])
            total = ref_count + alt_count
            if total > 0:
                vaf_values[i, 0] = alt_count / total

    # Biallelic indel: REF and ALT differ in length
    elif len(ref_base) != len(alt_base):
        try:
            tar_arr = variant.format("TAR")
            tir_arr = variant.format("TIR")
        except KeyError:
            return None

        if tar_arr is None or tir_arr is None:
            return None

        for i in range(n_samples):
            alt_count = float(tar_arr[i])
            ref_count = float(tir_arr[i])
            total = ref_count + alt_count
            if total > 0:
                vaf_values[i, 0] = alt_count / total

    else:
        return None

    return vaf_values


def _num_samples_from_samples(samples):
    """Determine number of samples from VCF samples list or variant object."""
    if isinstance(samples, (list, tuple)):
        return len(samples)
    elif hasattr(samples, '__len__'):
        return len(samples)
    return 0


def calculate_vaf_from_fields(variant, field, name, den_field=None, den_name=None, samples=None):
    """
    Calculates the VAF for each sample by dividing two numeric fields.

    Args:
        variant (cyvcf2.Variant): The variant record.
        field (str): Either "INFO" or "FORMAT".
        name (str): The field ID for the numerator.
        den_field (str or None): Either "INFO" or "FORMAT" for the denominator.
            Defaults to ``field`` when None.
        den_name (str): The field ID for the denominator. Must be provided by the caller.
        samples (list or None): List of sample names, or None to infer from variant.

    Returns:
        numpy.ndarray: shape (n_samples, n_alt_alleles), or None if fields are missing.
    """
    n_samples = _num_samples_from_samples(samples) if samples else 1

    if n_samples == 0:
        return np.array([], dtype=np.float32).reshape(0, 0)

    # Determine how many alternate alleles this variant has
    n_alt_alleles = len(variant.ALT)

    if n_alt_alleles == 0:
        return np.array([], dtype=np.float32).reshape(n_samples, 0)

    # Fetch numerator and denominator values from their respective field locations
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

    # Handle both scalar (INFO) and array (FORMAT) values
    if field == "INFO":
        # INFO fields: compute VAF per allele across all samples
        if isinstance(num_val, (list, np.ndarray)):
            # Per-allele values (Number=A or similar)
            if len(num_val) != n_alt_alleles:
                print(f"Warning: mismatch in {name} array length for variant "
                      f"{variant.CHROM}:{variant.POS}", file=sys.stderr)
                return None
            # All samples share the same numerator per-allele; we assume
            # total depth is the denominator for the first sample only
            # (or we use den as a scalar per allele).
            # General approach: if den is also per-allele, compute element-wise.
            if isinstance(den_val, (list, np.ndarray)) and len(den_val) == n_alt_alleles:
                vaf = np.array(den_val, dtype=np.float32)
                for i in range(n_alt_alleles):
                    if vaf[i] != 0:
                        vaf[i] = num_val[i] / vaf[i]
                    else:
                        vaf[i] = np.nan
                # Broadcast to all samples
                vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
                for i in range(n_samples):
                    vaf_values[i, :] = vaf
                return vaf_values
            else:
                # den is scalar — use it as total depth, numerator is per-allele count
                den_scalar = float(den_val) if isinstance(den_val, (int, float, np.integer, np.floating)) else 0
                vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
                for i in range(n_samples):
                    for j in range(n_alt_alleles):
                        if den_scalar > 0:
                            vaf_values[i, j] = float(num_val[j]) / den_scalar
                return vaf_values
        else:
            # Scalar numerator — broadcast across all alt alleles
            num_scalar = float(num_val)
            den_scalar = float(den_val) if isinstance(den_val, (int, float, np.integer, np.floating)) else 0
            vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
            for i in range(n_samples):
                for j in range(n_alt_alleles):
                    if den_scalar > 0:
                        vaf_values[i, j] = num_scalar / den_scalar
            return vaf_values

    else:  # FORMAT
        # FORMAT fields: shape is (n_samples, ... )
        num_arr = np.array(num_val, dtype=np.float32)
        den_arr = np.array(den_val, dtype=np.float32)

        # Handle different Number values
        if num_arr.shape == (n_samples,):
            # Number=1 — scalar per sample, assume scalar denominator
            if np.isscalar(den_val) or den_arr.shape == ():
                den_scalar = float(den_val) if np.isscalar(den_val) else float(den_arr)
            elif den_arr.shape == (n_samples,):
                den_scalar = den_arr
            else:
                den_scalar = float(den_val) if np.isscalar(den_val) else den_arr[0]

            vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
            for i in range(n_samples):
                d = den_scalar[i] if np.ndim(den_scalar) > 0 else den_scalar
                if isinstance(d, (int, float, np.integer, np.floating)) and d != 0:
                    for j in range(n_alt_alleles):
                        vaf_values[i, j] = float(num_arr[i]) / d
            return vaf_values

        elif num_arr.shape == (n_samples, n_alt_alleles):
            # Number=A — one value per sample per alt allele
            if den_arr.shape == (n_samples,):
                # Denominator is scalar per sample
                vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
                for i in range(n_samples):
                    d = float(den_arr[i])
                    if d != 0:
                        for j in range(n_alt_alleles):
                            vaf_values[i, j] = float(num_arr[i, j]) / d
                return vaf_values
            elif den_arr.shape == (n_samples, n_alt_alleles):
                # Both are per-sample per-allele
                vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
                for i in range(n_samples):
                    for j in range(n_alt_alleles):
                        n_val = float(num_arr[i, j])
                        d_val = float(den_arr[i, j])
                        if d_val != 0:
                            vaf_values[i, j] = n_val / d_val
                return vaf_values
            else:
                # Scalar denominator
                den_scalar = float(den_val) if np.isscalar(den_val) else 0
                vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)
                for i in range(n_samples):
                    for j in range(n_alt_alleles):
                        if den_scalar != 0:
                            vaf_values[i, j] = float(num_arr[i, j]) / den_scalar
                return vaf_values

        else:
            print(f"Warning: unexpected shape for FORMAT field {name}: "
                  f"{num_arr.shape} at variant {variant.CHROM}:{variant.POS}",
                  file=sys.stderr)
            return None


def add_vaf_to_vcf(input_vcf_path, output_vcf_path,
                   calculate_from_ad=True,
                   num_field=None, num_name=None,
                   den_field=None, den_name=None):
    """
    Reads an input VCF, calculates VAF, and writes to an output VCF file.

    When calculate_from_ad=True, the script first tries Strelka-specific
    FORMAT fields (AU/TU/GU/CU for SNVs, TAR/TIR for indels) and falls back
    to the standard AD FORMAT field if those are not present.
    """
    try:
        vcf_reader = VCF(input_vcf_path)
    except Exception as e:
        print(f"Error opening input VCF file '{input_vcf_path}': {e}", file=sys.stderr)
        sys.exit(1)

    # Add the new VAF FORMAT field definition to the header
    try:
        vcf_reader.add_format_to_header({
            'ID': 'VAF',
            'Description': 'Variant Allele Frequency',
            'Type': 'Float',
            'Number': 'A'
        })
    except ValueError:
        pass  # Field may already exist

    try:
        vcf_writer = Writer(output_vcf_path, vcf_reader)
    except Exception as e:
        print(f"Error creating output VCF file '{output_vcf_path}': {e}", file=sys.stderr)
        vcf_reader.close()
        sys.exit(1)

    print(f"Processing VCF: {input_vcf_path}", file=sys.stderr)
    print(f"Writing output to: {output_vcf_path}", file=sys.stderr)

    processed_count = 0
    for variant in vcf_reader:
        vaf_array = None

        if calculate_from_ad:
            # Try Strelka-specific fields first, then fall back to AD
            vaf_array = calculate_vaf_from_strelka(variant, vcf_reader.samples)
            if vaf_array is None:
                vaf_array = calculate_vaf_from_ad(variant, vcf_reader.samples)
        elif num_field and num_name and den_field and den_name:
            vaf_array = calculate_vaf_from_fields(variant, num_field, num_name, den_field=den_field, den_name=den_name, samples=vcf_reader.samples)

        if vaf_array is not None:
            try:
                variant.set_format('VAF', vaf_array)
            except Exception as e:
                print(f"Error setting VAF format for variant at "
                      f"{variant.CHROM}:{variant.POS}: {e}", file=sys.stderr)

        vcf_writer.write_record(variant)
        processed_count += 1
        if processed_count % 1000 == 0:
            print(f"Processed {processed_count} variants...", file=sys.stderr)

    vcf_reader.close()
    vcf_writer.close()
    print(f"Finished processing. Total variants processed: {processed_count}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Calculate VAF and add to VCF")
    parser.add_argument("input_vcf", help="Input VCF/BCF file")
    parser.add_argument("output_vcf", help="Output VCF/BCF file")
    parser.add_argument("--from-ad", action="store_true", default=False,
                        help="Calculate VAF from AD FORMAT field (default)")
    parser.add_argument("--num-field", type=str, default=None,
                        help="Field location for numerator: 'INFO' or 'FORMAT'")
    parser.add_argument("--num-name", type=str, default=None,
                        help="Field name for numerator")
    parser.add_argument("--den-field", type=str, default=None,
                        help="Field location for denominator: 'INFO' or 'FORMAT'")
    parser.add_argument("--den-name", type=str, default=None,
                        help="Field name for denominator (appends '_den' to num_name if not specified)")

    args = parser.parse_args()

    # If den_name not specified, derive from num_name
    if args.den_name is None and args.num_name is not None:
        args.den_name = f"{args.num_name}_den"

    calc_from_ad = args.from_ad
    calc_from_fields = (args.num_field is not None and args.num_name is not None
                        and args.den_field is not None and args.den_name is not None)

    if not calc_from_ad and not calc_from_fields:
        print("Error: specify --from-ad or provide --num-field/--num-name/--den-field/--den-name",
              file=sys.stderr)
        sys.exit(1)

    add_vaf_to_vcf(
        input_vcf_path=args.input_vcf,
        output_vcf_path=args.output_vcf,
        calculate_from_ad=calc_from_ad,
        num_field=args.num_field,
        num_name=args.num_name,
        den_field=args.den_field,
        den_name=args.den_name,
    )


if __name__ == "__main__":
    main()
