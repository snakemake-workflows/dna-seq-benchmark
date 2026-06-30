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
        print(f"Warning: AD field missing for variant at {variant.CHROM}:{variant.POS}. Skipping VAF calculation.", file=sys.stderr)
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


def _num_samples_from_samples(samples):
    """Determine number of samples from VCF samples list or variant object."""
    if isinstance(samples, (list, tuple)):
        return len(samples)
    elif hasattr(samples, '__len__'):
        return len(samples)
    return 0


def calculate_vaf_from_fields(variant, field, name, samples=None):
    """
    Calculates the VAF for each sample by dividing two numeric fields.

    Args:
        variant (cyvcf2.Variant): The variant record.
        field (str): Either "INFO" or "FORMAT".
        name (str): The field ID (e.g. "NUMERATOR" / "DENOMINATOR").
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

    # Fetch numerator and denominator values
    if field == "INFO":
        num_val = variant.INFO.get(name)
        den_val = variant.INFO.get(f"{name}_den")
        if num_val is None or den_val is None:
            return None
    elif field == "FORMAT":
        try:
            num_val = variant.format(name)
            den_val = variant.format(f"{name}_den")
        except KeyError:
            return None
    else:
        raise ValueError(f"Unsupported field location: {field}. Must be 'INFO' or 'FORMAT'.")

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

    Args:
        input_vcf_path: Path to the input VCF.
        output_vcf_path: Path for the output VCF.
        calculate_from_ad: If True, calculate VAF from AD FORMAT field.
        num_field/den_field: Field location ("INFO" or "FORMAT") for custom calc.
        num_name/den_name: Field name for numerator and denominator.
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
            vaf_array = calculate_vaf_from_ad(variant, vcf_reader.samples)
        elif num_field and num_name and den_field and den_name:
            vaf_array = calculate_vaf_from_fields(variant, num_field, num_name)

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
