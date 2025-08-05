from cyvcf2 import VCF, Writer
import sys
import argparse
import numpy as np
                                  
def get_snv_allele_freq(variant):
    refCounts = variant.format(variant.REF + "U")
    altCounts = variant.format(variant.ALT[0] + "U")

    # TODO: check which value is the correct one from the matrix (this leads to many zero VAF)
    tier1RefCounts = refCounts[0, 0]
    tier1AltCounts = altCounts[0, 0]

    vaf = tier1AltCounts / (tier1AltCounts + tier1RefCounts)

    return vaf

def get_indel_allele_freq(variant):
    tier1RefCounts = variant.format("TAR")[0,0]
    tier1AltCounts = variant.format("TIR")[0,0]

    vaf = tier1AltCounts / (tier1AltCounts + tier1RefCounts)
    return vaf

def calculate_vaf(variant, samples):
    """
    Calculates the Variant Allele Frequency (VAF) for each sample and
    each alternate allele based on the AD (Allelic Depth) FORMAT field.

    Args:
        variant (cyvcf2.Variant): The variant object.
        samples (list): List of sample names.

    Returns:
        numpy.ndarray: A numpy array of shape (n_samples, n_alt_alleles)
                       containing the VAFs. Returns None if AD field is missing.
                       Missing values or divisions by zero result in np.nan.
    """
    try:
        # Get Allelic Depths (AD) - shape: (n_samples, n_alleles)
        # n_alleles includes the reference allele
        ad = variant.format('AD')
    except KeyError:
        # AD field is not present for this variant
        print(f"Warning: AD field missing for variant at {variant.CHROM}:{variant.POS}. Skipping VAF calculation.", file=sys.stderr)
        return None

    n_samples = len(samples)
    n_alleles = len(variant.alleles) # Includes reference
    n_alt_alleles = n_alleles - 1

    # Initialize VAF array with NaNs
    # Shape: (n_samples, n_alt_alleles)
    vaf_values = np.full((n_samples, n_alt_alleles), np.nan, dtype=np.float32)

    for i in range(n_samples):
        sample_ad = ad[i]

        # Check for missing AD data for the sample (represented by negative values or could be None depending on VCF)
        # cyvcf2 often uses negative numbers for missing integers in FORMAT fields like AD
        if np.any(sample_ad < 0):
             # Keep VAF as NaN if AD is missing for this sample
            continue

        # Calculate total depth for this sample
        total_depth = np.sum(sample_ad)

        if total_depth == 0:
            # Avoid division by zero, keep VAFs as NaN
            continue

        # Calculate VAF for each alternate allele
        for j in range(n_alt_alleles):
            alt_depth = sample_ad[j + 1] # AD format is [ref_depth, alt1_depth, alt2_depth, ...]
            vaf = alt_depth / total_depth
            vaf_values[i, j] = vaf

    return vaf_values

def add_vaf_to_vcf(input_vcf_path, output_vcf_path):
    """
    Reads an input VCF, calculates VAF for each variant/sample, adds it
    as a new FORMAT field 'VAF', and writes to an output VCF file.
    """
    # Open the input VCF file
    try:
        vcf_reader = VCF(input_vcf_path)
    except Exception as e:
        print(f"Error opening input VCF file '{input_vcf_path}': {e}", file=sys.stderr)
        sys.exit(1)


    # Add the new VAF FORMAT field definition to the header
    # Number='A' means one value per alternate allele
    # Type='Float' for the frequency value
    try:
        vcf_reader.add_format_to_header({
            'ID': 'VAF',
            'Description': 'Variant Allele Frequency calculated from AD field (Alt Depth / Total Depth)',
            'Type': 'Float',
            'Number': 'A' # One value per alternate allele
        })
    except ValueError as e:
         # Catch error if the field already exists
         print(f"Warning: FORMAT field 'VAF' might already exist in header: {e}", file=sys.stderr)


    # Create a VCF writer object using the modified header
    try:
        vcf_writer = Writer(output_vcf_path, vcf_reader)
    except Exception as e:
        print(f"Error creating output VCF file '{output_vcf_path}': {e}", file=sys.stderr)
        vcf_reader.close()
        sys.exit(1)

    print(f"Processing VCF: {input_vcf_path}")
    print(f"Writing output to: {output_vcf_path}")

    processed_count = 0
    # Iterate through each variant in the VCF
    for variant in vcf_reader:
        # Calculate VAFs for all samples for the current variant
        vaf_array = calculate_vaf(variant, vcf_reader.samples)

        # Add the calculated VAFs to the variant's FORMAT fields
        # The vaf_array must have shape (n_samples, n_alt_alleles)
        if vaf_array is not None:
            try:
                # Use set_format to add/update the VAF field for all samples
                variant.set_format('VAF', vaf_array)
            except Exception as e:
                 print(f"Error setting VAF format for variant at {variant.CHROM}:{variant.POS}: {e}", file=sys.stderr)
                 # Decide if you want to skip writing this variant or write without VAF
                 # Here, we'll still write the variant but VAF might be missing/incorrect

        # Write the (potentially modified) variant record to the output file
        vcf_writer.write_record(variant)
        processed_count += 1
        if processed_count % 1000 == 0:
            print(f"Processed {processed_count} variants...", file=sys.stderr)

    # Close the VCF reader and writer
    vcf_reader.close()
    vcf_writer.close()
    print(f"Finished processing. Total variants processed: {processed_count}")