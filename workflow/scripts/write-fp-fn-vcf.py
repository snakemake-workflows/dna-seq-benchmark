import sys, os

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from pysam import VariantFile, VariantRecord

def load_fp_fn_table(file_path: str) -> pd.DataFrame:
    """Load TSV file into a pandas DataFrame."""
    return pd.read_csv(file_path, sep="\t")

def load_vcf(file_path: str) -> VariantFile:
    """Load a VCF file using pysam."""
    return VariantFile(file_path)

def get_variant_chr_pos(df: pd.DataFrame) -> pd.DataFrame:
    """Extract chromosome and position from the DataFrame."""
    if 'chromosome' not in df.columns or 'position' not in df.columns:
        raise ValueError("DataFrame must contain 'chromosome' and 'position' columns.")
    chr_pos = df[['chromosome', 'position']].drop_duplicates()
    return chr_pos

def get_variant_record(vcf: VariantFile, chr: str, pos: str) -> VariantRecord:
    """Get a specific variant record from the VCF file."""
    chr = str(chr)
    pos = int(pos)
    for record in vcf.fetch(chr, (pos-1), (pos+1)):
        if record.chrom == chr and record.pos == pos:
            return record
    raise ValueError(f"Variant not found in VCF: {chr}:{pos}")

def collect_records(vcf: VariantFile, df: pd.DataFrame) -> list:
    """Collect variant records from the VCF file based on the DataFrame."""
    records = []
    for _, row in df.iterrows():
        try:
            record = get_variant_record(vcf, row['chromosome'], row['position'])
            records.append(record)
        except ValueError as e:
            print(f"Error fetching record for {row['chromosome']}:{row['position']}: {e}", file=sys.stderr)
    print(f"Collected {len(records)} records.", file=sys.stderr)
    return records

def write_vcf(vcf_in: VariantFile, records: list, output_file: str):
    """Write collected records to a VCF file."""
    if not records:
        print("No records to write.", file=sys.stderr)
        return
    with VariantFile(output_file, 'w', header=vcf_in.header) as vcf_out:
        for record in records:
            vcf_out.write(record)
    print(f"VCF written to {output_file}", file=sys.stderr)


def process_variants(table_path: str, vcf_path: str, output_file: str):
    """Process variants from table and write to VCF."""
    table = load_fp_fn_table(table_path)
    vcf = load_vcf(vcf_path)

    if table.empty:
        print("No variants to process.", file=sys.stderr)
        # write an empty VCF file
        with VariantFile(output_file, 'w', header=vcf.header) as vcf_out:
            pass  # Just write the header
    else:
        variants = collect_records(vcf, get_variant_chr_pos(table))
        write_vcf(vcf, variants, output_file)

if snakemake.params.get("output") == "shared-fn" or snakemake.params.get("output") == "unique-fn":
    process_variants(
        snakemake.input.benchmark_table,
        snakemake.input.truth_vcf,
        snakemake.output[0]
    )

elif snakemake.params.get("output") == "unique-fp":
    process_variants(
        snakemake.input.benchmark_table,
        snakemake.input.callset_vcf,
        snakemake.output[0]
    )