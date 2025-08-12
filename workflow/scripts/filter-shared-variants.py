import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def load_variant_table(file_path: str, cls) -> tuple[pd.DataFrame, dict, int]:
     """Load TSV file into a pandas DataFrame."""
     variant_cols = ["chromosome", "position", "ref_allele", "alt_allele"]
     df = pd.read_csv(file_path, sep="\t")

     if df.empty:
         callset_variant_totals = 0
         total_callsets = 0
     else:
         callset_variant_totals = df.groupby("callset").size().to_dict()
         total_callsets = df["callset"].nunique()
         print(f"total number of {snakemake.wildcards.benchmark} callsets containing {cls}: {total_callsets}", file=sys.stderr)

         # Annotate with number of callsets per variant
         df = annotate_variant_callset_counts(df, variant_cols)

     return (df, callset_variant_totals, total_callsets)

def write_output(df: pd.DataFrame, filename: str):
    """Write DataFrame to given TSV file."""
    if df.empty:
        print(f"Warning: DataFrame is empty. Writing empty file to {filename}", file=sys.stderr)
    # Write DataFrame to TSV file
    df.to_csv(filename, sep="\t", index=False)
    print(f"Written: {filename}", file=sys.stderr)


def annotate_variant_callset_counts(df: pd.DataFrame, variant_cols: list) -> pd.DataFrame:
    """Annotate each variant with the number of unique callsets it appears in."""
    counts = df.groupby(variant_cols)["callset"].nunique().reset_index()
    counts.rename(columns={"callset": "callset_count"}, inplace=True)
    return df.merge(counts, on=variant_cols)


def filter_variants(df: pd.DataFrame, callset_count: int = None, total_callsets: int = None) -> pd.DataFrame:
    """Filter variants based on number or set of callsets."""
    # Return empty dataframe if input is empty or callset_count column doesn't exist
    if df.empty or "callset_count" not in df.columns:
        return df
        
    if callset_count is not None:
        df = df[df["callset_count"] == callset_count]
    elif total_callsets is not None:
        df = df[df["callset_count"] == total_callsets]

    return df


# Load fn variants
df_fn, callset_variant_totals_fn, total_callsets_fn = load_variant_table(snakemake.input.fn, "fn")

# variants present in all callsets
if df_fn.empty:
    write_output(df_fn, snakemake.output.shared_fn)
else:
    all_callsets_df_fn = filter_variants(df_fn, callset_count=total_callsets_fn)
    write_output(all_callsets_df_fn, snakemake.output.shared_fn)
