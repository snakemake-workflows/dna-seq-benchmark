import sys, os
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def load_variant_table(file_path: str, cls) -> pd.DataFrame:
    """Load TSV file into a pandas DataFrame."""
    variant_cols = ["chromosome", "position", "ref_allele", "alt_allele"]
    df = pd.read_csv(file_path, sep="\t")

    callset_variant_totals = df.groupby("callset").size().to_dict()
    total_callsets = df["callset"].nunique()
    print(f"total number {snakemake.wildcards.benchmark} callsets containing {cls}:", total_callsets, file=sys.stderr)

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


def write_per_callset_variants(df: pd.DataFrame, output_dir: str): #, benchmark: str, cls: str, callset_totals: dict):
    """Write one TSV per callset for variants that occur only in that callset and collect summary info."""
    os.makedirs(output_dir, exist_ok=True)  # create output folder if missing
    for callset, group_df in df.groupby("callset"):
        safe_name = callset.replace(" ", "_").replace("/", "_")
        filename = os.path.join(output_dir, f"unique_to_{safe_name}.tsv")
        if group_df.empty:
            print(f"Warning: DataFrame is empty. Not writing to {filename}", file=sys.stderr)
            continue
        write_output(group_df, filename)


def annotate_variant_callset_counts(df: pd.DataFrame, variant_cols: list) -> pd.DataFrame:
    """Annotate each variant with the number of unique callsets it appears in."""
    counts = df.groupby(variant_cols)["callset"].nunique().reset_index()
    counts.rename(columns={"callset": "callset_count"}, inplace=True)
    return df.merge(counts, on=variant_cols)


def filter_variants(df: pd.DataFrame, callset_count: int = None, total_callsets: int = None) -> pd.DataFrame:
    """Filter variants based on number or set of callsets."""
    variant_cols = ["chromosome", "position", "ref_allele", "alt_allele"]
    if callset_count is not None:
        df = df[df["callset_count"] == callset_count]
    elif total_callsets is not None:
        df = df[df["callset_count"] == total_callsets]

    return df


def filter_variants_by_exact_sets(df: pd.DataFrame, exact_sets: dict, variant_cols: list) -> dict:
    """Return a dict mapping labels to filtered DataFrames matching each exact callset set."""
    # First: compute set of callsets for each variant
    variant_callsets = df.groupby(variant_cols)["callset"].unique().reset_index()
    variant_callsets["callset_set"] = variant_callsets["callset"].apply(lambda x: frozenset(x))

    label_to_df = {}

    for callset_set, label in exact_sets.items():
        matches = variant_callsets[variant_callsets["callset_set"] == callset_set]
        matched_variants = df.merge(matches[variant_cols], on=variant_cols)
        label_to_df[label] = matched_variants

    return label_to_df


def compare_benchmarks(df: pd.DataFrame, benchmark1: str, benchmark2: str) -> pd.DataFrame:
    """Compare two callsets and return a DataFrame with variants unique to each callset."""
    # Filter for the two callsets
    df1 = df[df["callset"].str.contains(benchmark1)]
    df2 = df[df["callset"].str.contains(benchmark2)]

    # Find unique variants in each callset
    unique_to_benchmark1 = df1[~df1.set_index(["chromosome", "position", "ref_allele", "alt_allele"]).index.isin(df2.set_index(["chromosome", "position", "ref_allele", "alt_allele"]).index)]
    unique_to_benchmark2 = df2[~df2.set_index(["chromosome", "position", "ref_allele", "alt_allele"]).index.isin(df1.set_index(["chromosome", "position", "ref_allele", "alt_allele"]).index)]

    # Combine results
    result = pd.concat([unique_to_benchmark1, unique_to_benchmark2], ignore_index=True)
    return result


def collect_summary(df: pd.DataFrame, benchmark: str, cls: str, filter_type: str, callset_info, callset_totals: dict = None) -> pd.DataFrame:
    """Collect summary with variant count per callset or callset set, including total variant counts when applicable."""
    if isinstance(callset_info, str):
        total_variants = callset_totals.get(callset_info, None) if callset_totals else None
        # callset-specific (e.g. unique_to_one or in_all)
        result = pd.DataFrame([{
            "benchmark": benchmark,
            "classification": cls,
            "filter_type": filter_type,
            "callset": callset_info,
            "variant_count": len(df),
            "total_variant_count": total_variants
        }])
    else:
        # # callset set (e.g. for exact sets)
        result = pd.DataFrame([{
            "benchmark": benchmark,
            "classification": cls,
            "filter_type": filter_type,
            "callset": ",".join(sorted(callset_info)),
            "variant_count": len(df),
            "total_variant_count": None
        }])
    return result



# fp variants present only in one callset
df_fp, callset_variant_totals_fp, total_callsets_fp = load_variant_table(snakemake.input.fp, "fp")
one_callset_df_fp = filter_variants(df_fp, callset_count=1)
output_dir = snakemake.output.unique_fp
write_per_callset_variants(one_callset_df_fp, output_dir) #, snakemake.wildcards.benchmark, "fp", df_fp, callset_variant_totals_fp)



df_fn, callset_variant_totals_fn, total_callsets_fn = load_variant_table(snakemake.input.fn, "fn")
# Variants present in all callsets
all_callsets_df_fn = filter_variants(df_fn, callset_count=total_callsets_fn)
write_output(all_callsets_df_fn, snakemake.output.shared_fn)

# fn variants present only in one callset
one_callset_df_fn = filter_variants(df_fn, callset_count=1)
output_dir = snakemake.output.unique_fn
write_per_callset_variants(one_callset_df_fn, output_dir)#, snakemake.wildcards.benchmark, "fn", df_fn, callset_variant_totals_fn)