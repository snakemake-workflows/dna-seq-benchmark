# Configuration Guide

Please follow the instructions in the template `config.yaml`.

## VAF (Variant Allele Frequency) Calculation

If your variant-calling pipeline does not write a `VAF` FORMAT field to the VCF, you can have the workflow calculate it and inject it into the VCF automatically. This enables VAF-stratified precision/recall reporting.

### Option A: Calculate VAF from the AD (Allelic Depth) field

If your VCF contains a standard `AD` FORMAT field (i.e., `[ref_count, alt_count]` for biallelic sites), simply set:

```yaml
variant-calls:
  callset:
    my_callset:
      path: "path/to/variant-calls.vcf.gz"
      vaf-field: "tbc"          # ← "to-be-calculated"
      benchmark: giab-na12878
      genome-build: grch38
```

`vaf-field: "tbc"` tells the Snakemake workflow to run `calc-vaf.py --from-ad` on the VCF, computing VAF as:

```
VAF = alt_count / (ref_count + alt_count)
```

The calculated `VAF` FORMAT field (Type=Float, Number=A) is written to `results/calculate-vaf/{callset}.added-vaf.bcf`, which then feeds downstream precision/recall and FP/FN analyses.

### Option B: Calculate VAF from custom numerator / denominator fields

If your VCF lacks an `AD` field but has separate numerator and denominator FORMAT fields (e.g., `TUMOR_AF` and `TUMOR_DP`), specify both:

```yaml
variant-calls:
  callset:
    my_callset:
      path: "path/to/variant-calls.vcf.gz"
      vaf-field: "tbc"
      vaf-numerator:
        field: FORMAT
        name: TUMOR_AF
      vaf-denominator:
        field: FORMAT
        name: TUMOR_DP
      benchmark: giab-na12878
      genome-build: grch38
```

The workflow will run:

```
calc-vaf.py --num-field FORMAT --num-name TUMOR_AF \
            --den-field FORMAT --den-name TUMOR_DP
```

which computes `VAF = TUMOR_AF / TUMOR_DP` per sample.

The same `vaf-numerator`/`vaf-denominator` pattern works with `field: INFO` for INFO-level fields.

### Option C: Use a pre-existing VAF field

If your VCF already contains a VAF-like field, point to it directly instead of using `"tbc"`:

```yaml
variant-calls:
  callset:
    my_callset:
      path: "path/to/variant-calls.vcf.gz"
      vaf-field:
        field: FORMAT
        name: AF
      benchmark: giab-na12878
      genome-build: grch38
```

The workflow will **not** recalculate VAF; it reads the field directly during precision/recall analysis. You can verify the field name with:

```bash
bcftools view -h path/to/variant-calls.vcf.gz
```

### Priority & scope

- **Callset-level** `vaf-field` (under `variant-calls.callset`) takes priority over benchmark-level.
- **Benchmark-level** `vaf-field` (under `custom-benchmarks.<name>`) applies to all callsets that reference that benchmark and lack a callset-level override.
- The [presets](../workflow/resources/presets.yaml) file ships with example `vaf-field` settings for common SomaticVariant callers (e.g., `TVAF` from Mutect2).

### Required Snakemake command

VAF calculation is triggered automatically — no extra Snakemake arguments are needed. However, to ensure VAF-stratified outputs are generated, make sure your benchmark definition includes `vaf-field` (either at callset or benchmark level) and request VAF-stratified tables:

```bash
snakemake -s workflow/Snakefile --use-conda -j 8
```
