# Snakemake workflow: dna-seq-benchmark

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/dna-seq-benchmark/workflows/Tests/badge.svg?branch=main)](https://github.com/snakemake-workflows/dna-seq-benchmark/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for benchmarking variant calling approaches with Genome in a Bottle (GIAB) data (and other custom benchmark datasets). The workflow uses a combination of bedtools, mosdepth, rtg-tools, pandas and datavzrd.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/dna-seq-benchmark.html).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) benchmark-giabsitory and its DOI (see above).

## Output

The workflow writes both final deliverables and intermediate files under `results/`.

### Final results (usually what you want)

- `results/report/precision-recall/<benchmark>/<snvs|indels>/`: interactive precision/recall reports
- `results/report/fp-fn/callsets/<callset>/<fp|fn>/`: per-callset FP/FN reports
- `results/report/fp-fn/genomes/<genome>/<coverage>/<fp|fn>/`: genome-level FP/FN reports
- `results/annotated/tsv/<benchmark>/`: annotated shared FN tables
- `results/annotated/tsv/<benchmark>/<callset>.unique_<fp|fn>.annotated.tsv`: annotated unique FP/FN tables

### Supporting analysis tables

- `results/precision-recall/benchmarks/`: aggregated precision/recall benchmark tables
- `results/fp-fn/benchmarks/`: aggregated FP/FN benchmark tables and shared/unique splits
- `results/fp-fn/vcf/`: VCFs generated from shared/unique FP/FN tables

### Intermediates and automatic cleanup

- Raw somatic extraction tables are written to `results/intermediate/fp-fn/raw/callsets/`.
- Several per-coverage and per-callset aggregation inputs are marked as Snakemake `temp()` outputs and are removed automatically once downstream targets are finished.
- If you want to keep all intermediates for debugging, run Snakemake with `--notemp`.
