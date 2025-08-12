rule calc_precision_recall:
    input:
        calls="results/vcfeval/{callset}/{cov}/output.vcf.gz",
        idx="results/vcfeval/{callset}/{cov}/output.vcf.gz.tbi",
        common_src=common_src,
        truth=get_stratified_truth(),
        truth_idx=get_stratified_truth(".tbi"),
        query="results/stratified-variants/{callset}/{cov}.vcf.gz",
        query_index="results/stratified-variants/{callset}/{cov}.vcf.gz.tbi",
    output:
        "results/precision-recall/callsets/{callset}/{cov}.{vartype}.{mode}.tsv",
    log:
        "logs/calc-precision-recall/{callset}/{cov}/{vartype}.{mode}.log",
    params:
        vaf_fields=get_vaf_fields,
        vaf_status=get_vaf_status,
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/calc-precision-recall.py"


rule collect_stratifications:
    input:
        get_collect_stratifications_input,
    output:
        "results/precision-recall/callsets/{callset}.{vartype}.{mode}.tsv",
    params:
        coverages=get_nonempty_coverages,
        coverage_lower_bounds=get_coverages,
    log:
        "logs/collect-stratifications/{callset}/{vartype}.{mode}.log",
    conda:
        "../envs/stats.yaml"
    # We want this to be determined before FP/FN collection in order to avoid memory
    # issues with callsets that do not match the truth at all.
    priority: 2
    script:
        "../scripts/collect-stratifications.py"


rule collect_precision_recall:
    input:
        tables=get_collect_precision_recall_input,
    output:
        "results/precision-recall/benchmarks/{benchmark}.{vartype}.{mode}.tsv",
    params:
        callsets=lambda w: get_benchmark_callsets(w.benchmark),
        labels=get_collect_precision_recall_labels,
        vaf=get_vaf_status,
    log:
        "logs/collect-precision-recall/{benchmark}/{vartype}.{mode}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-precision-recall.py"


rule report_precision_recall:
    input:
        config=workflow.source_path(
            "../resources/datavzrd/precision-recall-config.yte.yaml"
        ),
        tables=get_report_precision_recall_input,
    output:
        report(
            directory("results/report/precision-recall/{benchmark}/{vartype}"),
            htmlindex="index.html",
            category="precision/recall",
            labels={
                "benchmark": "{benchmark}",
                "vartype": "{vartype}",
            },
        ),
    log:
        "logs/datavzrd/precision-recall/{benchmark}/{vartype}.log",
    params:
        somatic=get_somatic_status,
        vaf=get_vaf_status,
        high_coverage=get_high_coverage_status,
        genome=get_genome_name,
        version=get_genome_version,
    wrapper:
        "v7.0.0/utils/datavzrd"


rule extract_fp_fn:
    input:
        calls="results/vcfeval/{callset}/{cov}/output.vcf.gz",
        common_src=common_src,
    output:
        "results/fp-fn/callsets/{callset}/{cov}.{classification}.tsv",
    log:
        "logs/extract-fp-fn/{callset}/{cov}.{classification}.log",
    conda:
        "../envs/vembrane.yaml"
    script:
        "../scripts/extract-fp-fn.py"


# TODO: if one of the input callsets has all sites as FN, the resulting merged table
# becomes empty. This needs to be fixed, but is not super urgent because it is not realistic to
# happen in reality.
rule collect_fp_fn:
    input:
        tables=get_collect_fp_fn_input,
    output:
        main="results/fp-fn/genomes/{genome}/{cov}/{classification}/main.tsv",
        dependency_sorting=directory(
            "results/fp-fn/genomes/{genome}/{cov}/{classification}/dependency-sorting"
        ),
    params:
        callsets=get_collect_fp_fn_callsets,
        labels=get_collect_fp_fn_labels,
        label_names=lambda w: get_callsets_labels(get_genome_callsets(w.genome)),
        max_entries=config.get("max-fp-fn-entries", 300),
    log:
        "logs/collect-fp-fn/{genome}/{cov}/{classification}.log",
    conda:
        "../envs/stats.yaml"
    # This has to happen after precision/recall has been computed, otherwise we risk
    # extremely high memory usage if a callset does not match the truth at all.
    priority: 1
    script:
        "../scripts/collect-fp-fn.py"


rule collect_stratifications_fp_fn:
    input:
        get_collect_stratifications_fp_fn_input,
    output:
        "results/fp-fn/callsets/{callset}.{classification}.tsv",
    params:
        coverages=get_nonempty_coverages,
        coverage_lower_bounds=get_coverages,
    log:
        "logs/fp-fn/callsets/{callset}.{classification}.log",
    conda:
        "../envs/stats.yaml"
    # This has to happen after precision/recall has been computed, otherwise we risk
    # extremely high memory usage if a callset does not match the truth at all.
    priority: 1
    script:
        "../scripts/collect-stratifications-fp-fn.py"


rule collect_fp_fn_benchmark:
    input:
        tables=get_collect_fp_fn_benchmark_input,
    output:
        "results/fp-fn/benchmarks/{benchmark}.{classification}.tsv",
    params:
        callsets=lambda w: get_benchmark_callsets(w.benchmark),
    log:
        "logs/fp-fn/benchmarks/{benchmark}.{classification}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-fp-fn-benchmarks.py"


rule filter_shared_fn:
    input:
        fn="results/fp-fn/benchmarks/{benchmark}.fn.tsv",
    output:
        shared_fn="results/fp-fn/benchmarks/{benchmark}.shared_fn.tsv",
    log:
        "logs/filter-shared-fn/{benchmark}/{benchmark}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/filter-shared-variants.py"


# TODO: Add rule to filter shared FP variants


rule filter_unique_fn:
    input:
        fn="results/fp-fn/benchmarks/{benchmark}.fn.tsv",
    output:
        "results/fp-fn/benchmarks/{benchmark}/{callset}.unique_fn.tsv",
    params:
        variant_type="fn",
    log:
        "logs/filter-unique-variants/{benchmark}/{callset}.unique_fn.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/filter-unique-variants.py"


rule filter_unique_fp:
    input:
        fp="results/fp-fn/benchmarks/{benchmark}.fp.tsv",
    output:
        "results/fp-fn/benchmarks/{benchmark}/{callset}.unique_fp.tsv",
    params:
        variant_type="fp",
    log:
        "logs/filter-unique-variants/{benchmark}/{callset}.unique_fp.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/filter-unique-variants.py"


rule write_shared_fn_vcf:
    input:
        benchmark_table="results/fp-fn/benchmarks/{benchmark}.shared_fn.tsv",
        base_vcf=get_benchmark_renamed_truth,
        base_vcf_index=get_benchmark_renamed_truth_index,
    output:
        "results/fp-fn/vcf/{benchmark}/{benchmark}.shared_fn.vcf.gz",
    log:
        "logs/write-shared-fn-vcf/{benchmark}/{benchmark}.shared_fn.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/write-fp-fn-vcf.py"


rule write_unique_fn_vcf:
    input:
        benchmark_table="results/fp-fn/benchmarks/{benchmark}/{callset}.unique_fn.tsv",
        base_vcf=get_benchmark_renamed_truth,
        base_vcf_index=get_benchmark_renamed_truth_index,
    output:
        "results/fp-fn/vcf/{benchmark}/{callset}.unique_fn.vcf.gz",
    log:
        "logs/write-unique-fp-vcf/{benchmark}/{callset}.unique_fn.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/write-fp-fn-vcf.py"


rule write_unique_fp_vcf:
    input:
        benchmark_table="results/fp-fn/benchmarks/{benchmark}/{callset}.unique_fp.tsv",
        base_vcf="results/normalized-variants/{callset}.vcf.gz",
        base_vcf_index="results/normalized-variants/{callset}.vcf.gz.tbi",
    output:
        "results/fp-fn/vcf/{benchmark}/{callset}.unique_fp.vcf.gz",
    log:
        "logs/write-unique-fp-vcf/{benchmark}/{callset}.unique_fp.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/write-fp-fn-vcf.py"


rule report_fp_fn:
    input:
        main_dataset="results/fp-fn/genomes/{genome}/{cov}/{classification}/main.tsv",
        dependency_sorting_datasets="results/fp-fn/genomes/{genome}/{cov}/{classification}/dependency-sorting",
        config=workflow.source_path("../resources/datavzrd/fp-fn-config.yte.yaml"),
    output:
        report(
            directory("results/report/fp-fn/genomes/{genome}/{cov}/{classification}"),
            htmlindex="index.html",
            category="{classification} variants per genome",
            subcategory=lambda w: w.genome,
            labels=lambda w: {
                "coverage": w.cov,
                "genome": w.genome,
            },
        ),
    log:
        "logs/datavzrd/fp-fn/{genome}/{cov}/{classification}.log",
    params:
        labels=lambda w: get_callsets_labels(get_genome_callsets(w.genome)),
        version=get_genome_version,
    wrapper:
        "v7.0.0/utils/datavzrd"


rule report_fp_fn_callset:
    input:
        table="results/fp-fn/callsets/{callset}.{classification}.tsv",
        config=workflow.source_path(
            "../resources/datavzrd/fp-fn-per-callset-config.yte.yaml"
        ),
    output:
        report(
            directory("results/report/fp-fn/callsets/{callset}/{classification}"),
            htmlindex="index.html",
            category="{classification} variants per benchmark",
            subcategory=lambda w: config["variant-calls"][w.callset]["benchmark"],
            labels=lambda w: {
                "callset": w.callset,
            },
        ),
    log:
        "logs/datavzrd/fp-fn/{callset}/{classification}.log",
    params:
        labels=lambda w: get_callsets_labels(
            get_benchmark_callsets(config["variant-calls"][w.callset]["benchmark"])
        ),
        genome=get_genome_name,
        version=get_genome_version,
        somatic=get_somatic_status,
        high_coverage=get_high_coverage_status,
    wrapper:
        "v7.0.0/utils/datavzrd"


# TODO: Add rules to include unique and shared fp / fn variants in the report
