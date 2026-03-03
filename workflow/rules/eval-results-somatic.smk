#### Precision & Recall ####

rule calc_precision_recall:
    input:
        calls="results/vcfeval/{callset}/{cov}/output.vcf.gz",
        idx="results/vcfeval/{callset}/{cov}/output.vcf.gz.tbi",
        common_src=common_src,
        truth=get_stratified_truth(),
        truth_idx=get_stratified_truth(".tbi"),
        query="results/stratified-variants/{callset}/{cov}.vcf.gz",
        query_idx="results/stratified-variants/{callset}/{cov}.vcf.gz.tbi",
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
        "v8.0.3/utils/datavzrd"
