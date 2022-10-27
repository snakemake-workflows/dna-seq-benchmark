rule rename_contigs:
    input:
        calls=get_raw_callset,
        repl_file=get_rename_contig_file,
    output:
        "results/normalized-variants/{callset}.replaced-contigs.bcf",
    log:
        "logs/rename-contigs/{callset}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools annotate {input.calls} --rename-chrs {input.repl_file} "
        "-Ob -o {output} 2> {log}"


rule normalize_calls:
    input:
        get_callset,
        ref="resources/reference/genome.fasta",
        ref_index="resources/reference/genome.fasta.fai",
    output:
        "results/normalized-variants/{callset}.vcf.gz",
    params:
        extra=get_norm_params,
    log:
        "logs/normalize-calls/{callset}.log",
    wrapper:
        "v1.9.0/bio/bcftools/norm"


rule stratify_truth:
    input:
        variants=get_benchmark_truth,
        regions="resources/regions/{benchmark}/test-regions.cov-{cov}.bed",
    output:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz",
    log:
        "logs/stratify-truth.{benchmark}.{cov}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bedtools intersect -b {input.regions} -a "
        "<(bcftools view {input.variants}) -wa -f 1.0 -header | "
        "bcftools view -Oz > {output}) 2> {log}"


use rule stratify_truth as stratify_results with:
    input:
        variants="results/normalized-variants/{callset}.vcf.gz",
        regions=get_test_regions,
    output:
        "results/stratified-variants/{callset}/{cov}.vcf.gz",
    log:
        "logs/stratify-results/{callset}/{cov}.log",


rule index_stratified_truth:
    input:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz",
    output:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz.csi",
    log:
        "logs/bcftools-index/{benchmark}.truth.{cov}.log",
    wrapper:
        "v1.7.2/bio/bcftools/index"


checkpoint stat_truth:
    input:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz",
    output:
        "results/variants/{benchmark}.truth.cov-{cov}.stats.json",
    log:
        "logs/stat-truth/{benchmark}.{cov}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/stat-truth.py"


rule benchmark_variants:
    input:
        truth=get_stratified_truth(),
        truth_idx=get_stratified_truth(".csi"),
        query="results/stratified-variants/{callset}/{cov}.vcf.gz",
        genome="resources/reference/genome.fasta",
        genome_index="resources/reference/genome.fasta.fai",
    output:
        multiext(
            "results/happy/{callset}/{cov}/report",
            ".runinfo.json",
            ".vcf.gz",
            ".summary.csv",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
            ".roc.Locations.INDEL.csv.gz",
            ".roc.Locations.INDEL.PASS.csv.gz",
            ".roc.Locations.SNP.csv.gz",
        ),
    params:
        prefix=get_happy_prefix,
        engine="vcfeval",
    log:
        "logs/happy/{callset}/{cov}.log",
    wrapper:
        "v1.7.2/bio/hap.py/hap.py"


rule calc_precision_recall:
    input:
        "results/happy/{callset}/{cov}/report.vcf.gz",
    output:
        snvs="results/precision-recall/callsets/{callset}/{cov}.{vartype}.tsv",
    log:
        "logs/calc-precision-recall/{callset}/{cov}/{vartype}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/calc-precision-recall.py"


rule collect_stratifications:
    input:
        get_collect_stratifications_input,
    output:
        "results/precision-recall/callsets/{callset}.{vartype}.tsv",
    params:
        coverages=get_nonempty_coverages,
        coverage_lower_bounds=coverages,
    log:
        "logs/collect-stratifications/{callset}/{vartype}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-stratifications.py"


rule collect_precision_recall:
    input:
        get_collect_precision_recall_input,
    output:
        "results/precision-recall/benchmarks/{benchmark}.{vartype}.tsv",
    params:
        callsets=lambda w: get_benchmark_callsets(w.benchmark),
        labels=get_collect_precision_recall_labels,
    log:
        "logs/collect-precision-recall/{benchmark}/{vartype}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-precision-recall.py"


rule render_precision_recall_report_config:
    input:
        dataset="results/precision-recall/benchmarks/{benchmark}.{vartype}.tsv",
        template=workflow.source_path(
            "../resources/datavzrd/precision-recall-config.yte.yaml"
        ),
    output:
        "results/datavzrd-config/precision-recall/{benchmark}/{vartype}.config.yaml",
    log:
        "logs/yte/datavzrd-config/precision-recall/{benchmark}/{vartype}.log",
    template_engine:
        "yte"


rule report_precision_recall:
    input:
        config="results/datavzrd-config/precision-recall/{benchmark}/{vartype}.config.yaml",
        table="results/precision-recall/benchmarks/{benchmark}.{vartype}.tsv",
    output:
        report(
            directory("results/report/precision-recall/{benchmark}/{vartype}"),
            htmlindex="index.html",
            category="precision/recall",
            labels={"benchmark": "{benchmark}", "vartype": "{vartype}"},
        ),
    log:
        "logs/datavzrd/precision-recall/{benchmark}/{vartype}.log",
    wrapper:
        "v1.17.4/utils/datavzrd"


rule extract_fp_fn:
    input:
        calls="results/happy/{callset}/{cov}/report.vcf.gz",
    output:
        "results/fp-fn/callsets/{cov}/{callset}/{classification}.tsv",
    log:
        "logs/extract-fp-fn/{cov}/{callset}/{classification}.log",
    conda:
        "../envs/vembrane.yaml"
    script:
        "../scripts/extract-fp-fn.py"


rule collect_fp_fn:
    input:
        get_collect_fp_fn_input,
    output:
        main="results/fp-fn/genomes/{genome}/{cov}/{classification}/main.tsv",
        dependency_sorting=directory(
            "results/fp-fn/genomes/{genome}/{cov}/{classification}/dependency-sorting"
        ),
    params:
        callsets=get_collect_fp_fn_callsets,
        labels=get_collect_fp_fn_labels,
        label_names=lambda w: get_callset_labels(get_genome_callsets(w.genome)),
    log:
        "logs/collect-fp-fn/{genome}/{cov}/{classification}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-fp-fn.py"


rule render_fp_fn_report_config:
    input:
        main_dataset="results/fp-fn/genomes/{genome}/{cov}/{classification}/main.tsv",
        dependency_sorting_datasets="results/fp-fn/genomes/{genome}/{cov}/{classification}/dependency-sorting",
        template=workflow.source_path("../resources/datavzrd/fp-fn-config.yte.yaml"),
    output:
        "results/datavzrd-config/fp-fn/{genome}/{cov}/{classification}.config.yaml",
    params:
        labels=lambda w: get_callset_labels(get_genome_callsets(w.genome)),
    log:
        "logs/yte/datavzrd-config/fp-fn/{genome}/{cov}/{classification}.log",
    template_engine:
        "yte"


rule report_fp_fn:
    input:
        main_dataset="results/fp-fn/genomes/{genome}/{cov}/{classification}/main.tsv",
        dependency_sorting_datasets="results/fp-fn/genomes/{genome}/{cov}/{classification}/dependency-sorting",
        config="results/datavzrd-config/fp-fn/{genome}/{cov}/{classification}.config.yaml",
    output:
        report(
            directory("results/report/fp-fn/{genome}/{cov}/{classification}"),
            htmlindex="index.html",
            category="{classification} variants",
            subcategory=lambda w: w.genome,
            labels=lambda w: {"coverage": w.cov},
        ),
    log:
        "logs/datavzrd/fp-fn/{genome}/{cov}/{classification}.log",
    wrapper:
        "v1.17.4/utils/datavzrd"
