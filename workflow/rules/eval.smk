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
        "results/happy/{callset}/{cov}/report.vcf.gz"
    output:
        "results/precision-recall/callsets/{callset}/{cov}.tsv"
    log:
        "logs/calc-precision-recall/{callset}/{cov}.log"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/calc-precision-recall.py"


rule collect_stratifications:
    input:
        get_collect_stratifications_input,
    output:
        "results/precision-recall/callsets/{callset}.tsv",
    params:
        coverages=get_nonempty_coverages,
    log:
        "logs/collect-stratifications/{callset}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-stratifications.py"


rule collect_precision_recall:
    input:
        get_collect_precision_recall_input
    output:
        "results/precision-recall/benchmarks/{benchmark}.tsv"
    params:
        callsets=lambda w: get_benchmark_callsets(w.benchmark)
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-precision-recall.py"


rule render_precision_recall_report_config:
    input:
        dataset="results/precision-recall/benchmarks/{benchmark}.tsv",
        template=workflow.source_path("../resources/datavzrd/precision-recall-config.yte.yaml"),
    output:
        "results/datavzrd-config/precision-recall/{benchmark}/all.config.yaml",
    log:
        "logs/yte/datavzrd-config/precision-recall/{benchmark}.log",
    template_engine:
        "yte"


rule report_precision_recall:
    input:
        config="results/datavzrd-config/precision-recall/{benchmark}/all.config.yaml",
        table="results/precision-recall/benchmarks/{benchmark}.tsv"
    output:
        report(
            directory("results/report/precision-recall/{benchmark}"),
            htmlindex="index.html",
            category="precision/recall",
            labels={"benchmark": "{benchmark}"},
        ),
    log:
        "logs/datavzrd/precision-recall/{benchmark}.log",
    wrapper:
        "v1.17.0/utils/datavzrd"


rule collect_subsets:
    input:
        calls="results/happy/{callset}/{cov}/report.vcf.gz",
    output:
        "results/classified-subsets/{cov}/{callset}.{type,FP|FN}.tsv",
    log:
        "logs/vembrane/subsets/{cov}/{callset}.{type}.log",
    params:
        #filter=lambda w: '\'FORMAT["BD"]["QUERY"] == "FP"\'' if w.type == "FP" else '\'FORMAT["BD"]["TRUTH"] == "FN"\''
        filter=get_subset_filter,
    conda:
        "../envs/vembrane.yaml"
    shell:
        """
        (filtered=$(mktemp)
        bcftools norm -m-any {input.calls} | vembrane filter {params.filter:q} > $filtered
        vembrane table --header 'chromosome, position, ref_allele, alt_allele, true_genotype, predicted_genotype' \
        'CHROM, POS, REF, ALT, \
        "{{}}/{{}}".format(*sorted(FORMAT["GT"]["TRUTH"])) if FORMAT["GT"]["TRUTH"] is not NA else ".", \
        "{{}}/{{}}".format(*sorted(FORMAT["GT"]["QUERY"])) if FORMAT["GT"]["QUERY"] is not NA else "."' \
        $filtered > {output}
        rm $filtered) 2> {log}
        """


rule merge_subsets:
    input:
        get_merged_classified_subsets_input,
    output:
        "results/merged-classified-subsets/{genome}/{cov}/all.{type,FP|FN}.tsv",
    params:
        callsets=get_merged_classified_subsets_callsets,
    log:
        "logs/merge-substes/{genome}/{cov}/{type}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/merge-subsets.py"


rule render_subset_report_config:
    input:
        dataset="results/merged-classified-subsets/{genome}/{cov}/all.{type}.tsv",
        template=workflow.source_path("../resources/datavzrd/subset-config.yte.yaml"),
    output:
        "results/datavzrd-config/subsets/{genome}/{cov}/all.{type}.config.yaml",
    log:
        "logs/yte/datavzrd-config/subsets/{genome}/{cov}/{type}.log",
    template_engine:
        "yte"


rule report_subsets:
    input:
        dataset="results/merged-classified-subsets/{genome}/{cov}/all.{type}.tsv",
        config="results/datavzrd-config/subsets/{genome}/{cov}/all.{type}.config.yaml",
    output:
        report(
            directory("results/report/fp-fn/{genome}/{cov}/all.{type}"),
            htmlindex="index.html",
            category=lambda w: "false positives"
            if w.type == "FP"
            else "false negatives",
            subcategory=lambda w: w.genome,
            labels=lambda w: {"coverage": w.cov},
        ),
    log:
        "logs/datavzrd/fp-fn/{genome}/{cov}/{type}.log",
    wrapper:
        "v1.17.0/utils/datavzrd"
