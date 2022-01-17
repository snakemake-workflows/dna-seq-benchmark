rule replace_contig_names:
    input:
        calls=get_raw_callset,
        repl_file=get_rename_contig_file
    output:
        "results/normalized-variants/{callset}.replaced-contigs.bcf"
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools annotate {input.calls} --rename-chrs {input.repl_file} "
        "-Ob -o {output} 2> {log}"


rule normalize_calls:
    input:
        get_callset,
        genome="resources/reference/genome.fasta",
        genome_index="resources/reference/genome.fasta.fai",
    output:
        "results/normalized-variants/{callset}.vcf.gz",
    params:
        lambda w, input: f"--atomize -f {input.genome} --rm-dup exact -Oz",
    log:
        "logs/normalize-calls/{callset}.log",
    conda:
        "../envs/tools.yaml"
    wrapper:
        "0.80.2/bio/bcftools/norm"


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
        "bedtools intersect -b {input.regions} -a <(bcftools view {input.variants}) -wa -f 1.0 -header | bcftools view -Oz > {output} 2> {log}"


use rule stratify_truth as stratify_results with:
    input:
        variants="results/normalized-variants/{callset}.vcf.gz",
        regions=get_test_regions,
    output:
        "results/stratified-variants/{callset}/{cov}.vcf.gz",
    log:
        "logs/stratify-results/{callset}/{cov}.log",


rule bcftools_index:
    input:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz",
    output:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz.csi",
    log:
        "logs/bcftools-index/{benchmark}.truth.{cov}.log",
    wrapper:
        "0.80.1/bio/bcftools/index"


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
        prefix=get_io_prefix(lambda input, output: output[0]),
        engine="vcfeval",
    log:
        "logs/happy/{callset}/{cov}.log",
    wrapper:
        "0.80.1/bio/hap.py/hap.py"


rule collect_stratifications:
    input:
        expand("results/happy/{{callset}}/{cov}/report.summary.csv", cov=coverages),
    output:
        "results/report/{callset}.tsv",
    params:
        coverages=coverages,
    log:
        "logs/collect-stratifications/{callset}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-stratifications.py"


rule plot_precision_recall:
    input:
        "results/report/{callset}.tsv",
    output:
        "results/report/{callset}.plot.svg",
    params:
        cov_labels=get_plot_cov_labels(),
    log:
        "logs/plot-precision-recall/{callset}.log",
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/plot-precision-recall.py.ipynb"
