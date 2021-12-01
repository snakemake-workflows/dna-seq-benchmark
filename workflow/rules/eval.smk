rule normalize_calls:
    input:
        config["results"],
        genome="reference/reference.fasta",
        genome_index="reference/reference.fasta.fai",
    output:
        "normalized-results/all.vcf.gz"
    params:
        lambda w, input: f"--atomize -f {input.genome} --rm-dup exact -Oz"
    conda:
        "../envs/tools.yaml"
    wrapper:
        "0.80.2/bio/bcftools/norm"


rule stratify_truth:
    input:
        variants="benchmark/truth.vcf",
        regions="benchmark/test-regions.cov-{cov}.bed",
    output:
        "benchmark/truth.cov-{cov}.vcf.gz",
    log:
        "logs/stratify-truth.{cov}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "bedtools intersect -b {input.regions} -a <(bcftools view {input.variants}) -wa -f 1.0 -header | bcftools view -Oz > {output} 2> {log}"


use rule stratify_truth as stratify_results with:
    input:
        variants="normalized-results/all.vcf.gz",
        regions="benchmark/test-regions.cov-{cov}.bed",
    output:
        "stratified-results/{cov}.vcf.gz",
    log:
        "logs/stratify-results.{cov}.log",


rule bcftools_index:
    input:
        "benchmark/truth.cov-{cov}.vcf.gz",
    output:
        "benchmark/truth.cov-{cov}.vcf.gz.csi",
    wrapper:
        "v0.80.1/bio/bcftools/index"


rule benchmark_variants:
    input:
        truth="benchmark/truth.cov-{cov}.vcf.gz",
        truth_idx="benchmark/truth.cov-{cov}.vcf.gz.csi",
        query="stratified-results/{cov}.vcf.gz",
        genome="reference/reference.fasta",
        genome_index="reference/reference.fasta.fai",
    output:
        multiext(
            "report-cov-{cov}",
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
        "logs/happy.{cov}.log",
    wrapper:
        "0.80.1/bio/hap.py/hap.py"


rule collect_stratifications:
    input:
        expand("report-cov-{cov}.summary.csv", cov=coverages),
    output:
        "report.tsv",
    params:
        coverages=coverages,
    log:
        "logs/collect-stratifications.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-stratifications.py"


rule plot_precision_recall:
    input:
        "report.tsv"
    output:
        "report.plot.svg"
    params:
        cov_labels=get_plot_cov_labels()
    log:
        "logs/plot-precision-recall.log"
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/plot-precision-recall.py.ipynb"