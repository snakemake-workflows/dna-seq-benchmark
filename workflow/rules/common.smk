from urllib.parse import urlparse

benchmarks = dict(config.get("custom-benchmarks", dict()))

benchmarks["giab-na12878-exome"] = {
    "genome": "na12878",
    "bam-url": "ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam",
    "target-regions": "ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/TruSeq_exome_targeted_regions.hg19.bed",
    "grch37": True,
}

benchmarks["chm-eval"] = {
    "genome": "chm-eval",
    "bam-url": "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR134/ERR1341793/CHM1_CHM13_3.bam",
    "target-regions": None,
    "grch37": False,
}

genomes = {
    "na12878": {
        "truth": {
            "grch37": "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        },
        "confidence-regions": {
            "grch37": "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.bed",
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed",
        },
    },
    "chm-eval": {
        "archive": "https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar",
        "truth": {
            "grch38": "full.38.vcf.gz",
            "grch37": "full.37m.vcf.gz",
        },
        "confidence-regions": {
            "grch38": "full.38.bed.gz",
            "grch37": "full.37m.bed.gz",
        },
    },
}


if any(
    callset["benchmark"] == "giab-NA12878-exome"
    for _, callset in config.get("variant-calls", dict()).items()
) and config.get("grch37"):
    raise ValueError(
        "grch37 must be set to false in the config if giab-NA12878-exome benchmark is used"
    )


repl_chr = "s/chr//"


coverages = {
    "low": 1,
    "medium": 10,
    "high": 30,
}


def get_archive_input(wildcards):
    genome = genomes[wildcards.genome]
    if "archive" in genome:
        return f"resources/archives/{wildcards.genome}"
    else:
        return []


def get_benchmark_bam_url(wildcards):
    return get_benchmark(wildcards.benchmark)["bam-url"]


def get_bwa_input(wildcards):
    benchmark = get_benchmark(wildcards.benchmark)
    if "bam-url" in benchmark:
        return expand(
            "resources/reads/{benchmark}.{read}.fq",
            read=[1, 2],
            benchmark=wildcards.benchmark,
        )
    else:
        return benchmark["fastqs"]


def get_mosdepth_quantize():
    return ":".join(map(str, sorted(coverages.values()))) + ":"


def get_plot_cov_labels():
    def label(name):
        lower, upper = get_cov_interval(name)
        if upper:
            return f"{lower}-{upper-1}"
        return f"≥{lower}"

    return {name: label(name) for name in coverages}


def get_truth_url(wildcards, input):
    genome = genomes[wildcards.genome]
    truth = genome["truth"][wildcards.build]
    if input.archive:
        return f"{input.archive}/{truth}"
    else:
        return truth


def get_confidence_bed_cmd(wildcards, input):
    genome = genomes[wildcards.genome]
    bed = genome["confidence-regions"][wildcards.build]
    if input.archive:
        return f"cat {input.archive}/{bed}"
    else:
        return f"curl --insecure -L {bed}"


def get_io_prefix(getter):
    def inner(wildcards, input, output):

        return getter(input, output).split(".")[0]

    return inner


def get_happy_prefix(wildcards, output):
    runinfo_suffix = ".runinfo.json"
    for f in output:
        if f.endswith(runinfo_suffix):
            return f[: -len(runinfo_suffix)]


def get_cov_label(wildcards):
    lower, upper = get_cov_interval(wildcards.cov)
    if upper:
        return f"{lower}:{upper}"
    return f"{lower}:inf"


def get_cov_interval(name):
    threshold = coverages[name]
    upper_bound = None

    greater = [cov for cov in coverages.values() if cov > threshold]
    if greater:
        upper_bound = min(greater)

    return threshold, upper_bound


def get_callset(wildcards):
    callset = config["variant-calls"][wildcards.callset]
    if "rename-contigs" in callset:
        return "results/normalized-variants/{callset}.replaced-contigs.bcf"
    else:
        return get_raw_callset(wildcards)


def get_raw_callset(wildcards):
    callset = config["variant-calls"][wildcards.callset]
    return callset["path"]


def get_target_bed_statement(wildcards):
    target_bed = get_benchmark(wildcards.benchmark)["target-regions"]

    if urlparse(target_bed).scheme == "":
        return f"cat {target_bed}"
    else:
        return f"curl --insecure -L {target_bed}"


def get_liftover_statement(wildcards, input, output):
    benchmark = get_benchmark(wildcards.benchmark)

    if benchmark["grch37"] and not config.get("grch37"):
        return f"| liftOver /dev/stdin {input.liftover} {output} /dev/null"
    else:
        return f"> {output}"


def get_limit_regions():
    if config["limit-regions"]["activate"]:
        return config["limit-regions"]["bed"]
    else:
        return []


def get_read_limit_param(wildcards, input):
    if input.get("regions"):
        return f"-L {input.regions}"
    else:
        return ""


def get_limit_regions_intersect_statement(wildcards, input):
    if input.get("limit_regions"):
        return f"| bedtools intersect -a /dev/stdin -b {input.limit_regions}"
    else:
        return ""


def get_benchmark(benchmark):
    try:
        return benchmarks[benchmark]
    except KeyError:
        raise ValueError(
            f"Benchmark name {benchmark} does not occur in the custom-benchmarks section."
        )


def get_benchmark_truth(wildcards):
    genome = get_benchmark(wildcards.benchmark)["genome"]
    return f"resources/variants/{genome}.truth.vcf"


def get_stratified_truth(suffix=""):
    def inner(wildcards):
        benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
        return (f"results/variants/{benchmark}.truth.cov-{{cov}}.vcf.gz{suffix}",)

    return inner


def get_confidence_regions(wildcards):
    benchmark = get_benchmark(wildcards.benchmark)
    if benchmark["genome"] == "NA12878":
        return "resources/regions/NA12898.confidence-regions.bed"
    else:
        # TODO add CHM support
        raise ValueError(f"Unsupported genome {benchmark['genome']}")


def get_test_regions(wildcards):
    benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
    return f"resources/regions/{benchmark}/test-regions.cov-{{cov}}.bed"


def get_rename_contig_file(wildcards):
    return config["variant-calls"][wildcards.callset].get("rename-contigs")


if "variant-calls" in config:

    wildcard_constraints:
        callset="|".join(config["variant-calls"]),
