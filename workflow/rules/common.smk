from urllib.parse import urlparse

import yaml

with open(workflow.source_path("../resources/presets.yaml")) as presets:
    presets = yaml.load(presets, Loader=yaml.SafeLoader)

benchmarks = presets["benchmarks"]
genomes = presets["genomes"]

benchmarks.update(config.get("custom-benchmarks", dict()))

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


def get_archive_url(wildcards):
    genome = genomes[wildcards.genome]
    return genome["archive"]


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
    truth = genome["truth"][get_genome_build()]
    if isinstance(truth, dict):
        truth = truth[wildcards.truthset]
    if input.archive:
        return f"{input.archive}/{truth}"
    else:
        return truth


def get_truthsets(csi=False):
    def inner(wildcards):
        genome = genomes[wildcards.genome]
        truthsets = genome["truth"][get_genome_build()]
        return expand(
            "resources/variants/{genome}/{truthset}.truth.bcf",
            genome=wildcards.genome,
            truthset=truthsets,
        )

    return inner


def get_confidence_bed_cmd(wildcards, input):
    genome = genomes[wildcards.genome]
    bed = genome["confidence-regions"][get_genome_build()]

    if input.archive:
        return f"cat {input.archive}/{bed}"
    else:
        return f"curl --insecure -L {bed}"


def get_genome_build():
    if config.get("grch37"):
        return "grch37"
    else:
        return "grch38"


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


def get_target_regions(wildcards):
    benchmark = get_benchmark(wildcards.benchmark)
    if "target-regions" in benchmark:
        return "resources/regions/{benchmark}/target-regions.bed"
    else:
        return []


def get_target_regions_intersect_statement(wildcards, input):
    if input.target:
        return f"bedtools intersect -a /dev/stdin -b {input.target} |"
    else:
        return ""


def get_liftover_statement(wildcards, input, output):
    benchmark = get_benchmark(wildcards.benchmark)

    if benchmark["grch37"] and not config.get("grch37"):
        return f"| liftOver /dev/stdin {input.liftover} {output} /dev/null"
    else:
        return f"> {output}"


def get_read_limit_param(wildcards, input):
    if config.get("limit-reads"):
        return "| head -n 110000"  # a bit more than 100000 reads because we also have the header
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
    genome_name = get_benchmark(wildcards.benchmark)["genome"]
    genome = genomes[genome_name]

    truthset = genome["truth"][get_genome_build()]
    if isinstance(truthset, str):
        return f"resources/variants/{genome_name}/all.truth.bcf"
    else:
        return f"resources/variants/{genome_name}.merged.truth.bcf"


def get_stratified_truth(suffix=""):
    def inner(wildcards):
        benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
        return (f"results/variants/{benchmark}.truth.cov-{{cov}}.vcf.gz{suffix}",)

    return inner


def get_confidence_regions(wildcards):
    benchmark = get_benchmark(wildcards.benchmark)
    return f"resources/regions/{benchmark['genome']}.confidence-regions.bed"


def get_test_regions(wildcards):
    benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
    return f"resources/regions/{benchmark}/test-regions.cov-{{cov}}.bed"


def get_rename_contig_file(wildcards):
    return config["variant-calls"][wildcards.callset].get("rename-contigs")


def get_norm_params(wildcards, input):
    target = ""
    if config.get("limit-reads"):
        target = "--targets 1"
    return f"--atomize -f {input.genome} --check-ref s --rm-dup exact -Oz {target}"


def get_nonempty_coverages(wildcards):
    benchmark = config["variant-calls"][wildcards.callset]["benchmark"]

    def isempty(cov):
        with checkpoints.stat_truth.get(benchmark=benchmark, cov=cov).output[
            0
        ].open() as f:
            stat = json.load(f)
        return stat["isempty"]

    return [cov for cov in coverages if not isempty(cov)]


def get_collect_stratifications_input(wildcards):
    import json

    return expand(
        "results/happy/{{callset}}/{cov}/report.summary.csv",
        cov=get_nonempty_coverages(wildcards),
    )


if "variant-calls" in config:

    wildcard_constraints:
        callset="|".join(config["variant-calls"]),
