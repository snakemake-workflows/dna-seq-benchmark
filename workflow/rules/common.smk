from urllib.parse import urlparse

import yaml
import json

if "all" in config.get("variant-calls", dict()):
    raise ValueError(
        "A callset given in the variant-calls section of the config may not be called 'all'. "
        "Please choose a different name."
    )

with open(workflow.source_path("../resources/presets.yaml")) as presets:
    presets = yaml.load(presets, Loader=yaml.SafeLoader)

benchmarks = presets["benchmarks"]
genomes = presets["genomes"]
callsets = config.get("variant-calls", dict())

benchmarks.update(config.get("custom-benchmarks", dict()))
used_benchmarks = {callset["benchmark"] for callset in callsets.values()}

used_callsets = {callset for callset in callsets.keys()}

used_genomes = {benchmarks[benchmark]["genome"] for benchmark in used_benchmarks}

wildcard_constraints:
    benchmark="|".join(benchmarks)

if (
    any(callset["benchmark"] == "giab-NA12878-exome" for callset in callsets.values())
    and config["reference-genome"] != "grch38"
):
    raise ValueError(
        "grch38 must be set as reference-genome in the config if giab-NA12878-exome benchmark is used"
    )


repl_chr = "s/chr//"


low_coverages = {
    "lc_low": 1,
    "lc_medium": 10,
    "lc_high": 30,
}

high_coverages = {
    "hc_very_low": 1,
    "hc_low": 10,
    "hc_medium": 30,
    "hc_upper_medium": 50,
    "hc_high": 100,
    "hc_very_high": 300,
}


common_src = [
    workflow.source_path("../scripts/common/__init__.py"),
    workflow.source_path("../scripts/common/classification.py"),
]


def get_reference_genome_build():
    if "reference-genome" not in config:
        raise ValueError("Missing required configuration: reference-genome")
    if config["reference-genome"] == "grch37":
        return "GRCh37"
    elif config["reference-genome"] == "grch38":
        return "GRCh38"
    else:
        raise ValueError(
            f"Invalid reference genome build: {config['reference-genome']}. Must be one of: grch37, grch38"
        )


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


def get_mosdepth_quantize(wildcards):
    coverages = get_coverages(wildcards)
    return ":".join(map(str, sorted(coverages.values()))) + ":"


def get_plot_cov_labels():  # TODO check if ever used anywhere
    def label(name):
        lower, upper = get_cov_interval(name)
        if upper:
            return f"{lower}-{upper-1}"
        return f"≥{lower}"

    return {name: label(name) for name in low_coverages}


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
        if csi:
            return expand(
                "resources/variants/{genome}/{truthset}.truth.bcf.csi",
                genome=wildcards.genome,
                truthset=truthsets,
            )
        else:
            return expand(
                "resources/variants/{genome}/{truthset}.truth.bcf",
                genome=wildcards.genome,
                truthset=truthsets,
            )

    return inner


def get_confidence_bed_cmd(wildcards, input):
    genome = genomes[wildcards.genome]
    bed = genome["confidence-regions"][get_genome_build()]

    unpack_cmd = "| zcat " if bed.endswith(".gz") else ""

    if input.archive:
        return f"cat {input.archive}/{bed} {unpack_cmd}"
    if is_local_file(bed):
        return f"cat {bed} {unpack_cmd}"
    else:
        return f"curl --insecure -L {bed} {unpack_cmd}"


def get_genome_build():
    if "reference-genome" not in config:
        raise ValueError("Missing required configuration: reference-genome")
    return config["reference-genome"]


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
    coverages = get_coverages(wildcards)
    lower, upper = get_cov_interval(wildcards.cov, coverages)
    if upper:
        return f"{lower}:{upper}"
    return f"{lower}:inf"


def get_cov_interval(name, coverages):
    threshold = coverages[name]
    upper_bound = None

    greater = [cov for cov in coverages.values() if cov > threshold]
    if greater:
        upper_bound = min(greater)

    return threshold, upper_bound


def get_callset(wildcards):
    callset = config["variant-calls"][wildcards.callset]
    if get_somatic_status(wildcards):
        return "results/normalized-variants/{callset}.gt-added.vcf.gz"
    elif callset.get("rename-contigs", False) != False:
        return "results/normalized-variants/{callset}.replaced-contigs.vcf.gz"
    elif callset["genome-build"] == "grch37":
        return "results/normalized-variants/{callset}.lifted.vcf.gz"
    else:
        return get_raw_callset(wildcards)


def get_callset_correct_contigs(wildcards):
    callset = config["variant-calls"][wildcards.callset]
    if "rename-contigs" in callset:
        return "results/normalized-variants/{callset}.replaced-contigs.vcf.gz"
    elif callset["genome-build"] == "grch37":
        return "results/normalized-variants/{callset}.lifted.vcf.gz"
    else:
        return get_raw_callset(wildcards)


def get_callset_correct_contigs_liftover(wildcards):
    callset = config["variant-calls"][wildcards.callset]
    if callset["genome-build"] == "grch37":
        return "results/normalized-variants/{callset}.lifted.vcf.gz"
    elif callset.get("rename-contigs", False) != False:
        return "results/normalized-variants/{callset}.replaced-contigs.vcf.gz"
    else:
        return get_raw_callset(wildcards)


def get_raw_callset(wildcards):
    callset = config["variant-calls"][wildcards.callset]
    return callset["path"]


def is_local_file(path):
    return urlparse(path).scheme == ""


def get_target_bed_input(wildcards):
    target_bed = get_benchmark(wildcards.benchmark)["target-regions"]
    if is_local_file(target_bed):
        return target_bed
    else:
        return []


def get_target_bed_statement(wildcards):
    target_bed = get_benchmark(wildcards.benchmark)["target-regions"]

    unpack_cmd = "| zcat " if target_bed.endswith(".gz") else ""

    if is_local_file(target_bed):
        return f"cat {target_bed} {unpack_cmd}"
    else:
        return f"curl --insecure -L {target_bed} {unpack_cmd}"


def get_target_regions(wildcards):
    if hasattr(wildcards, "benchmark"):
        benchmark_dict = get_benchmark(wildcards.benchmark)
        if "target-regions" in benchmark_dict:
            return f"resources/regions/{wildcards.benchmark}/target-regions.bed"
        else:
            return []
    else:
        benchmark_name = config["variant-calls"][wildcards.callset]["benchmark"]
        benchmark = get_benchmark(benchmark_name)
        if "target-regions" in benchmark:
            return f"resources/regions/{benchmark_name}/target-regions.bed"
        else:
            return []


def intersect_calls(wildcards):
    if get_target_regions(wildcards) == []:
        return False
    else:
        return True


def get_target_regions_intersect_statement(wildcards, input):
    if input.target:
        return f"bedtools intersect -a /dev/stdin -b {input.target} |"
    else:
        return ""


def get_liftover_statement(wildcards, input, output):
    benchmark = get_benchmark(wildcards.benchmark)
    if benchmark["grch37"] and not config["reference-genome"] == "grch37":
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


def get_genome_truth(wildcards):
    genome_name = wildcards.genome
    genome = genomes[genome_name]

    truthset = genome["truth"][get_genome_build()]
    if isinstance(truthset, str):
        return f"resources/variants/{genome_name}/all.truth.bcf"
    else:
        return f"resources/variants/{genome_name}.merged.truth.bcf"


def get_benchmark_truth(wildcards):
    genome = get_benchmark(wildcards.benchmark)["genome"]
    if get_somatic_status(wildcards):
        return f"resources/variants/{genome}/all.truth.format-added.vcf.gz"
    else:
        return f"resources/variants/{genome}/all.truth.norm.bcf"


def get_stratified_truth(suffix=""):
    def inner(wildcards):
        benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
        return f"results/variants/{benchmark}.truth.cov-{{cov}}.vcf.gz{suffix}"

    return inner


def get_confidence_regions(wildcards):
    benchmark = get_benchmark(wildcards.benchmark)
    return f"resources/regions/{benchmark['genome']}.confidence-regions.bed"


def get_test_regions(wildcards):
    benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
    return f"resources/regions/{benchmark}/test-regions.cov-{{cov}}.bed"


def get_rename_contig_file(wildcards):
    if config["variant-calls"][wildcards.callset][
        "genome-build"
    ] == "grch37" and config["variant-calls"][wildcards.callset].get(
        "rename-contigs", False
    ):
        return workflow.source_path(
            "../resources/rename-contigs/grch37_ucsc2ensembl.txt"
        )
    if config["variant-calls"][wildcards.callset][
        "genome-build"
    ] == "grch38" and config["variant-calls"][wildcards.callset].get(
        "rename-contigs", False
    ):
        return workflow.source_path(
            "../resources/rename-contigs/grch38_ucsc2ensembl.txt"
        )
    else:
        return config["variant-calls"][wildcards.callset].get("rename-contigs", False)


def get_callset_subcategory(wildcards):
    return config["variant-calls"][wildcards.callset].get("subcategory")


def get_norm_params(wildcards):
    target = ""
    if config.get("limit-reads"):
        target = "--targets 1"
    return f"--atomize --check-ref s --rm-dup exact {target}"


def get_mosdepth_input(bai=False):
    ext = ".bai" if bai else ""

    def inner(wildcards):
        benchmark = get_benchmark(wildcards.benchmark)
        bam = benchmark.get("bam")
        if bam:
            if "bai" in benchmark and bai:
                return benchmark["bai"]
            return bam + ext
        else:
            return f"results/read-alignments/{wildcards.benchmark}.dedup.bam{ext}"

    return inner


def _get_nonempty_coverages(callset):
    benchmark = config["variant-calls"][callset]["benchmark"]
    benchmark_dict = get_benchmark(benchmark)
    if benchmark_dict.get("high-coverage", False):
        coverages = high_coverages
    else:
        coverages = low_coverages

    def isempty(cov):
        with checkpoints.stat_truth.get(benchmark=benchmark, cov=cov).output[
            0
        ].open() as f:
            stat = json.load(f)
        return stat["isempty"]

    return [cov for cov in coverages if not isempty(cov)]


def get_nonempty_coverages(wildcards):
    return _get_nonempty_coverages(wildcards.callset)


def get_coverages(wildcards):
    if hasattr(wildcards, "benchmark"):
        high_cov_status = benchmarks[wildcards.benchmark].get("high-coverage")
    else:
        benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
        high_cov_status = benchmarks[benchmark].get("high-coverage")
    if high_cov_status:
        coverages = high_coverages
    else:
        coverages = low_coverages
    return coverages


def get_coverages_of_callset(callset):
    benchmark = config["variant-calls"][callset]["benchmark"]
    high_cov_status = benchmarks[benchmark].get("high-coverage")
    if high_cov_status:
        coverages = high_coverages
    else:
        coverages = low_coverages
    return coverages


def get_somatic_status(wildcards):
    if hasattr(wildcards, "benchmark"):
        return genomes[benchmarks[wildcards.benchmark]["genome"]].get("somatic")
    else:
        benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
        return genomes[benchmarks[benchmark]["genome"]].get("somatic")


def get_somatic_sample_name(wildcards):
    return config["variant-calls"][wildcards.callset]["tumor_sample_name"]


def get_somatic_flag(wildcards):
    if get_somatic_status(wildcards):
        sample_name_baseline = "truth"
        sample_name_callset = config["variant-calls"][wildcards.callset][
            "tumor_sample_name"
        ]  # get name tumor via config -> from dict
        somatic_flag = (
            f"--squash-ploidy --sample {sample_name_baseline},{sample_name_callset}"
        )
    else:
        somatic_flag = ""
    return somatic_flag


def get_vaf_fields(wildcards):
    vaf_callset = config["variant-calls"][wildcards.callset].get("vaf-field")

    benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
    vaf_benchmark = benchmarks[benchmark].get("vaf-field")

    # can return (None, None) if param not set
    return (vaf_callset, vaf_benchmark)


def get_vaf_status(wildcards):
    vaf_benchmark = benchmarks[wildcards.benchmark].get("vaf-field")
    if vaf_benchmark is None:
        return False
    else:
        callsets = get_benchmark_callsets(wildcards.benchmark)
        vaf_callsets = [
            config["variant-calls"][callset].get("vaf-field") for callset in callsets
        ]
        if any(vaf_callset is not None for vaf_callset in vaf_callsets):
            return True
        else:
            return False


def get_high_coverage_status(wildcards):
    if hasattr(wildcards, "benchmark"):
        benchmark = get_benchmark(wildcards.benchmark)
        return benchmark.get("high-coverage", False)
    else:
        benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
        return benchmarks[benchmark].get("high-coverage", False)


def get_collect_stratifications_input(wildcards):
    import json

    return expand(
        "results/precision-recall/callsets/{{callset}}/{cov}.{{vartype}}.tsv",
        cov=get_nonempty_coverages(wildcards),
    )


def get_collect_stratifications_fp_fn_input(wildcards):
    return expand(
        "results/fp-fn/callsets/{{callset}}/{cov}.{{classification}}.tsv",
        cov=get_nonempty_coverages(wildcards),
    )


def get_fp_fn_reports(wildcards):
    for genome in used_genomes:
        yield from expand(
            "results/report/fp-fn/genomes/{genome}/{cov}/{classification}",
            genome=genome,
            cov={
                cov
                for callset in get_genome_callsets(genome)
                for cov in _get_nonempty_coverages(callset)
            },
            classification=["fp", "fn"],
        )


def get_fp_fn_reports_benchmarks(wildcards):
    for genome in used_genomes:
        yield from expand(
            "results/report/fp-fn/benchmarks/{benchmark}/{classification}",
            benchmark={benchmark for benchmark in used_benchmarks},
            classification=["fp", "fn"],
        )


def get_benchmark_callsets(benchmark):
    return [
        callset
        for callset, entries in config["variant-calls"].items()
        if entries["benchmark"] == benchmark
    ]


def get_collect_precision_recall_input(wildcards):
    callsets = get_benchmark_callsets(wildcards.benchmark)
    return expand(
        "results/precision-recall/callsets/{callset}.{{vartype}}.tsv", callset=callsets
    )


def get_collect_fp_fn_benchmark_input(wildcards):
    callsets = get_benchmark_callsets(wildcards.benchmark)
    return expand(
        "results/fp-fn/callsets/{callset}.{{classification}}.tsv", callset=callsets
    )


def get_genome_name(wildcards):
    if hasattr(wildcards, "benchmark"):
        return get_benchmark(wildcards.benchmark).get("genome")
    if hasattr(wildcards, "callset"):
        benchmark = config["variant-calls"][wildcards.callset]["benchmark"]
        return get_benchmark(benchmark).get("genome")
    else:
        return wildcards.genome


def get_genome_version(wildcards):
    return genomes[get_genome_name(wildcards)].get("version")


def get_genome_callsets(genome):
    return sorted(
        callset
        for callset, entries in config["variant-calls"].items()
        if benchmarks[entries["benchmark"]]["genome"] == genome
    )


def get_callsets_labels(callsets):
    return sorted(
        set(
            label
            for callset in callsets
            for label in config["variant-calls"][callset].get("labels", [])
        )
    )


def get_callset_label_entries(callsets):
    labels = get_callsets_labels(callsets)

    def get_label_row(label):
        def labelval(label, callset):
            return config["variant-calls"][callset].get("labels", dict()).get(label)

        return {
            callset: labelval(label, callset)
            for callset in callsets
            if labelval(label, callset)
        }

    return [get_label_row(label) for label in labels]


def get_collect_fp_fn_callsets(wildcards):
    callsets = get_genome_callsets(wildcards.genome)
    callsets = [
        callset
        for callset in callsets
        if wildcards.cov in get_coverages_of_callset(callset)
    ]
    return callsets


def get_collect_fp_fn_input(wildcards):
    callsets = get_collect_fp_fn_callsets(wildcards)
    return expand(
        "results/fp-fn/callsets/{callset}/{{cov}}.{{classification}}.tsv",
        callset=callsets,
    )


def get_collect_fp_fn_labels(wildcards):
    callsets = get_collect_fp_fn_callsets(wildcards)
    return get_callset_label_entries(callsets)


def get_collect_precision_recall_labels(wildcards):
    callsets = get_benchmark_callsets(wildcards.benchmark)
    return get_callset_label_entries(callsets)


if "variant-calls" in config:

    wildcard_constraints:
        callset="|".join(config["variant-calls"]),
        classification="fp|fn",
        comparison="genotype|existence",
        vartype="snvs|indels",
