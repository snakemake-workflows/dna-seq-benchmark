repl_chr = "s/chr//"

if config["custom-reads"]["activate"]:
    reads = config["custom-reads"]["fastqs"]
    assert len(reads) == 2, "Expecting paired end custom reads in two FASTQ files"
else:
    reads = expand("resources/reads/reads.{read}.fq", read=[1, 2])

coverages = {
    "low": 1,
    "medium": 10,
    "high": 30,
}


def get_mosdepth_quantize():
    return ":".join(map(str, sorted(coverages.values()))) + ":"


def get_plot_cov_labels():
    def label(name):
        lower, upper = get_cov_interval(name)
        if upper:
            return f"{lower}-{upper-1}"
        return f"≥{lower}"

    return {name: label(name) for name in coverages}


def get_io_prefix(getter):
    def inner(wildcards, input, output):
        return getter(input, output).split(".")[0]

    return inner


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
    return config["variant-calls"][wildcards.callset]


def get_target_bed_statement():
    if config["custom-reads"]["activate"]:
        bed = config["custom-reads"]["target-regions"]
        return f"cat {bed}"
    else:
        return (
            "curl --insecure -L"
            " ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/TruSeq_exome_targeted_regions.hg19.bed"
        )


def get_liftover_statement(wildcards, input, output):
    if not config["custom-reads"]["activate"] or config["custom-reads"]["hg19"]:
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
        return f"bedtools -a /dev/stdin -b {input.limit_regions} |"
    else:
        return ""


wildcard_constraints:
    callset="|".join(config["variant-calls"]),
