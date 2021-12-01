repl_chr = "s/chr//"
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


wildcard_constraints:
    callset="|".join(config["variant-calls"]),
