import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

def load_data(f, coverage):
    d = pd.read_csv(f)
    d.insert(0, "coverage", coverage)
    return d

report = pd.concat(load_data(f, cov) for cov, f in zip(snakemake.params.coverages, snakemake.input))

report.to_csv(snakemake.output[0], sep="\t", index=False)
