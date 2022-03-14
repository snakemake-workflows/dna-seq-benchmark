import json
from pysam import VariantFile


with VariantFile(snakemake.input[0]) as infile, open(
    snakemake.output[0], "w"
) as outfile:
    json.dump({"isempty": not any(infile)}, outfile)
