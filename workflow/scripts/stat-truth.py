import json
from pysam import VariantFile


with VariantFile(open(snakemake.input[0], "rb")) as infile, open(
    snakemake.output[0], "w"
) as outfile:
    json.dump({"isempty": not any(infile)}, outfile)
