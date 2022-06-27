import os
from subprocess import check_output
import subprocess
import pysam

outdir = snakemake.output[0]
os.makedirs(outdir, exist_ok=True)

with pysam.VariantFile(snakemake.input.vcf[0]) as invcf:
    for record in invcf:
        bd = record.samples["TRUTH"]["BD"][0]
        if bd == "FN":
            region = f"{record.contig}:{record.pos - 100}-{record.pos + 100}"
            highlight = f"{record.pos}-{record.pos + len(record.alleles[1])}"
            subprocess.run(
                f"alignoth -h {highlight} -r {snakemake.input.genome} -g {region} | "
                f"vl2svg > {outdir}/{record.contig}:{record.pos}:{':'.join(record.alleles)}.svg",
                check=True,
                shell=True,
            )
