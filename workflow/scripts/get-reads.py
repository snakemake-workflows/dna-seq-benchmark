import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam
import dnaio


def aln_to_fq(aln):
    return dnaio.SequenceRecord(
        name=aln.query_name,
        sequence=aln.get_forward_sequence(),
        quality=aln.query_qualities,
    )


limit = snakemake.params.limit
bam = pysam.AlignmentFile(snakemake.params.bam_url)

buffer = {}
n_written = 0
with dnaio.open(snakemake.output[0], snakemake.output[1], "w") as fqwriter:
    for aln in bam:
        if limit is not None and n_written >= limit:
            break
        if aln.is_secondary or aln.is_supplementary:
            continue

        qname = aln.query_name.removesuffix("/1").removesuffix("/2")

        mate_aln = buffer.get(aln.query_name)
        if mate_aln is None:
            buffer[aln.query_name] = aln
        else:
            if aln.is_read2:
                aln, mate_aln = mate_aln, aln

            fqwriter.write(aln_to_fq(aln), aln_to_fq(mate_aln))
            n_written += 1
