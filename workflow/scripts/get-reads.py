import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam
import dnaio


def aln_to_fq(qname, aln):
    return dnaio.SequenceRecord(
        name=qname,
        sequence=aln.get_forward_sequence(),
        qualities="".join(
            map(lambda qual: chr(qual + 33), aln.get_forward_qualities())
        ),
    )


limit = snakemake.params.limit
bam = pysam.AlignmentFile(snakemake.params.bam_url)

buffer = {}
n_written = 0
with dnaio.open(snakemake.output[0], snakemake.output[1], mode="w") as fqwriter:
    for aln in bam:
        if limit is not None and n_written >= limit:
            break
        if aln.is_secondary or aln.is_supplementary:
            continue

        # Some aligners (e.g. minimap2) add /1 and /2 to the read name.
        # We remove them here to get the same name for both reads of a pair.
        qname = aln.query_name.removesuffix("/1").removesuffix("/2")

        mate_aln = buffer.get(qname)
        if mate_aln is None:
            buffer[qname] = aln
        else:
            if aln.is_read2:
                aln, mate_aln = mate_aln, aln
            del buffer[qname]

            fqwriter.write(aln_to_fq(qname, aln), aln_to_fq(qname, mate_aln))
            n_written += 1
