if not config["custom-reads"]["activate"]:

    rule get_reads:
        input:
            regions=get_limit_regions(),
        output:
            r1=public_reads[0],
            r2=public_reads[1],
        params:
            limit=get_read_limit_param,
        log:
            "logs/download-reads.log",
        conda:
            "../envs/tools.yaml"
        shell:
            "(samtools view -f3 -u"
            " ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam"
            " {params.limit} |"
            " samtools sort -n -u | samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -) 2> {log}"


rule get_truth:
    output:
        "resources/variants/truth.vcf",
    log:
        "logs/get-truth.log",
    params:
        repl_chr=repl_chr,
        url=get_truth_url(),
    conda:
        "../envs/tools.yaml"
    shell:
        "(bcftools view {params.url}"
        " | sed {params.repl_chr} > {output}"
        ") 2> {log}"


rule get_confidence_bed:
    output:
        "resources/regions/confidence-regions.bed",
    log:
        "logs/get-confidence-regions.log",
    params:
        repl_chr=repl_chr,
        url=get_confidence_bed_url(),
    conda:
        "../envs/tools.yaml"
    shell:
        "(curl --insecure -L {params.url}"
        " | sed {params.repl_chr} > {output}"
        ") 2> {log}"


rule get_liftover_track:
    output:
        "resources/reference/liftover.chain.gz",
    log:
        "logs/get-liftover-track.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "curl --insecure -L http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz > {output} 2> {log}"


rule get_target_bed:
    input:
        liftover="resources/reference/liftover.chain.gz",
    output:
        "resources/regions/target-regions.raw.bed",
    params:
        get_bed=get_target_bed_statement(),
        liftover=get_liftover_statement,
    log:
        "logs/get-target-bed.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "({params.get_bed} {params.liftover}) 2> {log}"


rule postprocess_target_bed:
    input:
        "resources/regions/target-regions.raw.bed",
    output:
        "resources/regions/target-regions.bed",
    log:
        "logs/fix-target-bed.log",
    params:
        repl_chr=repl_chr,
    conda:
        "../envs/tools.yaml"
    shell:
        "sed {params.repl_chr} {input} > {output} 2> {log}"


rule get_reference:
    output:
        "resources/reference/genome.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh37" if config["grch37"] else "GRCh38",
        release="104",
    log:
        "logs/get-genome.log",
    wrapper:
        "0.79.0/bio/reference/ensembl-sequence"


rule samtools_faidx:
    input:
        "resources/reference/genome.fasta",
    output:
        "resources/reference/genome.fasta.fai",
    log:
        "logs/samtools-faidx.log",
    wrapper:
        "0.79.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/reference/genome.fasta",
    output:
        multiext("reference/reference/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa-index.log",
    params:
        prefix=get_io_prefix(lambda input, output: output[0]),
    wrapper:
        "0.79.0/bio/bwa/index"


rule bwa_mem:
    input:
        reads=get_bwa_input(),
        index=rules.bwa_index.output,
    output:
        "results/read-alignments/all.bam",
    log:
        "logs/bwa-mem.log",
    params:
        index=get_io_prefix(lambda input, output: input.index[0]),
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
    threads: 8
    wrapper:
        "0.79.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/read-alignments/all.bam",
    output:
        bam="results/read-alignments/all.dedup.bam",
        metrics="results/read-alignments/dedup.metrics.txt",
    log:
        "logs/picard-dedup.log",
    params:
        extra="REMOVE_DUPLICATES=true",
    resources:
        mem_mb=1024,
    wrapper:
        "0.79.0/bio/picard/markduplicates"


rule samtools_index:
    input:
        "results/read-alignments/all.dedup.bam",
    output:
        "results/read-alignments/all.dedup.bam.bai",
    log:
        "logs/samtools-index.log",
    wrapper:
        "0.79.0/bio/samtools/index"


rule mosdepth:
    input:
        bam="results/read-alignments/all.dedup.bam",
        bai="results/read-alignments/all.dedup.bam.bai",
    output:
        "results/coverage/coverage.mosdepth.global.dist.txt",
        "results/coverage/coverage.quantized.bed.gz",
        summary="results/coverage/coverage.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth.log",
    params:
        extra="--no-per-base --mapq 59",  # we do not want low MAPQ regions end up being marked as high coverage
        quantize=get_mosdepth_quantize(),
    wrapper:
        "0.77.0/bio/mosdepth"


rule stratify_regions:
    input:
        confidence="resources/regions/confidence-regions.bed",
        target="resources/regions/target-regions.bed",
        coverage="results/coverage/coverage.quantized.bed.gz",
        limit_regions=get_limit_regions(),
    output:
        "resources/regions/test-regions.cov-{cov}.bed",
    params:
        intersect_limit=get_limit_regions_intersect_statement,
        cov_label=get_cov_label,
    log:
        "logs/stratify-regions/{cov}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        # The intersection can generate duplicate or overlapping entries if
        # the target regions bed contains overlapping entries.
        # Hence it is important to have a final bedtools merge to ensure that
        # these overlapping lines are removed.
        "(bedtools intersect"
        " -a {input.confidence}"
        " -b <(zcat {input.coverage} | grep '{params.cov_label}') |"
        " bedtools intersect -a /dev/stdin -b {input.target}"
        " {params.intersect_limit} |"
        " sort -k1,1 -k2,2n |"
        " bedtools merge -i /dev/stdin"
        ") > {output} 2> {log}"
