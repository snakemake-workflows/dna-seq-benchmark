rule get_reads:
    output:
        r1="resources/reads/{benchmark}.1.fq",
        r2="resources/reads/{benchmark}.2.fq",
    params:
        limit=branch(
            lookup("limit-reads", within=config, default=False),
            then=100000,
            otherwise=None,
        ),
        bam_url=get_benchmark_bam_url,
    log:
        "logs/download-reads/{benchmark}.log",
    conda:
        "../envs/pysam.yaml"
    retries: 3
    script:
        "../scripts/get-reads.py"


rule get_archive:
    output:
        directory("resources/archives/{genome}"),
    params:
        url=get_archive_url,
    log:
        "logs/get-archive/{genome}.log",
    conda:
        "../envs/tools.yaml"
    retries: 3
    shell:
        "(mkdir -p {output}; curl -L {params.url} | tar -x -C {output} --strip-components 1) 2> {log}"


rule get_truth:
    input:
        archive=get_archive_input,
    output:
        "resources/variants/{genome}/{truthset}.truth.bcf",
    log:
        "logs/get-truth/{genome}/{truthset}.log",
    params:
        repl_chr=repl_chr,
        url=get_truth_url,
    conda:
        "../envs/tools.yaml"
    shell:
        "ls {input.archive}; (bcftools view {params.url}"
        " | sed {params.repl_chr} | bcftools view -Ob - > {output}"
        ") 2> {log}"


rule merge_truthsets:
    input:
        bcf=get_truthsets(),
        csi=get_truthsets(csi=True),
    output:
        "resources/variants/{genome}.merged.truth.bcf",
    log:
        "logs/merge-truthsets/{genome}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools concat -O b --allow-overlap {input.bcf} > {output} 2> {log}"


rule normalize_truth:
    input:
        get_genome_truth,
        ref="resources/reference/genome.fasta",
        ref_index="resources/reference/genome.fasta.fai",
    output:
        "resources/variants/{genome}/all.truth.norm.bcf",
    params:
        extra=get_norm_params,
    log:
        "logs/normalize-truth/{genome}.log",
    wrapper:
        "v1.9.0/bio/bcftools/norm"


rule get_confidence_bed:
    input:
        archive=get_archive_input,
    output:
        "resources/regions/{genome}.confidence-regions.bed",
    log:
        "logs/get-confidence-regions/{genome}.log",
    params:
        repl_chr=repl_chr,
        cmd=get_confidence_bed_cmd,
    conda:
        "../envs/tools.yaml"
    shell:
        "({params.cmd} | sed {params.repl_chr} > {output}) 2> {log}"


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
        target_bed=get_target_bed_input,
        liftover="resources/reference/liftover.chain.gz",
    output:
        "resources/regions/{benchmark}/target-regions.raw.bed",
    params:
        get_bed=get_target_bed_statement,
        liftover=get_liftover_statement,
    log:
        "logs/get-target-bed/{benchmark}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "({params.get_bed} {params.liftover}) 2> {log}"


rule postprocess_target_bed:
    input:
        "resources/regions/{benchmark}/target-regions.raw.bed",
    output:
        "resources/regions/{benchmark}/target-regions.bed",
    log:
        "logs/fix-target-bed/{benchmark}.log",
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
        build=get_reference_genome_build(),
        release="104",
        chromosome="1" if config.get("limit-reads") else None,
    log:
        "logs/get-genome.log",
    wrapper:
        "v1.7.2/bio/reference/ensembl-sequence"


rule get_liftover_chain:
    output:
        "resources/liftover/GRCh37_to_GRCh38.chain.gz",
    log:
        "logs/get_liftover_chain.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "curl http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz -o {output}"


rule samtools_faidx:
    input:
        "resources/reference/genome.fasta",
    output:
        "resources/reference/genome.fasta.fai",
    log:
        "logs/samtools-faidx.log",
    wrapper:
        "v6.2.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/reference/genome.fasta",
    output:
        idx=multiext(
            "resources/reference/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "logs/bwa_index/genome.log",
    wrapper:
        "v7.0.0/bio/bwa/index"


rule bwa_mem:
    input:
        reads=get_bwa_input,
        idx=rules.bwa_index.output,
    output:
        "results/read-alignments/{benchmark}.bam",
    log:
        "logs/bwa-mem/{benchmark}.log",
    params:
        extra=r"-R '@RG\tID:{benchmark}\tSM:{benchmark}'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
    threads: 8
    wrapper:
        "v7.0.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        bams="results/read-alignments/{benchmark}.bam",
    output:
        bam="results/read-alignments/{benchmark}.dedup.bam",
        metrics="results/read-alignments/{benchmark}.dedup.metrics.txt",
    log:
        "logs/picard-dedup/{benchmark}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v7.0.0/bio/picard/markduplicates"


rule samtools_index:
    input:
        "results/read-alignments/{benchmark}.dedup.bam",
    output:
        "results/read-alignments/{benchmark}.dedup.bam.bai",
    log:
        "logs/samtools-index/{benchmark}.log",
    wrapper:
        "v1.7.2/bio/samtools/index"


rule mosdepth:
    input:
        bam=get_mosdepth_input(),
        bai=get_mosdepth_input(bai=True),
    output:
        "results/coverage/{benchmark}/coverage.mosdepth.global.dist.txt",
        "results/coverage/{benchmark}/coverage.quantized.bed.gz",
        summary="results/coverage/{benchmark}/coverage.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/{benchmark}.log",
    params:
        extra="--no-per-base --mapq 59",  # we do not want low MAPQ regions end up being marked as high coverage
        quantize=get_mosdepth_quantize,
    wrapper:
        "v1.7.2/bio/mosdepth"


rule stratify_regions:
    input:
        confidence=get_confidence_regions,
        target=get_target_regions,
        coverage="results/coverage/{benchmark}/coverage.quantized.bed.gz",
    output:
        "resources/regions/{benchmark}/test-regions.cov-{cov}.bed",
    params:
        cov_label=get_cov_label,
        intersect_target_regions=get_target_regions_intersect_statement,
    log:
        "logs/stratify-regions/{benchmark}/{cov}.log",
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
        " {params.intersect_target_regions}"
        " sort -k1,1 -k2,2n |"
        " bedtools merge -i /dev/stdin "
        ") > {output} 2> {log}"
