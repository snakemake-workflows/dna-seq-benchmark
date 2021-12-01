rule get_reads:
    output:
        r1=reads[0],
        r2=reads[1],
    log:
        "logs/download-reads.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(samtools view -f3 -u"
        " ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam |"
        " samtools sort -n -u | samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -) 2> {log}"


rule get_truth:
    output:
        "benchmark/truth.vcf",
    log:
        "logs/get-truth.log",
    params:
        repl_chr=repl_chr,
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools view "
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "
        " | sed {params.repl_chr} > {output} 2> {log}"


rule get_confidence_bed:
    output:
        "benchmark/confidence-regions.bed",
    log:
        "logs/get-confidence-regions.log",
    params:
        repl_chr=repl_chr,
    conda:
        "../envs/tools.yaml"
    shell:
        "curl --insecure -L "
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed | "
        "sed {params.repl_chr} > {output} 2> {log}"


rule get_liftover_track:
    output:
        "benchmark/liftover.chain.gz",
    log:
        "logs/get-liftover-track.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "curl --insecure -L http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz > {output} 2> {log}"


rule get_target_bed:
    input:
        liftover="benchmark/liftover.chain.gz",
    output:
        pipe("benchmark/target-regions.raw.bed"),
    log:
        "logs/get-target-bed.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(curl --insecure -L"
        " ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/TruSeq_exome_targeted_regions.hg19.bed |"
        " liftOver /dev/stdin {input.liftover} {output} /dev/null) 2> {log}"


rule postprocess_target_bed:
    input:
        "benchmark/target-regions.raw.bed",
    output:
        "benchmark/target-regions.bed",
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
        "reference/reference.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="104",
    log:
        "logs/get-genome.log",
    wrapper:
        "0.79.0/bio/reference/ensembl-sequence"


rule samtools_faidx:
    input:
        "reference/reference.fasta",
    output:
        "reference/reference.fasta.fai",
    log:
        "logs/samtools-faidx.log",
    wrapper:
        "0.79.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "reference/reference.fasta",
    output:
        multiext("reference/reference", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa-index.log",
    params:
        prefix=get_io_prefix(lambda input, output: output[0]),
    wrapper:
        "0.79.0/bio/bwa/index"


rule bwa_mem:
    input:
        reads=reads,
        index=multiext("reference/reference", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "reads/mapped.bam",
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
        "reads/mapped.bam",
    output:
        bam="reads/mapped.dedup.bam",
        metrics="reads/dedup.metrics.txt",
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
        "reads/mapped.dedup.bam",
    output:
        "reads/mapped.dedup.bam.bai",
    log:
        "logs/samtools-index.log",
    wrapper:
        "0.79.0/bio/samtools/index"


rule mosdepth:
    input:
        bam="reads/mapped.dedup.bam",
        bai="reads/mapped.dedup.bam.bai",
    output:
        "coverage/coverage.mosdepth.global.dist.txt",
        "coverage/coverage.quantized.bed.gz",
        summary="coverage/coverage.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth.log",
    params:
        extra="--no-per-base --mapq 59",  # we do not want low MAPQ regions end up being marked as high coverage
        quantize=get_mosdepth_quantize(),
    wrapper:
        "0.77.0/bio/mosdepth"


rule stratify_regions:
    input:
        confidence="benchmark/confidence-regions.bed",
        target="benchmark/target-regions.bed",
        coverage="coverage/coverage.quantized.bed.gz",
    output:
        "benchmark/test-regions.cov-{cov}.bed",
    log:
        "logs/stratify-regions/{cov}.log",
    params:
        cov_label=get_cov_label,
    conda:
        "../envs/tools.yaml"
    shell:
        "(bedtools intersect"
        " -a {input.confidence}"
        " -b <(zcat {input.coverage} | grep '{params.cov_label}') |"
        " bedtools intersect -a /dev/stdin -b {input.target}"
        ") > {output} 2> {log}"
