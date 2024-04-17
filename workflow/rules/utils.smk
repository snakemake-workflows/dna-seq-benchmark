rule index_vcf_variants:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    log:
        "logs/index-vcf-variants/{prefix}.log",
    wrapper:
        "v1.9.0/bio/bcftools/index"
