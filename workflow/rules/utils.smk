rule index_variants:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi",
    log:
        "logs/index-variants/{prefix}.log",
    wrapper:
        "v1.2.0/bio/bcftools/index"


rule index_vcf_variants:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    params:
        extra="--tbi",
    log:
        "logs/index-vcf-variants/{prefix}.log",
    wrapper:
        "v1.2.0/bio/bcftools/index"
