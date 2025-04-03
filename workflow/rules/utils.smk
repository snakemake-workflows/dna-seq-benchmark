rule index_vcf:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    log:
        "logs/bcftools-index-vcf/{prefix}.log",
    wrapper:
        "v1.9.0/bio/bcftools/index"


rule index_bcf:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi",
    log:
        "logs/bcftools-index-bcf/{prefix}.log",
    wrapper:
        "v1.9.0/bio/bcftools/index"
