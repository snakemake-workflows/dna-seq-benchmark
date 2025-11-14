rule index_vcf:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    log:
        "logs/bcftools-index-vcf/{prefix}.log",
    wrapper:
        "v8.0.2/bio/bcftools/index"


rule index_bcf:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi",
    log:
        "logs/bcftools-index-bcf/{prefix}.log",
    wrapper:
        "v8.0.2/bio/bcftools/index"


rule sort_vcf:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.sorted.vcf.gz",
    log:
        "logs/bcftools-sort-vcf/{prefix}.log",
    wrapper:
        "v8.0.2/bio/bcftools/sort"
