rule index_variants:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi",
    log:
        "logs/index-variants/{prefix}.log",
    wrapper:
        "v1.7.2/bio/bcftools/index"
