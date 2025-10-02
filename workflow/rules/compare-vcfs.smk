rule get_reference_dict:
    input:
        reference="resources/reference/genome.fasta",
    output:
        "resources/reference/genome.fasta.dict",
    log:
        "logs/get-reference-dict.log",
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CreateSequenceDictionary  -R {input.reference} -O {input.reference}.dict &> {log}"


rule merge_callsets:
    input:
        snv_vcf=lambda wildcards: get_raw_callset(wildcards)["snvs"],
        indel_vcf=lambda wildcards: get_raw_callset(wildcards)["indels"],
        snv_tbi=lambda wildcards: get_raw_callset_index(wildcards)["snvs"],
        indel_tbi=lambda wildcards: get_raw_callset_index(wildcards)["indels"],
    output:
        "results/merge-callsets/{callset}.merged.vcf.gz",
    log:
        "logs/merge-callsets/{callset}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools concat -O z --allow-overlap {input.snv_vcf} {input.indel_vcf} > {output} 2> {log}"


rule liftover_callset:
    input:
        callset=get_callset_correct_contigs,
        liftover_chain="resources/liftover/GRCh37_to_GRCh38.chain.gz",
        reference="resources/reference/genome.fasta",
        reference_dict="resources/reference/genome.fasta.dict",
    output:
        lifted="results/normalized-variants/{callset}.lifted.vcf.gz",
        rejected="results/normalized-variants/{callset}.lifted_rejected_variants.vcf.gz",
    log:
        "logs/liftover_callset/{callset}.log",
    conda:
        "../envs/picard.yaml"
    resources:
        mem_mb=64000,
    shell:
        "picard LiftoverVcf -Xmx64g --MAX_RECORDS_IN_RAM 100000 \
         -I {input.callset} \
         -O {output.lifted} \
         --CHAIN {input.liftover_chain} \
         --REJECT {output.rejected} \
         -R {input.reference} &> {log}"


rule rename_contigs:
    input:
        calls=get_callset_correct_contigs_liftover_merge,
        repl_file=get_rename_contig_file,
    output:
        "results/normalized-variants/{callset}.replaced-contigs.vcf.gz",
    log:
        "logs/rename-contigs/{callset}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools annotate {input.calls} --rename-chrs {input.repl_file} "
        "-Oz -o {output} 2> {log}"


rule add_genotype_field:
    input:
        get_callset_correct_contigs_liftover,
    output:
        "results/normalized-variants/{callset}.gt-added.vcf.gz",
    log:
        "logs/add_genotype_field/{callset}.log",
    params:
        get_somatic_sample_name,
    conda:
        "../envs/vatools.yaml"
    shell:
        # part after || gets executed if vcf-genotype-annotater fails because GT field is already present
        # bcftools convert makes sure that input for vcf-genotype-annotator is in vcf format
        "vcf-genotype-annotator <(bcftools convert -Ov {input}) {params} 0/1 -o {output} &> {log} || bcftools view {input} -Oz > {output}"


rule add_format_field:
    input:
        bcf="resources/variants/{genome}/all.truth.norm.bcf",
    output:
        "resources/variants/{genome}/all.truth.format-added.vcf.gz",
    log:
        "logs/add_format_field/{genome}.log",
    conda:
        "../envs/vatools.yaml"
    shell:
        """
        if bcftools view -h {input.bcf} | grep -q FORMAT; then
            bcftools reheader -s <(echo 'truth') {input.bcf} | bcftools view -Oz > {output}
        else
            vcf-genotype-annotator <(bcftools convert -Ov {input.bcf}) truth 0/1 -o {output} &> {log}
        fi
        """


rule remove_non_pass:
    input:
        get_callset,
    output:
        "results/filtered-variants/{callset}.bcf",
    log:
        "logs/filter/{callset}.log",
    params:
        extra="-f 'PASS,.'",
    wrapper:
        "v7.6.0/bio/bcftools/view"


rule intersect_calls_with_target_regions:
    input:
        bcf="results/filtered-variants/{callset}.bcf",
        regions=get_target_regions,
    output:
        pipe("results/normalized-variants/{callset}_intersected.vcf"),
    log:
        "logs/intersect-calls/{callset}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bedtools intersect -b {input.regions} -a "
        "<(bcftools view {input.bcf}) -wa -f 1.0 -header > {output}) 2> {log}"


rule restrict_to_reference_contigs:
    input:
        calls="results/filtered-variants/{callset}.bcf",
        calls_index="results/filtered-variants/{callset}.bcf.csi",
        ref_index="resources/reference/genome.fasta.fai",
    output:
        "results/filtered-variants/{callset}_restricted.bcf",
    log:
        "logs/restrict-to-reference-contigs/{callset}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bcftools view --regions $(cut -f1 {input.ref_index} | tr '\\n' ',') {input.calls} |"
        " bcftools reheader -f {input.ref_index} > {output}) 2> {log}"


rule normalize_calls:
    input:
        calls=branch(
            intersect_calls,
            then="results/normalized-variants/{callset}_intersected.vcf",
            otherwise="results/filtered-variants/{callset}_restricted.bcf",
        ),
        ref="resources/reference/genome.fasta",
        ref_index="resources/reference/genome.fasta.fai",
    output:
        "results/normalized-variants/{callset}.vcf.gz",
    params:
        extra=get_norm_params,
    log:
        "logs/normalize-calls/{callset}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bcftools norm {params.extra} --fasta-ref {input.ref} {input.calls} | "
        "bcftools view -Oz > {output}) 2> {log}"


rule stratify_truth:
    input:
        variants=get_benchmark_truth,
        regions="resources/regions/{benchmark}/test-regions.cov-{cov}.bed",
    output:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz",
    log:
        "logs/stratify-truth/{benchmark}.{cov}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bedtools intersect -b {input.regions} -a "
        "<(bcftools view {input.variants} | bcftools reheader -s <(echo 'truth')) -wa -f 1.0 -header | "
        "bcftools view -Oz > {output}) 2> {log}"


rule stratify_results:
    input:
        variants="results/normalized-variants/{callset}.vcf.gz",
        regions=get_test_regions,
    output:
        "results/stratified-variants/{callset}/{cov}.vcf.gz",
    log:
        "logs/stratify-results/{callset}/{cov}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bedtools intersect -b {input.regions} -a "
        "<(bcftools view {input.variants}) -wa -f 1.0 -header | "
        "bcftools view -Oz > {output}) 2> {log}"


rule index_stratified_truth:
    input:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz",
    output:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz.csi",
    log:
        "logs/bcftools-index/{benchmark}.truth.{cov}.log",
    wrapper:
        "v7.6.0/bio/bcftools/index"


checkpoint stat_truth:
    input:
        "results/variants/{benchmark}.truth.cov-{cov}.vcf.gz",
    output:
        "results/variants/{benchmark}.truth.cov-{cov}.stats.json",
    log:
        "logs/stat-truth/{benchmark}.{cov}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/stat-truth.py"


rule generate_sdf:
    input:
        genome="resources/reference/genome.fasta",
        genome_index="resources/reference/genome.fasta.fai",
    output:
        directory("resources/reference/genome-sdf"),
    log:
        "logs/rtg-tools/sdf.log",
    conda:
        "../envs/rtg-tools.yaml"
    shell:
        "rtg format --output {output} {input.genome} &> {log}"


rule benchmark_variants:
    input:
        truth=get_stratified_truth(),
        truth_idx=get_stratified_truth(".tbi"),
        query="results/stratified-variants/{callset}/{cov}.vcf.gz",
        query_index="results/stratified-variants/{callset}/{cov}.vcf.gz.tbi",
        genome="resources/reference/genome-sdf",
    output:
        "results/vcfeval/{callset}/{cov}/output.vcf.gz",
    log:
        "logs/vcfeval/{callset}/{cov}.log",
    params:
        output=lambda w, output: os.path.dirname(output[0]),
        somatic=get_somatic_flag,
    conda:
        "../envs/rtg-tools.yaml"
    threads: 32
    shell:
        "rm -r {params.output}; rtg vcfeval --threads {threads} --ref-overlap --all-records --no-roc "
        "--output-mode ga4gh --baseline {input.truth} --calls {input.query} "
        "--output {params.output} --template {input.genome} {params.somatic} &> {log}"
