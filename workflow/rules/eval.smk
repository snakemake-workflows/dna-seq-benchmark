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


rule liftover_callset:
    input:
        callset=get_callset_correct_contigs,
        liftover_chain="resources/liftover/GRCh37_to_GRCh38.chain.gz",
        reference="resources/reference/genome.fasta",
        reference_dict="resources/reference/genome.fasta.dict",
    output:
        "results/normalized-variants/{callset}.lifted.vcf.gz",
    log:
        "logs/liftover_callset/{callset}.log",
    conda:
        "../envs/picard.yaml"
    resources:
        mem_mb=64000,
    shell:
        "picard LiftoverVcf -Xmx64g --MAX_RECORDS_IN_RAM 100000 -I {input.callset} -O {output} --CHAIN {input.liftover_chain} --REJECT {output}_rejected_variants.vcf -R {input.reference} &> {log}"


rule rename_contigs:
    input:
        calls=get_raw_callset,
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
        "resources/variants/{genome}/all.truth.norm.bcf",
    output:
        "resources/variants/{genome}/all.truth.format-added.vcf.gz",
    log:
        "logs/add_format_field/{genome}.log",
    conda:
        "../envs/vatools.yaml"
    shell:
        # TODO: Optional - Check first if FORMAT field is present for example with
        # TODO: bcftools view -h out.vcf.gz | grep FORMAT oder bcftools query -l all.bcf
        # bcftools convert makes sure that input for vcf-genotype-annotator is in vcf format
        # adds FORMAT field with GT field and sample name 'truth'
        "vcf-genotype-annotator <(bcftools convert -Ov {input}) truth 0/1 -o {output} &> {log}"


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
        "v3.3.6/bio/bcftools/view"


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
        "v1.7.2/bio/bcftools/index"


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
        "rm -r {params.output}; rtg vcfeval --threads {threads} --ref-overlap --all-records "
        "--output-mode ga4gh --baseline {input.truth} --calls {input.query} "
        "--output {params.output} --template {input.genome} {params.somatic} &> {log}"


rule calc_precision_recall:
    input:
        calls="results/vcfeval/{callset}/{cov}/output.vcf.gz",
        idx="results/vcfeval/{callset}/{cov}/output.vcf.gz.tbi",
        common_src=common_src,
        truth=get_stratified_truth(),
        truth_idx=get_stratified_truth(".tbi"),
        query="results/stratified-variants/{callset}/{cov}.vcf.gz",
        query_index="results/stratified-variants/{callset}/{cov}.vcf.gz.tbi",
    output:
        snvs="results/precision-recall/callsets/{callset}/{cov}.{vartype}.tsv",
    log:
        "logs/calc-precision-recall/{callset}/{cov}/{vartype}.log",
    params:
        vaf_fields=get_vaf_fields,
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/calc-precision-recall.py"


rule collect_stratifications:
    input:
        get_collect_stratifications_input,
    output:
        "results/precision-recall/callsets/{callset}.{vartype}.tsv",
    params:
        coverages=get_nonempty_coverages,
        coverage_lower_bounds=get_coverages,
    log:
        "logs/collect-stratifications/{callset}/{vartype}.log",
    conda:
        "../envs/stats.yaml"
    # We want this to be determined before FP/FN collection in order to avoid memory
    # issues with callsets that do not match the truth at all.
    priority: 2
    script:
        "../scripts/collect-stratifications.py"


rule collect_precision_recall:
    input:
        tables=get_collect_precision_recall_input,
    output:
        "results/precision-recall/benchmarks/{benchmark}.{vartype}.tsv",
    params:
        callsets=lambda w: get_benchmark_callsets(w.benchmark),
        labels=get_collect_precision_recall_labels,
        vaf=get_vaf_status,
    log:
        "logs/collect-precision-recall/{benchmark}/{vartype}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/collect-precision-recall.py"


rule report_precision_recall:
    input:
        config=workflow.source_path(
            "../resources/datavzrd/precision-recall-config.yte.yaml"
        ),
        table="results/precision-recall/benchmarks/{benchmark}.{vartype}.tsv",
    output:
        report(
            directory("results/report/precision-recall/{benchmark}/{vartype}"),
            htmlindex="index.html",
            category="precision/recall",
            labels={
                "benchmark": "{benchmark}",
                "vartype": "{vartype}",
            },
        ),
    log:
        "logs/datavzrd/precision-recall/{benchmark}/{vartype}.log",
    params:
        somatic=get_somatic_status,
        vaf=get_vaf_status,
        high_coverage=get_high_coverage_status,
        genome=get_genome_name,
        version=get_genome_version,
    wrapper:
        "v3.13.7/utils/datavzrd"


rule extract_fp_fn:
    input:
        calls="results/vcfeval/{callset}/{cov}/output.vcf.gz",
        common_src=common_src,
    output:
        "results/fp-fn/callsets/{cov}/{callset}/{classification}.tsv",
    log:
        "logs/extract-fp-fn/{cov}/{callset}/{classification}.log",
    conda:
        "../envs/vembrane.yaml"
    script:
        "../scripts/extract-fp-fn.py"


# TODO: if one of the input callsets has all sites as FN, the resulting merged table
# becomes empty. This needs to be fixed, but is not super urgent because it is not realistic to
# happen in reality.
rule collect_fp_fn:
    input:
        tables=get_collect_fp_fn_input,
    output:
        main="results/fp-fn/genomes/{genome}/{cov}/{classification}/main.tsv",
        dependency_sorting=directory(
            "results/fp-fn/genomes/{genome}/{cov}/{classification}/dependency-sorting"
        ),
    params:
        callsets=get_collect_fp_fn_callsets,
        labels=get_collect_fp_fn_labels,
        label_names=lambda w: get_callsets_labels(get_genome_callsets(w.genome)),
        max_entries=config.get("max-fp-fn-entries", 300),
    log:
        "logs/collect-fp-fn/{genome}/{cov}/{classification}.log",
    conda:
        "../envs/stats.yaml"
    # This has to happen after precision/recall has been computed, otherwise we risk
    # extremely high memory usage if a callset does not match the truth at all.
    priority: 1
    script:
        "../scripts/collect-fp-fn.py"


rule report_fp_fn:
    input:
        main_dataset="results/fp-fn/genomes/{genome}/{cov}/{classification}/main.tsv",
        dependency_sorting_datasets="results/fp-fn/genomes/{genome}/{cov}/{classification}/dependency-sorting",
        config=workflow.source_path("../resources/datavzrd/fp-fn-config.yte.yaml"),
    output:
        report(
            directory("results/report/fp-fn/{genome}/{cov}/{classification}"),
            htmlindex="index.html",
            category="{classification} variants",
            subcategory=lambda w: w.genome,
            labels=lambda w: {"coverage": w.cov},
        ),
    log:
        "logs/datavzrd/fp-fn/{genome}/{cov}/{classification}.log",
    params:
        labels=lambda w: get_callsets_labels(get_genome_callsets(w.genome)),
        version=get_genome_version,
    wrapper:
        "v3.13.7/utils/datavzrd"
