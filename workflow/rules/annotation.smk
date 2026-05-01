rule get_downsampled_vep_cache:
    input:
        workflow.source_path(
            "../resources/ci-test-references/vep_cache_113_GRCh38_chr22.tar.gz",
        ),
    output:
        directory("resources/vep/cache_downsampled"),
    log:
        "logs/vep/downsampled_cache.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(mkdir -p {output}; tar -xzf {input} -C {output} --strip-components 1) 2> {log}"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species="homo_sapiens",
        build=get_reference_genome_build(),
        release="115",
    log:
        "logs/vep/cache.log",
    cache: "omit-software"
    wrapper:
        "v9.4.2/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release="115",
    log:
        "logs/vep/plugins.log",
    wrapper:
        "v9.5.0/bio/vep/plugins"


rule download_revel:
    output:
        temp("resources/revel_scores.zip"),
    log:
        "logs/vep_plugins/download_revel.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "curl https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip -o {output} &> {log}"


rule process_revel_scores:
    input:
        "resources/revel_scores.zip",
    output:
        "resources/revel_scores.tsv.gz",
    params:
        build=get_reference_genome_build(),
    log:
        "logs/vep_plugins/process_revel_scores.log",
    resources:
        tmpdir=temp("tmpdir"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        tmpfile=$(mktemp "${{TMPDIR:-/tmp}}"/revel_scores.XXXXXX)
        trap "rm -f $tmpfile" EXIT
        unzip -p {input} | tr "," "\t" | sed '1s/.*/#&/' | bgzip -c > $tmpfile
        if [ "{params.build}" == "GRCh38" ] ; then
            zgrep -h -v ^#chr $tmpfile | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat <(zcat $tmpfile | head -n1) - | bgzip -c > {output}
        elif [ "{params.build}" == "GRCh37" ] ; then
            cat $tmpfile > {output}
        else
            echo "Annotation of REVEL scores only supported for GRCh37 or GRCh38" > {log}
            exit 125
        fi
        """


rule tabix_revel_scores:
    input:
        "resources/revel_scores.tsv.gz",
    output:
        "resources/revel_scores.tsv.gz.tbi",
    params:
        get_tabix_revel_params(),
    log:
        "logs/tabix/revel.log",
    wrapper:
        "v9.4.1/bio/tabix/index"


rule annotate_shared_fn:
    input:
        calls="results/fp-fn/vcf/{benchmark}/{benchmark}.shared_fn.sorted.vcf.gz",
        cache=get_vep_cache_dir(),
        plugins=access.random("resources/vep/plugins"),
        revel=lambda wc: get_plugin_aux("REVEL"),
        revel_tbi=lambda wc: get_plugin_aux("REVEL", True),
        fasta=access.random("resources/reference/genome.fasta"),
        fai="resources/reference/genome.fasta.fai",
    output:
        calls="results/annotated/vcf/{benchmark}/{benchmark}.shared_fn.annotated.vcf.gz",
        stats="results/annotated/stats/{benchmark}/{benchmark}.shared_fn.stats.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["REVEL"],
        extra="--everything --check_existing --vcf_info_field ANN --hgvsg --sift b --polyphen b ",
    log:
        "logs/vep/fp-fn/{benchmark}/{benchmark}.shared_fn.annotate.log",
    threads: 4
    group:
        "annotation"
    wrapper:
        "v9.4.0/bio/vep/annotate"


rule annotate_unique_fp_fn:
    input:
        calls="results/fp-fn/vcf/{benchmark}/{callset}.unique_{classification}.sorted.vcf.gz",
        cache=get_vep_cache_dir(),
        plugins=access.random("resources/vep/plugins"),
        revel=lambda wc: get_plugin_aux("REVEL"),
        revel_tbi=lambda wc: get_plugin_aux("REVEL", True),
        fasta=access.random("resources/reference/genome.fasta"),
        fai="resources/reference/genome.fasta.fai",
    output:
        calls="results/annotated/vcf/{benchmark}/{callset}.unique_{classification}.annotated.vcf.gz",
        stats="results/annotated/stats/{benchmark}/{callset}.unique_{classification}.stats.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["REVEL"],
        extra="--everything --check_existing --vcf_info_field ANN --hgvsg --sift b --polyphen b ",
    log:
        "logs/vep/fp-fn/{benchmark}/{callset}.unique_{classification}.annotate.log",
    threads: 4
    group:
        "annotation"
    wrapper:
        "v9.4.0/bio/vep/annotate"


rule vembrane_table_shared_fn:
    input:
        "results/annotated/vcf/{benchmark}/{benchmark}.shared_fn.annotated.norm.vcf",
    output:
        "results/annotated/tsv/{benchmark}/{benchmark}.shared_fn.annotated.tsv",
    params:
        expression='CHROM, POS, REF, ALT, ANN["SYMBOL"], ANN["VARIANT_CLASS"], ANN["IMPACT"], ANN["Consequence"], (lambda d: next(iter(d.keys())) if d else "")(ANN["SIFT"]), (lambda d: next(iter(d.values())) if d else None)(ANN["SIFT"]), (lambda d: next(iter(d.keys())) if d else "")(ANN["PolyPhen"]), (lambda d: next(iter(d.values())) if d else None)(ANN["PolyPhen"]), ANN["REVEL"]',
        extra="--header 'CHROM,POS,REF,ALT,SYMBOL,VARIANT_CLASS,IMPACT,Consequence,SIFT,SIFT_SCORE,PolyPhen,PolyPhen_SCORE,REVEL'",
    log:
        "logs/vembrane/{benchmark}/{benchmark}.shared_fn.annotate.log",
    wrapper:
        "v9.4.0/bio/vembrane/table"


rule vembrane_table_unique_fp_fn:
    input:
        "results/annotated/vcf/{benchmark}/{callset}.unique_{classification}.annotated.norm.vcf",
    output:
        "results/annotated/tsv/{benchmark}/{callset}.unique_{classification}.annotated.tsv",
    params:
        expression='CHROM, POS, REF, ALT, ANN["SYMBOL"], ANN["VARIANT_CLASS"], ANN["IMPACT"], ANN["Consequence"], (lambda d: next(iter(d.keys())) if d else "")(ANN["SIFT"]), (lambda d: next(iter(d.values())) if d else None)(ANN["SIFT"]), (lambda d: next(iter(d.keys())) if d else "")(ANN["PolyPhen"]), (lambda d: next(iter(d.values())) if d else None)(ANN["PolyPhen"]), ANN["REVEL"]',
        extra="--header 'CHROM,POS,REF,ALT,SYMBOL,VARIANT_CLASS,IMPACT,Consequence,SIFT,SIFT_SCORE,PolyPhen,PolyPhen_SCORE,REVEL'",
    log:
        "logs/vembrane/{benchmark}/{callset}.unique_{classification}.annotate.log",
    wrapper:
        "v9.4.0/bio/vembrane/table"
