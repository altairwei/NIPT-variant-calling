rule Module_1_Raw_FastQC:
    input:
        "data/{sample_id}.fq.gz"
    output:
        html=temp(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}", "{sample_id}_fastqc.html")),
        zip=temp(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}", "{sample_id}_fastqc.zip"))
    params:
        outdir=temp(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}"))
    log: get_log_path("{sample_id}")
    threads: 1
    shell:
        """
        fastqc -t {threads} \
            --outdir {params.outdir} {input} > {log} 2> {log}
        """

rule Module_1_QualityControl:
    input:
        "data/{sample_id}.fq.gz"
    output:
        fq=temp(os.path.join(OUTPUT_DIR, "clean", "{sample_id}.clean.fq.gz")),
        html=temp(os.path.join(OUTPUT_DIR, "clean", "{sample_id}.html")),
        json=temp(os.path.join(OUTPUT_DIR, "clean", "{sample_id}.json"))
    log:
        get_log_path("{sample_id}")
    shell:
        """
        fastp -i {input} -o {output.fq} \
            --qualified_quality_phred=5 \
            --unqualified_percent_limit=50 \
            --n_base_limit=10 \
            --adapter_sequence="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA" \
            --adapter_sequence_r2="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG" \
            --disable_trim_poly_g \
            --thread=16 \
            -j {output.json} \
            -h {output.html} \
            -R {wildcards.sample_id} 2> {log}
        """

rule Module_1_Raw_MultiQC:
    input:
        expand(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}",
            "{sample_id}_fastqc.html"), sample_id=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}",
            "{sample_id}_fastqc.zip"), sample_id=SAMPLES)
    output:
        protected(os.path.join(OUTPUT_DIR, "report", "raw.multiqc.html"))
    params:
        searchdirs=os.path.join(OUTPUT_DIR, "report", "raw"),
        outdir=os.path.join(OUTPUT_DIR, "report")
    log:
        get_log_path("multiqc_raw")
    shell:
        """
        multiqc \
            --force \
            -o {params.outdir} \
            -n raw.multiqc.html \
            {params.searchdirs} > {log} 2> {log}
        """

rule Module_1_Clean_MultiQC:
    input:
        expand(os.path.join(OUTPUT_DIR, "clean",
            "{sample_id}.json"), sample_id=SAMPLES)
    output:
        protected(os.path.join(OUTPUT_DIR, "report", "clean.multiqc.html"))
    params:
        searchdirs=os.path.join(OUTPUT_DIR, "clean"),
        outdir=os.path.join(OUTPUT_DIR, "report")
    log:
        get_log_path("multiqc_clean")
    shell:
        """
        multiqc \
            --force \
            -o {params.outdir} \
            -n clean.multiqc.html \
            {params.searchdirs} > {log} 2> {log}
        """
