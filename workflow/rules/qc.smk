rule Module_1_Raw_FastQC:
    input:
        config["file_pattern"]
    output:
        html=temp(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}",
            config["name_pattern"] + "_fastqc.html")),
        zip=temp(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}",
            config["name_pattern"] + "_fastqc.zip"))
    params:
        outdir=temp(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}"))
    log: get_log_path("{sample_id}")
    threads: 1
    resources:
        mem_mb=2048
    shell:
        """
        fastqc -t {threads} \
            --memory {resources.mem_mb} \
            --outdir {params.outdir} {input} > {log} 2> {log}
        """

rule Module_1_Raw_MultiQC:
    input:
        expand(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}",
            config["name_pattern"] + "_fastqc.html"), sample_id=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "report", "raw", "{sample_id}",
            config["name_pattern"] + "_fastqc.zip"), sample_id=SAMPLES)
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
