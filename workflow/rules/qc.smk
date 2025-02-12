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

rule Module_1_Statistics_Step_1:
    """
    Use Samtools to calculate alignment statistics for the alignment files.
    """
    input:
        bam=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam"),
        bai=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam.bai")
    output:
        temp(os.path.join(OUTPUT_DIR, "report", "alignment", "{sample_id}.sorted.rmdup.BQSR.bamstats"))
    log:
        get_log_path("{sample_id}")
    shell:
        "samtools stats {input.bam} > {output} 2>> {log}"

rule Module_1_Statistics_Step_2:
    """
    Use Bedtools to calculate alignment statistics for the alignment files.
    """
    input:
        bam=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam"),
        bai=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam.bai")
    output:
        bgzip=os.path.join(OUTPUT_DIR, "report", "alignment", "{sample_id}.sorted.rmdup.BQSR.cvg.bed.gz"),
        tabix=os.path.join(OUTPUT_DIR, "report", "alignment", "{sample_id}.sorted.rmdup.BQSR.cvg.bed.gz.tbi")
    log:
        get_log_path("{sample_id}")
    shell:
        "bedtools genomecov -ibam {input.bam} -bga -split 2>> {log} | bgzip > {output.bgzip} && tabix -p bed {output.bgzip}"

rule Module_1_Alignment_MultiQC:
    input:
        expand(os.path.join(OUTPUT_DIR, "report", "alignment", "{sample_id}.sorted.rmdup.BQSR.bamstats"), sample_id=SAMPLES)
    output:
        protected(os.path.join(OUTPUT_DIR, "report", "alignment.multiqc.html"))
    log:
        get_log_path("multiqc_alignment")
    params:
        searchdirs=os.path.join(OUTPUT_DIR, "report", "alignment"),
        outdir=os.path.join(OUTPUT_DIR, "report"),
        filename="alignment.multiqc.html"
    shell:
        """
        multiqc \
            --force \
            -o {params.outdir} \
            -n {params.filename} \
            --export \
            {params.searchdirs} > {log} 2> {log}
        """
