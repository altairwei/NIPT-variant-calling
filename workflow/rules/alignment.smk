rule Module_1_Alignment_Step_1:
    """
    Used the BWA single-end alignment model to map the single-end reads 
    (typically 35 bp) to the latest human genome reference. The -e 10 parameter
    makes indel detection more sensitive, though its impact is minimal.
    """
    input:
        config["file_pattern"]
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sai"))
    log:
        get_log_path("{sample_id}")
    threads: 4
    shell:
        "bwa aln -e 10 -t {threads} -i 5 -q 0 {config[ref]} {input} > {output} 2>> {log}"

rule Module_1_Alignment_Step_2:
    """
    Use the BWA single-end alignment model to map the single-end reads 
    (typically 35 bp) to the latest human genome reference.
    """
    input:
        fq=config["file_pattern"],
        sai=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sai")
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.bam"))
    log:
        get_log_path("{sample_id}")
    params:
        read_group=f"@RG\\tID:{{sample_id}}\\tPL:{PLATFORM}\\tSM:{{sample_id}}"
    shell:
        """
        bwa samse -r {params.read_group:q} {config[ref]} {input.sai} {input.fq} 2>> {log} \
            | samtools view -h -Sb - > {output}
        """

rule Module_1_PostAlignment_Step_1:
    """
    The alignment reads were then sorted using Samtools.
    """
    input:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.bam")
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.bam"))
    threads: 8
    log:
        get_log_path("{sample_id}")
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input} 2>> {log}"

rule Module_1_PostAlignment_Step_2:
    """
    Removal of potential PCR duplicates.
    """
    input:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.bam")
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.bam"))
    log:
        get_log_path("{sample_id}")
    shell:
        "samtools rmdup {input} {output} 2>> {log}"

rule Module_1_PostAlignment_Step_3:
    input:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.bam")
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.bam.bai"))
    shell:
        "samtools index {input}"

rule Module_1_Recalibration_Step_1:
    """
    Recalibrate base quality in the NIPT reads, using known SNPs and indels as references.
    dbsnp_146.hg38.vcf.gz is available in http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/
    """
    input:
        bam=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.bam"),
        bai=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.bam.bai")
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.recal_data.table"))
    resources:
        mem_mb=2*1024
    log:
        get_log_path("{sample_id}")
    benchmark:
        BENCH_DIR + "/gatk.BaseRecalibrator/{sample_id}.benchmark.txt"
    shell:
        """
        gatk BaseRecalibrator \
            --java-options "-Xmx{resources.mem_mb}m -XX:ConcGCThreads=1" \
            -R {config[ref]} \
            -I {input.bam} \
            --known-sites {config[dbSNP]} \
            --known-sites {config[gatk_bundle_dir]}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --known-sites {config[gatk_bundle_dir]}/Homo_sapiens_assembly38.known_indels.vcf.gz \
            -O {output} > {log} 2>&1
        """

rule Module_1_Recalibration_Step_2:
    """
    Recalibrate base quality in the NIPT reads, using known SNPs and indels as references.
    """
    input:
        tbl=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.recal_data.table"),
        bam=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.bam"),
        bai=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.bam.bai")
    output:
        protected(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam"))
    log:
        get_log_path("{sample_id}")
    shell:
        """
        gatk ApplyBQSR \
            --java-options "-XX:ConcGCThreads=1" \
            -R {config[ref]} \
            --bqsr-recal-file {input.tbl} \
            -I {input.bam} \
            --create-output-bam-index false \
            -O {output} > {log} 2>&1
        """

rule Module_1_PostRecalibration:
    input:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam")
    output:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam.bai")
    shell:
        "samtools index {input}"

rule Module_1_ListSamples:
    """
    Create a list of all the bam files.
    """
    input:
        bam=expand(os.path.join(OUTPUT_DIR, "alignments", 
            "{sample_id}.sorted.rmdup.BQSR.bam"), sample_id=SAMPLES),
        bai=expand(os.path.join(OUTPUT_DIR, "alignments",
            "{sample_id}.sorted.rmdup.BQSR.bam.bai"), sample_id=SAMPLES)
    output:
        bamlist=os.path.join(OUTPUT_DIR, "all.bam.list"),
        snlist=os.path.join(OUTPUT_DIR, "all.samplename.list"),
        sexlist=os.path.join(OUTPUT_DIR, "all.samplesex.list")
    run:
        with open(output.bamlist, "w") as f_bam:
            for bam_file in input.bam:
                f_bam.write(bam_file + "\n")
        with open(output.snlist, "w") as f_sample, open(output.sexlist, "w") as f_sex:
            for sample_id in SAMPLES:
                f_sample.write(sample_id + "\n")
                f_sex.write(sample_id + "\t" + "Female\n")

