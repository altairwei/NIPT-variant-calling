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

rule Module_1_Alignment_Step_1:
    """
    Used the BWA single-end alignment model to map the single-end reads 
    (typically 35 bp) to the latest human genome reference. The -e 10 parameter
    makes indel detection more sensitive, though its impact is minimal.
    """
    input:
        os.path.join(OUTPUT_DIR, "clean", "{sample_id}.clean.fq.gz")
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sai"))
    log:
        get_log_path("{sample_id}")
    params:
        ref_index_prefix=REF_INDEX_PREFIX
    threads: 4
    shell:
        "bwa aln -e 10 -t {threads} -i 5 -q 0 {params.ref_index_prefix} {input} > {output} 2>> {log}"

rule Module_1_Alignment_Step_2:
    """
    Use the BWA single-end alignment model to map the single-end reads 
    (typically 35 bp) to the latest human genome reference.
    """
    input:
        fq=os.path.join(OUTPUT_DIR, "clean", "{sample_id}.clean.fq.gz"),
        sai=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sai")
    output:
        temp(os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.bam"))
    log:
        get_log_path("{sample_id}")
    params:
        ref_index_prefix=REF_INDEX_PREFIX,
        read_group=f"@RG\\tID:{{sample_id}}\\tPL:{PLATFORM}\\tSM:{{sample_id}}"
    shell:
        """
        bwa samse -r {params.read_group:q} {params.ref_index_prefix} {input.sai} {input.fq} 2>> {log} \
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
    params:
        **common_params
    resources:
        mem_mb=2*1024
    log:
        get_log_path("{sample_id}")
    shell:
        """
        gatk BaseRecalibrator \
            --java-options "-Xmx{resources.mem_mb}m" \
            -R {params.ref} \
            -I {input.bam} \
            --known-sites {params.gatk_bundle_dir}/dbsnp_146.hg38.vcf.gz \
            -O {output} 2>> {log}
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
    params:
        **common_params
    log:
        get_log_path("{sample_id}")
    shell:
        """
        gatk ApplyBQSR \
            -R {params.ref} \
            --bqsr-recal-file {input.tbl} \
            -I {input.bam} \
            -O {output} 2>> {log}
        """

rule Module_1_PostRecalibration:
    input:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam")
    output:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam.bai")
    shell:
        "samtools index {input}"

rule Module_1_Statistics_Step_1:
    """
    Use Samtools to calculate alignment statistics for the alignment files.
    """
    input:
        bam=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam"),
        bai=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bam.bai")
    output:
        os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.bamstats")
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
        bgzip=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.cvg.bed.gz"),
        tabix=os.path.join(OUTPUT_DIR, "alignments", "{sample_id}.sorted.rmdup.BQSR.cvg.bed.gz.tbi")
    log:
        get_log_path("{sample_id}")
    shell:
        "bedtools genomecov -ibam {input.bam} -bga -split 2>> {log} | bgzip > {output.bgzip} && tabix -p bed {output.bgzip}"

