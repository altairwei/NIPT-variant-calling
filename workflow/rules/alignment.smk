import datetime

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

configfile: "config/config.yml"

output_dir = config["output_dir"]
gatk_bundle_dir = config["gatk_bundle_dir"]
ref = config["ref"]
ref_index_prefix = config["ref_index_prefix"]

common_params = {
    "ref": ref,
    "gatk_bundle_dir": gatk_bundle_dir
}

def get_log_path(sample_id):
    return os.path.join(config["log_dir"], timestamp, f"{sample_id}.log")

### Step 1: Read alignment using bwa
### - Adjust the number of threads (-t) for bwa according to your cluster.
### - The -e 10 parameter makes indel detection more sensitive, though its impact is minimal.

rule Module_1_Alignment_Step_1_1:
    """
    Used the BWA single-end alignment model to map the single-end reads 
    (typically 35 bp) to the latest human genome reference.
    """
    input:
        "data/{sample_id}.fq.gz"
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sai")
    log:
        get_log_path("{sample_id}")
    params:
        ref_index_prefix=ref_index_prefix
    threads: config.get("threads", 8)
    shell:
        "bwa aln -e 10 -t {threads} -i 5 -q 0 {params.ref_index_prefix} {input} > {output} 2>> {log}"

rule Module_1_Alignment_Step_1_2:
    """
    Use the BWA single-end alignment model to map the single-end reads 
    (typically 35 bp) to the latest human genome reference.
    """
    input:
        fq="data/{cid}_{lid}_{snn}.fq.gz",
        sai=os.path.join(output_dir, "alignment", "{cid}_{lid}_{snn}.sai")
    output:
        os.path.join(output_dir, "alignment", "{cid}_{lid}_{snn}.bam")
    log:
        get_log_path("{cid}_{lid}_{snn}")
    params:
        ref_index_prefix=ref_index_prefix,
        read_group="@RG\\tID:{cid}_{lid}\\tPL:COMPLETE\\tSM:{cid}_{lid}_{snn}"
    shell:
        """
        bwa samse -r "{params.read_group}" {params.ref_index_prefix} {input.sai} {input.fq} 2>> {log} \
            | samtools view -h -Sb - > {output}
        """

rule Module_1_Alignment_Step_1_3:
    """
    The alignment reads were then sorted using Samtools.
    """
    input:
        os.path.join(output_dir, "alignment", "{sample_id}.bam")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.bam")
    threads: config.get("threads", 8)
    log:
        get_log_path("{sample_id}")
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input} 2>> {log}"

rule Module_1_Alignment_Step_1_4:
    """
    Removal of potential PCR duplicates.
    """
    input:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.bam")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.bam")
    log:
        get_log_path("{sample_id}")
    shell:
        "samtools rmdup {input} {output} 2>> {log}"

rule Module_1_Alignment_Step_1_5:
    input:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.bam")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.bam.bai")
    shell:
        "samtools index {input}"

### Step 2: Re-alignment with GATK
### realignment, adjust the memory usage according to data amount

rule Module_1_Alignment_Step_2_1:
    """
    Use GATK to realign indels in the NIPT reads based on known indel information from prior studies.
    """
    input:
        bam=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.bam"),
        bai=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.bam.bai")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.indel_target_intervals.list")
    params:
        **common_params
    log:
        get_log_path("{sample_id}")
    shell:
        """
        java -Xmx15g -jar $CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar \
            -T RealignerTargetCreator \
            -R {params.ref} \
            -I {input.bam} \
            -known {params.gatk_bundle_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            -known {params.gatk_bundle_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
            -o {output} 2>> {log}
        """

rule Module_1_Alignment_Step_2_2:
    input:
        bam=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.bam"),
        bai=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.bam.bai"),
        indel=os.path.join(output_dir, "alignment", "{sample_id}.indel_target_intervals.list")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.bam")
    params:
        **common_params
    log:
        get_log_path("{sample_id}")
    shell:
        """
        java -Xmx15g -jar $CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar \
            -T IndelRealigner \
            -R {params.ref} \
            -I {input.bam} \
            -known {params.gatk_bundle_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            -known {params.gatk_bundle_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
            --targetIntervals {input.indel} \
            --disable_bam_indexing \
            -o {output} 2>> {log}
        """

rule Module_1_Alignment_Step_2_3:
    input:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.bam")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.bam.bai")
    shell:
        "samtools index {input}"

### Step 3: BQSR base quality score recalibration

rule Module_1_Alignment_Step_3_1:
    """
    Recalibrate base quality in the NIPT reads, using known SNPs and indels as references.
    dbsnp_146.hg38.vcf.gz is available in http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/
    """
    input:
        bam=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.bam"),
        bai=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.bam.bai")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.recal_data.table")
    params:
        **common_params
    log:
        get_log_path("{sample_id}")
    shell:
        """
        java -jar $CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar \
            -T BaseRecalibrator \
            -nct 8 \
            -R {params.ref} \
            -I {input.bam} \
            --knownSites {params.gatk_bundle_dir}/dbsnp_146.hg38.vcf.gz \
            --knownSites {params.gatk_bundle_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --knownSites {params.gatk_bundle_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
            -o {output} 2>> {log}
        """

rule Module_1_Alignment_Step_3_2:
    """
    Recalibrate base quality in the NIPT reads, using known SNPs and indels as references.
    """
    input:
        tbl=os.path.join(output_dir, "alignment", "{sample_id}.recal_data.table"),
        bam=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.bam"),
        bai=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.bam.bai")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bam")
    params:
        **common_params
    log:
        get_log_path("{sample_id}")
    shell:
        """
        java -jar $CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar \
            -T PrintReads \
            -nct 8 \
            -R {params.ref} \
            --BQSR {input.tbl} \
            -I {input.bam} \
            -o {output} 2>> {log}
        """

rule Module_1_Alignment_Step_3_3:
    input:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bam")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bam.bai")
    shell:
        "samtools index {input}"

### Step 4. Bam statistics with samtools and bedtools

rule Module_1_Statistics_Step_4_1:
    """
    Use Samtools to calculate alignment statistics for the alignment files.
    """
    input:
        bam=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bam"),
        bai=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bam.bai")
    output:
        os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bamstats")
    log:
        get_log_path("{sample_id}")
    shell:
        "samtools stats {input.bam} > {output} 2>> {log}"

rule Module_1_Statistics_Step_4_2:
    """
    Use Bedtools to calculate alignment statistics for the alignment files.
    """
    input:
        bam=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bam"),
        bai=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.bam.bai")
    output:
        bgzip=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz"),
        tabix=os.path.join(output_dir, "alignment", "{sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz.tbi")
    log:
        get_log_path("{sample_id}")
    shell:
        "bedtools genomecov -ibam {input.bam} -bga -split 2>> {log} | bgzip > {output.bgzip} && tabix -p bed {output.bgzip}"

