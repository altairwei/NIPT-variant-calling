rule Module_1_Alignment:
    input:
        fq=lambda wildcards: [s[0] for s in samples if s[3] == wildcards.samid][0],
        cid=lambda wildcards: [s[1] for s in samples if s[3] == wildcards.samid][0],
        lid=lambda wildcards: [s[2] for s in samples if s[3] == wildcards.samid][0]
    output:
        os.path.join(output_dir, "alignment/finaloutdir", "{samid}.sai"),
        os.path.join(output_dir, "alignment/finaloutdir", "{sample_id}.bam"),
        os.path.join(output_dir, "alignment/finaloutdir", "{sample_id}.sorted.bam"),
        os.path.join(output_dir, "alignment/finaloutdir", "{sample_id}.sorted.rmdup.bam"),
        os.path.join(output_dir, "alignment/finaloutdir", "{samid}.sorted.rmdup.realign.BQSR.bam")
    params:
        ref=ref,
        ref_index_prefix=ref_index_prefix,
        gatk_bundle_dir=gatk_bundle_dir,
        tmpoutdir=os.path.join(output_dir, "alignment/tmpoutdir"),
        finaloutdir=os.path.join(output_dir, "alignment/finaloutdir")
    shell:
        """
        ../bwa_alignment/alignment_workflow.sh \
            -f {input.fq} -c {input.cid} -l {input.lid} -s {wildcards.samid} \
            -g $CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar -j java \
            -a bwa -t samtools -e bedtools -z bgzip -x tabix \
            -r {params.ref} -i {params.ref_index_prefix} \
            -b {params.gatk_bundle_dir} \
            -o {params.tmpoutdir} -u {params.finaloutdir}
        """"