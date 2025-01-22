rule Module_3_ListSamples:
    """
    Create a list of all the bam files.
    """
    input:
        bam=expand(os.path.join(OUTPUT_DIR, "alignments", 
            "{sample_id}.sorted.rmdup.BQSR.bam"), sample_id=[s[3] for s in SAMPLES]),
        bai=expand(os.path.join(OUTPUT_DIR, "alignments",
            "{sample_id}.sorted.rmdup.BQSR.bam.bai"), sample_id=[s[3] for s in SAMPLES])
    output:
        bamlist=temp(os.path.join(OUTPUT_DIR, "all.bam.list")),
        snlist=temp(os.path.join(OUTPUT_DIR, "all.samplename.list"))
    run:
        with open(output.bamlist, "w") as f_bam:
            for bam_file in input.bam:
                f_bam.write(bam_file + "\n")
        with open(output.snlist, "w") as f_sample:
            for sample in SAMPLES:
                sample_id = sample[3]
                f_sample.write(sample_id + "\n")

rule Module_3_GeneratePosfile:
    input:
        "data/refpanel/20220422_3202_phased_SNV_INDEL_SV/"
        "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    output:
        os.path.join(OUTPUT_DIR, "imputation", "{chr}", "Panel.{chr}.pos.txt")
    params:
        scr=r'{if(pre!=$2 && length($4) ==1 && length($5) ==1 ){print $1"\t"$2"\t"$4"\t"$5} pre=$2}'
    shell:
        """
        zcat {input} | grep -v "#" | awk '{params.scr}' > {output}
        """

rule Module_3_RunSTITCH_Step_1:
    """
    Given the large sample size used in STITCH, failing to segment chromosome 
    intervals in advance may lead to insufficient computing resources. Therefore,
    we need to divide chromosomes into smaller chunks.
    """
    input:
        bamlist=os.path.join(OUTPUT_DIR, "all.bam.list"),
        snlist=os.path.join(OUTPUT_DIR, "all.samplename.list"),
        posfile=os.path.join(OUTPUT_DIR, "imputation", "{chr}", "Panel.{chr}.pos.txt")
    output:
        vcf=os.path.join(OUTPUT_DIR, "imputation", "{chr}", "{chr}_{start}_{end}", "stitch.{chr}.{start}.{end}.vcf.gz"),
        idx=os.path.join(OUTPUT_DIR, "imputation", "{chr}", "{chr}_{start}_{end}", "stitch.{chr}.{start}.{end}.vcf.gz.tbi")
    params:
        ref=REF,
        region="--regionStart {start} --regionEnd {end} --chr {chr}",
        outdir=os.path.join(OUTPUT_DIR, "imputation", "{chr}", "{chr}_{start}_{end}")
    log:
        get_log_path("{chr}_{start}_{end}")
    shell:
        """
        ./bin/STITCH.R \
            --outputdir {params.outdir} \
            --bamlist {input.bamlist} \
            --sampleNames_file {input.snlist} \
            --posfile {input.posfile} \
            --reference {params.ref} \
            --K 10 --nGen 16000 --nCores 1 \
            {params.region} \
            --buffer 250000 2> {log} >> {log}
        
        bcftools index -t {output.vcf}
        """

CHROM_SEGMENTS=list(generate_chromosome_segments())

rule Module_3_RunSTITCH_Step_2:
    """
    Merge all genotype imputation results using BCFtools.

    Note: We ignored chromosome Y and Mitochondrial genome, since the imputation 
    on these genomes is complex.
    """
    input:
        vcf=[os.path.join(OUTPUT_DIR, "imputation",
                f"{chrom}", f"{chrom}_{start}_{end}", f"stitch.{chrom}.{start}.{end}.vcf.gz")
                    for chrom, start, end in CHROM_SEGMENTS],
        tbi=[os.path.join(OUTPUT_DIR, "imputation",
                f"{chrom}", f"{chrom}_{start}_{end}", f"stitch.{chrom}.{start}.{end}.vcf.gz.tbi")
                    for chrom, start, end in CHROM_SEGMENTS]
    output:
        vcf=protected(os.path.join(OUTPUT_DIR, "imputation", "STITCH.vcf.gz")),
        tbi=protected(os.path.join(OUTPUT_DIR, "imputation", "STITCH.vcf.gz.tbi"))
    threads: 4
    log: get_log_path("bcftools_concat_STITCH")
    shell:
        """
        bcftools concat \
            -Oz \
            --threads {threads} \
            -o {output.vcf} \
            {input.vcf} 2> {log}
        
        tabix -p vcf {output.vcf}
        """
