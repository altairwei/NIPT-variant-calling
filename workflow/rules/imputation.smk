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

rule Module_3_ConcatPosfiles:
    input:
        expand(os.path.join(OUTPUT_DIR, "imputation", "{chrom}", "Panel.{chrom}.pos.txt"),
            chrom=[f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        os.path.join(OUTPUT_DIR, "imputation", "Panel.allchrom.pos.txt")
    shell:
        "cat {input} > {output}"

rule Module_3_STITCH_PreImputation_Step_1:
    """
    STITCH itself can be used to filter SNPs.

    In the first step, all discovered sites in the sample are imputed without
    filtration. In the second step, only variants that pass the quality control 
    filters from the first step are reimputed.
    """
    input:
        bamlist=os.path.join(OUTPUT_DIR, "all.bam.list"),
        snlist=os.path.join(OUTPUT_DIR, "all.samplename.list"),
        posfile=os.path.join(OUTPUT_DIR, "calls", "{chr}", "BaseVar.{chr}.pos.txt")
    output:
        vcf=os.path.join(OUTPUT_DIR, "preimputation", "{chr}", "{chr}_{start}_{end}", "stitch.{chr}.{start}.{end}.vcf.gz"),
        idx=os.path.join(OUTPUT_DIR, "preimputation", "{chr}", "{chr}_{start}_{end}", "stitch.{chr}.{start}.{end}.vcf.gz.tbi")
    params:
        region="--regionStart {start} --regionEnd {end} --chr {chr}",
        outdir=os.path.join(OUTPUT_DIR, "preimputation", "{chr}", "{chr}_{start}_{end}")
    log:
        get_log_path("{chr}_{start}_{end}")
    benchmark:
        BENCH_DIR + "/STITCH.PreImputation/{chr}_{start}_{end}.benchmark.txt"
    shell:
        """
        ./bin/STITCH.R \
            --outputdir {params.outdir} \
            --bamlist {input.bamlist} \
            --sampleNames_file {input.snlist} \
            --posfile {input.posfile} \
            --reference {config[ref]} \
            --K 10 --nGen 16000 --nCores 1 \
            {params.region} \
            --buffer 250000 2> {log} >> {log}

        bcftools index -t {output.vcf}
        """

rule Module_3_STITCH_PreImputation_Step_2:
    """Generate chromosome-level position files."""
    input:
        vcf=lambda wildcards: [os.path.join(OUTPUT_DIR, "preimputation",
            f"{chrom}", f"{chrom}_{start}_{end}", f"stitch.{chrom}.{start}.{end}.vcf.gz")
                for chrom, start, end in list(generate_chromosome_segments(chroms=[wildcards.chr]))]
    output:
        os.path.join(OUTPUT_DIR, "preimputation", "{chr}", "stitch.preimputation.{chr}.pos.txt")
    log:
        get_log_path("Module_3_STITCH_PreImputation_Step_2_{chr}")
    params:
        # These cutoffs come from the STITCH article
        filter='INFO/INFO_SCORE>0.4 && INFO/HWE>0.000001'
    shell:
        """
        bcftools concat -Ov {input.vcf} 2> {log} \
            | bcftools view -i '{params.filter}' -H -Ov - 2>> {log} \
            | cut -f '1,2,4,5' > {output} 2>> {log}
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
        posfile=os.path.join(OUTPUT_DIR, "preimputation", "{chr}", "stitch.preimputation.{chr}.pos.txt")
    output:
        vcf=os.path.join(OUTPUT_DIR, "imputation", "{chr}", "{chr}_{start}_{end}", "stitch.{chr}.{start}.{end}.vcf.gz"),
        idx=os.path.join(OUTPUT_DIR, "imputation", "{chr}", "{chr}_{start}_{end}", "stitch.{chr}.{start}.{end}.vcf.gz.tbi")
    params:
        region="--regionStart {start} --regionEnd {end} --chr {chr}",
        outdir=os.path.join(OUTPUT_DIR, "imputation", "{chr}", "{chr}_{start}_{end}")
    log:
        get_log_path("{chr}_{start}_{end}")
    benchmark:
        BENCH_DIR + "/STITCH.R/{chr}_{start}_{end}.benchmark.txt"
    shell:
        """
        ./bin/STITCH.R \
            --outputdir {params.outdir} \
            --bamlist {input.bamlist} \
            --sampleNames_file {input.snlist} \
            --posfile {input.posfile} \
            --reference {config[ref]} \
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
        vcf=os.path.join(OUTPUT_DIR, "imputation", "STITCH.vcf.gz"),
        tbi=os.path.join(OUTPUT_DIR, "imputation", "STITCH.vcf.gz.tbi")
    threads: 16
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

rule Module_3_Generate_STITCH_Posfile:
    """
    These position files were used in the two-stage approach of STITCH imputation.
    We don't do any filtration, because the first round of STITCH imputation will
    filter SNPs directly.
    """
    input:
        os.path.join(OUTPUT_DIR, "imputation", "STITCH.vcf.gz")
    output:
        os.path.join(OUTPUT_DIR, "imputation", "{chr_id}", "STITCH.imputation.{chr_id}.pos.txt")
    log:
        get_log_path("Module_3_Generate_STITCH_Posfile_{chr_id}")
    params:
        # These cutoffs come from the STITCH article
        filter='INFO/INFO_SCORE>0.4 && INFO/HWE>0.000001'
    shell:
        """
        bcftools view -r {wildcards.chr_id} -i '{params.filter}' -H -Ov {input} 2>> {log} \
            | cut -f '1,2,4,5' > {output} 2>> {log}
        """
