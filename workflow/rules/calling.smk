rule Module_2_Calling_Step_2:
    """
    In Module 2, we begin conducting some analyses in parallel. In the example,
    we process the data and perform variant detection and allele frequency 
    estimation in 5 million basepair non-overlapping windows.
    """
    input:
        os.path.join(OUTPUT_DIR, "all.bam.list")
    output:
        vcf=os.path.join(OUTPUT_DIR, "calls", "{chr_id}", "{chr_id}_{start}_{end}", "basevar.{chr_id}_{start}_{end}.vcf.gz"),
        tbi=os.path.join(OUTPUT_DIR, "calls", "{chr_id}", "{chr_id}_{start}_{end}", "basevar.{chr_id}_{start}_{end}.vcf.gz.tbi"),
        cvg=os.path.join(OUTPUT_DIR, "calls", "{chr_id}", "{chr_id}_{start}_{end}", "basevar.{chr_id}_{start}_{end}.cvg.gz")
    params:
        region="{chr_id}:{start}-{end}",
        outprefix=os.path.join(OUTPUT_DIR, "calls", "{chr_id}", "{chr_id}_{start}_{end}", "basevar.{chr_id}_{start}_{end}")
    threads: 1 # Multiple threads will cause std::bad_alloc error.
    resources:
        # Each thread (-t/--thread) typically necessitates
        # only 3GB to 4GB of memory if the -B (--batch-count)
        # parameter is set to 200.
        # See https://github.com/ShujiaHuang/basevar/issues/14
        mem_mb=4*1024
    log:
        get_log_path("{chr_id}_{start}_{end}")
    benchmark:
        BENCH_DIR + "/BaseVarC.basetype/{chr_id}_{start}_{end}.benchmark.txt"
    shell:
        """
        ./bin/basevar basetype \
            -t {threads} \
            -L {input} \
            --filename-has-samplename \
            -R {config[ref]} \
            -r {params.region:q} \
            --min-af=0.001 \
            --output-vcf {output.vcf} \
            --output-cvg {output.cvg} > {log} 2> {log}
        """

CHROM_SEGMENTS_BASEVAR=list(generate_chromosome_segments(exclude=False))

rule Module_2_Calling_Step_3:
    input:
        vcf=[os.path.join(OUTPUT_DIR, "calls",
                f"{chr_id}", f"{chr_id}_{start}_{end}", f"basevar.{chr_id}_{start}_{end}.vcf.gz")
                    for chr_id, start, end in CHROM_SEGMENTS_BASEVAR],
        tbi=[os.path.join(OUTPUT_DIR, "calls",
                f"{chr_id}", f"{chr_id}_{start}_{end}", f"basevar.{chr_id}_{start}_{end}.vcf.gz.tbi")
                    for chr_id, start, end in CHROM_SEGMENTS_BASEVAR]
    output:
        vcf=os.path.join(OUTPUT_DIR, "calls", "BaseVar.vcf.gz"),
        tbi=os.path.join(OUTPUT_DIR, "calls", "BaseVar.vcf.gz.tbi")
    threads: 12
    log: get_log_path("bcftools_concat_basevar")
    shell:
        """
        bcftools concat \
            --threads {threads} \
            -O z \
            -o {output.vcf} \
            {input.vcf} 2> {log}

        tabix -p vcf {output.vcf}
        """

rule Module_2_Generate_BaseVar_Posfile:
    """
    These position files were used in the two-stage approach of STITCH imputation.
    We don't do any filtration, because the first round of STITCH imputation will
    filter SNPs directly.
    """
    input:
        os.path.join(OUTPUT_DIR, "calls", "BaseVar.vcf.gz")
    output:
        os.path.join(OUTPUT_DIR, "calls", "{chr_id}", "BaseVar.{chr_id}.pos.txt")
    params:
        scr=r'{if(pre!=$2 && length($4) ==1 && length($5) ==1 && $4~/^[ACGT]$/ && $5~/^[ACGT]$/){print $1"\t"$2"\t"$4"\t"$5} pre=$2}'
    log:
        get_log_path("Module_2_Concat_BaseVar_Posfiles_{chr_id}")
    shell:
        """
        bcftools norm -r {wildcards.chr_id} --force -m- {input} 2> {log} | grep -v "#" | awk '{params.scr}' > {output}
        """
