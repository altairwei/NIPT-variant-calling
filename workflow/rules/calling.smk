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
        ./bin/BaseVarC basetype \
            --input {input} \
            --reference {config[ref]} \
            --region {params.region:q} \
            --output {params.outprefix} \
            --batch 200 \
            --thread {threads} > {log} 2> {log}
        
        tabix -p vcf {output.vcf}
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
        vcf=protected(os.path.join(OUTPUT_DIR, "calls", "BaseVarC.vcf.gz")),
        tbi=protected(os.path.join(OUTPUT_DIR, "calls", "BaseVarC.vcf.gz.tbi"))
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
