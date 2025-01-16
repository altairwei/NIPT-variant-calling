rule Module_2_Calling_Step_2:
    """
    In Module 2, we begin conducting some analyses in parallel. In the example,
    we process the data and perform variant detection and allele frequency 
    estimation in 5 million basepair non-overlapping windows. To facilitate this
    parallelization, we use the pipeline generator create_pipeline.py, which
    distributes the computational tasks based on the --delta parameter across a 
    specific chromosome defined by the -c parameter.
    """
    input:
        os.path.join(OUTPUT_DIR, "calls", "all.bam.list")
    output:
        vcf=temp(os.path.join(OUTPUT_DIR, "calls", "{chr_id}_{start}_{end}.vcf.gz")),
        vcf_tbi=temp(os.path.join(OUTPUT_DIR, "calls", "{chr_id}_{start}_{end}.vcf.gz.tbi")),
        cvg=temp(os.path.join(OUTPUT_DIR, "calls", "{chr_id}_{start}_{end}.cvg.tsv.gz")),
        cvg_tbi=temp(os.path.join(OUTPUT_DIR, "calls", "{chr_id}_{start}_{end}.cvg.tsv.gz.tbi"))
    params:
        ref=REF,
        region="{chr_id}:{start}-{end}"
    threads: 1 # Multiple threads will cause std::bad_alloc error.
    resources:
        # Each thread (-t/--thread) typically necessitates
        # only 3GB to 4GB of memory if the -B (--batch-count)
        # parameter is set to 200.
        # See https://github.com/ShujiaHuang/basevar/issues/14
        mem_mb=4*1024
    log:
        get_log_path("{chr_id}_{start}_{end}")
    shell:
        """
        ./bin/basevar basetype -t {threads} \
            -L {input} \
            --filename-has-samplename \
            -R {params.ref} \
            -r {params.region:q} \
            --min-af=0.001 \
            --output-vcf {output.vcf} \
            --output-cvg {output.cvg} > {log}
        """

rule Module_2_Calling_Step_3:
    input:
        vcf=lambda wildcards: generate_regional_vcf_files(
            wildcards.chr_id, os.path.join(OUTPUT_DIR, "calls"), ".vcf.gz"),
        tbi=lambda wildcards: generate_regional_vcf_files(
            wildcards.chr_id, os.path.join(OUTPUT_DIR, "calls"), ".vcf.gz.tbi")
    output:
        vcf=protected(os.path.join(OUTPUT_DIR, "calls", "{chr_id}.vcf.gz")),
        tbi=protected(os.path.join(OUTPUT_DIR, "calls", "{chr_id}.vcf.gz.tbi"))
    threads: 12
    shell:
        """
        bcftools concat \
            --threads {threads} \
            -a --rm-dups all -O z \
            -o {output.vcf} \
            {input.vcf}

        tabix -p vcf {output.vcf}
        """
