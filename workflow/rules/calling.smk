rule Module_2_Calling_Step_1:
    input:
        bam=expand(os.path.join(OUTPUT_DIR, "alignments", 
            "{sample_id}.sorted.rmdup.realign.BQSR.bam"), sample_id=[s[3] for s in SAMPLES]),
        bai=expand(os.path.join(OUTPUT_DIR, "alignments",
            "{sample_id}.sorted.rmdup.realign.BQSR.bam.bai"), sample_id=[s[3] for s in SAMPLES])
    output:
        temp(os.path.join(OUTPUT_DIR, "calls", "all.bam.list"))
    run:
        with open(output[0], "w") as f:
            for bam_file in input.bam:
                f.write(bam_file + "\n")

rule Module_2_Calling_Step_2:
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

def generate_regional_vcf_files(wildcards, suffix=".vcf.gz"):
    delta = 5000000
    in_fai = config["ref_fai"]
    chroms = wildcards.chr_id

    with open(in_fai) as fh:
        for r in fh:
            col = r.strip().split()
            chr_id = col[0]
            chr_length = int(col[1])

            if chroms is not None and len(chroms):
                if chr_id not in chroms:
                    continue

            for i in range(0, chr_length, delta):
                start = i + 1
                end = i + delta if i + delta <= chr_length else chr_length
                yield os.path.join(OUTPUT_DIR, "calls", f"{chr_id}_{start}_{end}{suffix}")


rule Module_2_Calling_Step_3:
    input:
        vcf=lambda wildcards: generate_regional_vcf_files(wildcards, ".vcf.gz"),
        tbi=lambda wildcards: generate_regional_vcf_files(wildcards, ".vcf.gz.tbi")
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

rule Module_2_Calling_Step_4:
    input:
        cvg=lambda wildcards: generate_regional_vcf_files(wildcards, ".cvg.tsv.gz"),
        tbi=lambda wildcards: generate_regional_vcf_files(wildcards, ".cvg.tsv.gz.tbi")
    output:
        cvg=protected(os.path.join(OUTPUT_DIR, "calls", "{chr_id}.cvg.tsv.gz")),
        tbi=protected(os.path.join(OUTPUT_DIR, "calls", "{chr_id}.cvg.tsv.gz.tbi"))
    shell:
        "bedtools unionbedg -i {input} > {output.cvg} && tabix -p bed {output.cvg}"
