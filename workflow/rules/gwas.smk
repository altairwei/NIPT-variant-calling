rule Module_4_RunPCA:
    input:
        vcf=os.path.join(OUTPUT_DIR, "calls", "BaseVar.vcf.gz"),
        tbi=os.path.join(OUTPUT_DIR, "calls", "BaseVar.vcf.gz.tbi"),
        sexlist=os.path.join(OUTPUT_DIR, "all.samplesex.list")
    output:
        os.path.join(OUTPUT_DIR, "GWAS", "gwas.pca")
    threads: 16
    log: get_log_path("plink2_run_pca")
    shell:
        """
        plink2 \
            --pca 5 \
            --vcf {input.vcf} \
            --vcf-half-call haploid \
            --split-par hg38 \
            --update-sex {input.sexlist} \
            --out {output} \
            --threads {threads} 2> {log}
        """
