rule Module_4_RunPCA:
    input:
        vcf=os.path.join(OUTPUT_DIR, "calls", "BaseVarC.vcf.gz"),
        tbi=os.path.join(OUTPUT_DIR, "calls", "BaseVarC.vcf.gz.tbi")
    output:
        os.path.join(OUTPUT_DIR, "GWAS", "gwas.pca")
    threads: 16
    shell:
        """
        plink2 \
            --pca 5 \
            --vcf {input.vcf} \
            --out {output} \
            --threads {threads}
        """
