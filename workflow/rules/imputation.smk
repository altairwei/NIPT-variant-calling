## Step 1: Preparing the reference panel and imputation chunks

rule Module_3_Imputation_Step_1_1:
    """Download 1KGP reference panel"""
    output:
        vcf=protected(os.path.join(OUTPUT_DIR, "imputation", "refpanel", "CCDG_14151_B01_GRM_WGS_2020-08-05_{chr_id}.filtered.shapeit2-duohmm-phased.vcf.gz")),
        tbi=protected(os.path.join(OUTPUT_DIR, "imputation", "refpanel", "CCDG_14151_B01_GRM_WGS_2020-08-05_{chr_id}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"))
    params:
        url_prefix="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/"
        target_path=os.path.join(OUTPUT_DIR, "imputation", "refpanel")
    shell:
        """
        wget -c {params.url_prefix}{output.vcf} {params.url_prefix}{output.tbi} -P {params.target_path}
        """

rule Module_3_Imputation_Step_1_2:
    """Conduct normalization and filtration of the reference panel"""
    