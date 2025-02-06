## Prepare Data

### Human reference genome hg38 (GRCh38)

```shell
gsutil -m cp \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
  data/gatk
```

### GATK resource bundle

https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle

```shell
gsutil -m cp \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
  data/gatk
```

```shell
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > Homo_sapiens_assembly38.dbsnp138.vcf.gz
tabix -p vcf Homo_sapiens_assembly38.dbsnp138.vcf.gz
```

### 1000G Reference panel

Links:

- Used as reference panel: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
- Used to generate posfile: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/

Download all files under a remote folder:

```shell
wget -r -np -nH --cut-dirs=5 -R "index.html*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/

wget -r -np -nH --cut-dirs=5 -R "index.html*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
```

Rename files:

```shell
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
```

## Install Dependencies

### Download or Compile Tools

Put these tools in the `./bin` folder:

- https://github.com/rwdavies/STITCH/blob/1.6.8/STITCH.R
- https://github.com/Zilong-Li/BaseVarC

### Create Conda Environment

```shell
conda env create -f environment.yml
```

