## Prepare Data

### Human reference genome hg38 (GRCh38)

Download the analysis set of hg38:

```shell
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.{gz,fai,bwa_index.tar.gz}
```

Create GATK dictionary:

```shell
gatk CreateSequenceDictionary -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

### GATK resource bundle

Download the latest dbSNP:

```shell
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz{,.md5}
md5sum -c GCF_000001405.40.gz.md5
```

Convert contig IDs to UCSC using [chromToUcsc](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc):

```shell
chromToUcsc --get hg38
gunzip -c GCF_000001405.40.gz | chromToUcsc -a hg38.chromAlias.tsv | bgzip -c > dbsnp_156.hg38.vcf.gz
tabix -p vcf dbsnp_156.hg38.vcf.gz
```

https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle

```shell
gsutil -m cp \
  gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
  gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
  gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
  data/gatk
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

Put or compile these tools in the `./bin` folder:

- https://github.com/rwdavies/STITCH/blob/1.7.3/STITCH.R
- https://github.com/altairwei/basevar
- http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc

### Create Conda Environment

```shell
conda env create -f environment.yml
```

