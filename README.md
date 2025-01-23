## Prepare Data

- GATK resource bundle: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
- 1000G Reference panel: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
- 1000G Reference panel: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/

```shell
wget -r -np -nH --cut-dirs=5 -R "index.html*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/

wget -r -np -nH --cut-dirs=5 -R "index.html*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
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

