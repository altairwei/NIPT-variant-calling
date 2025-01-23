## Prepare Data

- https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
- https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
- https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/

## Install Dependencies

### Download or Compile Tools

Put these tools in the `./bin` folder:

- https://github.com/rwdavies/STITCH/blob/1.6.8/STITCH.R
- https://github.com/Zilong-Li/BaseVarC

### Create Conda Environment

```shell
conda env create -f environment.yml
```

