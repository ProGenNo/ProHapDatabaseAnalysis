# ProHapDatabaseAnalysis
Code related to the publication "ProHap enables proteomic database generation accounting for population diversity". Pipeline to analyze protein databases created by ProHap.

## Requirements and Usage
Required software is Snakemake and Conda, remaining libraries are included in the provided Conda environment, created automatically by Snakemake.

Steps for reproducing results:
- Download and unpack supplementary data from [https://zenodo.org/records/10149278](https://zenodo.org/records/10149278) \([231220_ProHap_ALL_proteindb.tar.gz](https://zenodo.org/records/10149278/files/231220_ProHap_ALL_proteindb.tar.gz?download=1)\)
- Clone this repository
- Create a file called `config.yaml` in the root of this repository following instructions in `config_example.yaml`
- Run the pipeline using `snakemake -c<# cores> -p --use-conda`
