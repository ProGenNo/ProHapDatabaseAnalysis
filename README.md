# ProHapDatabaseAnalysis
Code related to the publication "ProHap enables proteomic database generation accounting for population diversity". Pipeline to analyze protein databases created by ProHap.

## Requirements and Usage
Required software is Snakemake and Conda, remaining libraries are included in the provided Conda environment, created automatically by Snakemake.

Steps for reproducing results:
- Clone this repository
```
git clone https://github.com/ProGenNo/ProHapDatabaseAnalysis.git ;
cd ProHapDatabaseAnalysis ;
```
- Download and unpack supplementary data from [https://zenodo.org/records/10149278](https://zenodo.org/records/10149278)
```
cd data ;
wget https://zenodo.org/records/10149278/files/231220_ProHap_ALL_proteindb.tar.gz ; gunzip 231220_ProHap_ALL_proteindb.tar.gz ; tar xf 231220_ProHap_ALL_proteindb.tar ;
mkdir AFR ; cd AFR ;
wget https://zenodo.org/records/10149278/files/231221_ProHap_AFR_proteindb.tar.gz ; gunzip 231221_ProHap_AFR_proteindb.tar.gz ; tar xf 231221_ProHap_AFR_proteindb.tar ;
cd .. ; mkdir AMR ; cd AMR ;
wget https://zenodo.org/records/10149278/files/231221_ProHap_AMR_proteindb.tar.gz ; gunzip 231221_ProHap_AMR_proteindb.tar.gz ; tar xf 231221_ProHap_AMR_proteindb.tar ;
cd .. ; mkdir EUR ; cd EUR ;
wget https://zenodo.org/records/10149278/files/231221_ProHap_EUR_proteindb.tar.gz ; gunzip 231221_ProHap_EUR_proteindb.tar.gz ; tar xf 231221_ProHap_EUR_proteindb.tar ;
cd .. ; mkdir EAS ; cd EAS ;
wget https://zenodo.org/records/10149278/files/231221_ProHap_EAS_proteindb.tar.gz ; gunzip 231221_ProHap_EAS_proteindb.tar.gz ; tar xf 231221_ProHap_EAS_proteindb.tar ;
cd .. ; mkdir SAS ; cd SAS ;
wget https://zenodo.org/records/10149278/files/231221_ProHap_SAS_proteindb.tar.gz ; gunzip 231221_ProHap_SAS_proteindb.tar.gz ; tar xf 231221_ProHap_SAS_proteindb.tar ;
cd ../.. ;
```
- Run the pipeline using `snakemake -c<# cores> -p --use-conda`
