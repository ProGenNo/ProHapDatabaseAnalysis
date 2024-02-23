# ProHapDatabaseAnalysis
Code related to the publication "ProHap enables proteomic database generation accounting for population diversity": [doi.org/10.1101/2023.12.24.572591](https://doi.org/10.1101/2023.12.24.572591). Pipeline to analyze protein databases created by ProHap.

## Requirements and Usage
Required software is Snakemake and Conda \(install Conda / Mamba and Snakemake using [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba)\), remaining libraries are included in the provided Conda environment, created automatically by Snakemake.

### Steps for reproducing results:
- Clone this repository
```
git clone https://github.com/ProGenNo/ProHapDatabaseAnalysis.git ;
cd ProHapDatabaseAnalysis ;
```
- Download and unpack supplementary data from [https://zenodo.org/records/10688618](https://zenodo.org/records/10688618)
```
cd data ;
wget https://zenodo.org/records/10688618/files/240216_ProHap_ALL_proteindb.tar.gz ; gunzip 240216_ProHap_ALL_proteindb.tar.gz ; tar xf 240216_ProHap_ALL_proteindb.tar ;
mkdir AFR ; cd AFR ;
wget https://zenodo.org/records/10688618/files/240216_ProHap_AFR_proteindb.tar.gz ; gunzip 240216_ProHap_AFR_proteindb.tar.gz ; tar xf 240216_ProHap_AFR_proteindb.tar ;
cd .. ; mkdir AMR ; cd AMR ;
wget https://zenodo.org/records/10688618/files/240216_ProHap_AMR_proteindb.tar.gz ; gunzip 240216_ProHap_AMR_proteindb.tar.gz ; tar xf 240216_ProHap_AMR_proteindb.tar ;
cd .. ; mkdir EUR ; cd EUR ;
wget https://zenodo.org/records/10688618/files/240216_ProHap_EUR_proteindb.tar.gz ; gunzip 240216_ProHap_EUR_proteindb.tar.gz ; tar xf 240216_ProHap_EUR_proteindb.tar ;
cd .. ; mkdir EAS ; cd EAS ;
wget https://zenodo.org/records/10688618/files/240216_ProHap_EAS_proteindb.tar.gz ; gunzip 240216_ProHap_EAS_proteindb.tar.gz ; tar xf 240216_ProHap_EAS_proteindb.tar ;
cd .. ; mkdir SAS ; cd SAS ;
wget https://zenodo.org/records/10688618/files/240216_ProHap_SAS_proteindb.tar.gz ; gunzip 240216_ProHap_SAS_proteindb.tar.gz ; tar xf 240216_ProHap_SAS_proteindb.tar ;
cd ../.. ;
```
- By default, the pipeline uses at most 20 CPU cores per job. To change this, edit the `config.yaml` file.
- Run the pipeline using `snakemake --cores 20 -p --use-conda`. Change the `--cores` parameter to set the maximum amount of CPU cores.

### Result files
- `results/peptide_list_ALL.tsv`: list of all peptides included in the database (created using all samples in the 1000 Genomes panel)
- `results/all_discoverable_variants.tsv`: list of all variants discoverable in peptides
- `results/all_included_variants.csv`: list of all genetic variants included in the database
- `results/variant_stats.txt`: discoverability of variants in peptides
- `results/peptide_coverage_stats.tsv`: coverage of the proteome by canonical / variant peptides (used for Figure 2A)
- `results/peptide_coverage_freq_log.tsv`: coverage of the proteome by canonical / variant peptides based on haplotype frequency threshold (log. scale, used for Figure 2B)
- `results/transcript_freqs_by_superpop.tsv`: list of haplotype frequencies per superpopulation for each transcript (used for Figure 2C-G)

### Reproducing plots in Figure 2
Due to the computational demands, we expect the Snakemake pipeline to be run on a high-performance server. To reproduce the plots after running the entire pipeline, clone this repository again to a computer allowing work with Jupyter notebooks.

After cloning the repository, create a directory called `results`. Copy the following files to this directory: `results/peptide_coverage_stats.tsv`, `results/peptide_coverage_freq_log.tsv`, `results/transcript_freqs_by_superpop.tsv`.

To reproduce the plots in Figure 2, run the cells of the notebook in `src/Populations_comparison_plots.ipynb`.
