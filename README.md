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
- Download and unpack supplementary data from [https://doi.org/10.5281/zenodo.10149277](https://doi.org/10.5281/zenodo.10149277) \(1000 Genomes\), [https://doi.org/10.5281/zenodo.12671301](https://doi.org/10.5281/zenodo.12671301) \(HRC\), and [https://doi.org/10.5281/zenodo.12686818](https://doi.org/10.5281/zenodo.12686818) (Pangenome)
```
cd data ;
wget https://zenodo.org/records/12671237/files/240527_ProHap_ALL.tar.gz ; tar xf 240527_ProHap_ALL.tar.gz ;
mkdir AFR ; cd AFR ;
wget https://zenodo.org/records/12671237/files/240530_ProHap_AFR.tar.gz ; tar xf 240530_ProHap_AFR.tar.gz ;
cd .. ; mkdir AMR ; cd AMR ;
wget https://zenodo.org/records/12671237/files/240530_ProHap_AMR.tar.gz ; tar xf 240530_ProHap_AMR.tar.gz ;
cd .. ; mkdir EUR ; cd EUR ;
wget https://zenodo.org/records/12671237/files/240530_ProHap_EUR.tar.gz ; tar xf 240530_ProHap_EUR.tar.gz ;
cd .. ; mkdir EAS ; cd EAS ;
wget https://zenodo.org/records/12671237/files/240530_ProHap_EAS.tar.gz ; tar xf 240530_ProHap_EAS.tar.gz ;
cd .. ; mkdir SAS ; cd SAS ;
wget https://zenodo.org/records/12671237/files/240530_ProHap_SAS.tar.gz ; tar xf 240530_ProHap_SAS.tar.gz ;
cd .. ; mkdir HRC ; cd HRC ;
wget https://zenodo.org/records/12671302/files/240703_HRC1.1_GRCh38.tar.gz ; tar xf 240703_HRC1.1_GRCh38.tar.gz ;
cd .. ; mkdir pangenome ; cd pangenome ;
wget https://zenodo.org/records/12686819/files/240703_ProHap_pangenome_ALL.tar.gz ; tar xf 240703_ProHap_pangenome_ALL.tar.gz ;
cd ../.. ;
```
- By default, the pipeline uses at most 12 CPU cores per job. To change this, edit the `config.yaml` file.
- Run the pipeline using `snakemake --cores 12 -p --use-conda`. Change the `--cores` parameter to set the maximum amount of CPU cores.
- Hardware requirements: Up to 80 GB RAM, 50 GB disk space, configurable number of CPU cores. Note that the RAM requirements decrease the fewer CPU cores are used, as data are copied internally by Python when running functions in parallel.

### Result files
- `results/peptide_list_ALL.tsv`: list of all peptides included in the database (created using all samples in the 1000 Genomes panel)
- `results/all_discoverable_variants.tsv`: list of all variants discoverable in peptides
- `results/all_included_variants.csv`: list of all genetic variants included in the database
- `results/variant_stats.txt`: discoverability of variants in peptides
- `results/peptide_coverage_stats.tsv`: coverage of the proteome by canonical / variant peptides (used for Figure 2A)
- `results/peptide_coverage_freq_log.tsv`: coverage of the proteome by canonical / variant peptides based on haplotype frequency threshold (log. scale, used for Figure 2B)
- `results/transcript_freqs_by_superpop.tsv`: list of haplotype frequencies per superpopulation for each transcript (used for Figure 2C-G)
- `results/peptide_indiv_coverage.tsv`: coverage of the proteome by canonical / variant peptides per individual sample (1000 Genomes)
- `results/uniprot_comparison_stats.tsv`: size of overlap between tryptic peptides in the SwissProt, UniProt Isoforms and all provided ProHap databases

### Reproducing figures
Due to the computational demands, we expect the Snakemake pipeline to be run on a high-performance server. We provide three Jupyter notebooks to reproduce the figures using the result files as above.

To reproduce the plots in Figure 2, run the cells of the notebook in `notebooks/Populations_comparison_plots.ipynb`.
To reproduce the plots in Extended Data Figure 1 and 3, run the cells of the notebook in `notebooks/variant_types_and_overlap_plots.ipynb`.
To reproduce the plots in Extended Data Figure 4 and 5, run the cells of the notebook in `notebooks/spectra_mirror_plots.ipynb`.
