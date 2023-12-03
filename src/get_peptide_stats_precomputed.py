from common import read_fasta
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Reads a precomputed peptide stats file, gives summary statistics on protein coverage by canonical/non-canonical peptides.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file, peptide stats by protein")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="input FASTA file")

args = parser.parse_args()

all_proteins = read_fasta(args.fasta_file)

pep_df = pd.read_csv(args.input_file, sep='\t', header=0)

total_aa = [0, 0, 0, 0, 0]

for index, row in pep_df.iterrows():
    coverage_str = row['coverage']

    for region_str in coverage_str.split(';'):
        region = list(map(lambda x: int(x), region_str.split('_')))
        region_len = region[1] - region[0]
        region_type = region[2] + 1
        total_aa[region_type] += region_len

total_aa_sum = 0

for protein in all_proteins.values():
    if ('ensref' in protein['tag']):
        total_aa_sum += len(protein['sequence'].replace('*', ''))

print ("Proteome length:", total_aa_sum)
print ("Canonical proteome: %d AAs - %.2f %%,\npossible single-variant peptides: %d AAs - %.2f %%,\npossible multi-variant peptides: %d AAs - %.2f %%,\npossible frameshift peptides: %d AAs - %.2f %%,\nsequences not matching to peptides: %d AAs - %.2f %%" % (total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100,total_aa[4], (total_aa[4] / total_aa_sum) * 100, (total_aa_sum - sum(total_aa[1:])), ((total_aa_sum - sum(total_aa[1:])) / total_aa_sum) * 100))
