from common import read_fasta
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Reads a precomputed peptide coverage file, gives summary statistics on protein coverage by canonical/non-canonical peptides.')

parser.add_argument("-i", dest="input_files", required=True,
                    help="input TSV files (comma-separated list), peptide coverage by protein")

parser.add_argument("-ref_fa", dest="ref_fasta", required=True,
                    help="reference proteome (Ensembl) fasta file") 

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ("Reading", args.ref_fasta)
all_proteins = read_fasta(args.ref_fasta)

inputs = args.input_files.split(",")

results_concat = pd.concat([ pd.read_table(infile) for infile in inputs ])
results_concat = results_concat[results_concat['coverage_aa'] > 0]

total_aa_sum = 0

for protein in all_proteins.values():
    total_aa_sum += len(protein['sequence'].replace('*', ''))

results_concat['coverage_percent'] = results_concat['coverage_aa'].apply(lambda x: (x / total_aa_sum) * 100)
results_concat.to_csv(args.output_file, sep='\t', index=False)
