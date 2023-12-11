from common import read_fasta
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Reads a precomputed peptide coverage file, gives summary statistics on protein coverage by canonical/non-canonical peptides.')

parser.add_argument("-i", dest="input_files", required=True,
                    help="input TSV files (comma-separated list), peptide coverage by protein")

parser.add_argument("-pop", dest="populations", required=True,
                    help="population codes (comma-separated list), need to match the list of input files")

parser.add_argument("-ref_fa", dest="ref_fasta", required=True,
                    help="reference proteome (Ensembl) fasta file") 

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ("Reading", args.ref_fasta)
all_proteins = read_fasta(args.ref_fasta)

inputs = args.input_filenames.split(",")
populations = args.populations.split(',')

if (len(inputs) != len(populations)):
    raise Exception("Please provide the same number of populations as input files. Currently provided input files: " + str(inputs) + "; currently provided population codes: " + str(populations))

total_aa_sum = 0

for protein in all_proteins.values():
    total_aa_sum += len(protein['sequence'].replace('*', ''))

print ("Proteome length:", total_aa_sum)

aggregated_results = []

for i,pop in enumerate(populations):
    print ("Reading", inputs[i])
    pep_df = pd.read_csv(inputs[i], sep='\t', header=0)

    total_aa = [0, 0, 0, 0, 0]

    for index, row in pep_df.iterrows():
        coverage_str = row['coverage']

        for region_str in coverage_str.split(';'):
            region = list(map(lambda x: int(x), region_str.split('_')))
            region_len = region[1] - region[0]
            region_type = region[2] + 1
            total_aa[region_type] += region_len

    aggregated_results.append([pop, total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100,total_aa[4], (total_aa[4] / total_aa_sum) * 100, (total_aa_sum - sum(total_aa[1:])), ((total_aa_sum - sum(total_aa[1:])) / total_aa_sum) * 100])

    print ('Population:', pop)
    print ("Canonical proteome: %d AAs - %.2f %%,\npossible single-variant peptides: %d AAs - %.2f %%,\npossible multi-variant peptides: %d AAs - %.2f %%,\npossible frameshift peptides: %d AAs - %.2f %%,\nsequences not matching to peptides: %d AAs - %.2f %%" % (total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100,total_aa[4], (total_aa[4] / total_aa_sum) * 100, (total_aa_sum - sum(total_aa[1:])), ((total_aa_sum - sum(total_aa[1:])) / total_aa_sum) * 100))
    print ('-------------------------------')

results_df = pd.DataFrame(data=aggregated_results, columns=['population', 'canonical_aa', 'canonical_perc', 'singlevar_aa', 'singlevar_perc', 'multivar_aa', 'multivar_perc', 'frameshift_aa', 'frameshift_perc', 'not_covered_aa', 'not_covered_perc'])
results_df.to_csv(args.output_file, sep='\t', index=False)
