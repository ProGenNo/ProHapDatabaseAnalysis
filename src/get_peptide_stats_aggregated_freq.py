from common import read_fasta
from math import log
import argparse
import bisect
import pandas as pd
import numpy as np

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

inputs = args.input_files.split(",")
populations = args.populations.split(',')

if (len(inputs) != len(populations)):
    raise Exception("Please provide the same number of populations as input files. Currently provided input files: " + str(inputs) + "; currently provided population codes: " + str(populations))

total_aa_sum = 0

for protein in all_proteins.values():
    total_aa_sum += len(protein['sequence'].replace('*', ''))

print ("Proteome length:", total_aa_sum)

aggregated_results = []
FREQ_STEPS = [ x for x in np.logspace(0, 1, 1000) ]
STEP_SUM = sum(FREQ_STEPS)
FREQ_POINTS = [0] + list(np.cumsum([ x / STEP_SUM for x in FREQ_STEPS ]))

for i,pop in enumerate(populations):
    print ("Reading", inputs[i])
    pep_df = pd.read_csv(inputs[i], sep='\t', header=0)

    total_aa = [ 0 for freq in FREQ_POINTS ]

    for index, row in pep_df.iterrows():
        coverage_str = row['coverage']

        for region_str in coverage_str.split(';'):
            region = list(map(lambda x: float(x), region_str.split('_')))
            region_len = region[1] - region[0]
            region_freq = region[2]
            
            if (region_freq >= 0):
                idx = bisect.bisect_right(FREQ_POINTS, region_freq)
                for i in range(idx):
                    total_aa[i] += region_len

    aggregated_results.append([ pop, ';'.join([ str(freq) for freq in FREQ_POINTS ]), ';'.join([ str(aa_len) for aa_len in total_aa ]), ';'.join([ str((aa_len / total_aa_sum) * 100) for aa_len in total_aa ]) ])
    
results_df = pd.DataFrame(data=aggregated_results, columns=['population', 'freq_trheshold', 'var_region_length', 'var_region_percent'])
results_df.to_csv(args.output_file, sep='\t', index=False)
