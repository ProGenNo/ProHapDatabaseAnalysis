from common import read_fasta
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Reads a precomputed peptide coverage file, gives summary statistics on protein coverage by canonical/non-canonical peptides.')

parser.add_argument("-i", dest="input_files", required=True,
                    help="input TSV files (comma-separated list), peptide coverage by protein")

parser.add_argument("-pop", dest="populations", required=True,
                    help="population codes (comma-separated list), need to match the list of input files")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

inputs = args.input_files.split(",")
populations = args.populations.split(',')

if (len(inputs) != len(populations)):
    raise Exception("Please provide the same number of populations as input files. Currently provided input files: " + str(inputs) + "; currently provided population codes: " + str(populations))

aggregated_results = []

for i,pop in enumerate(populations):
    print ("Reading", inputs[i])
    pep_df = pd.read_csv(inputs[i], sep='\t', header=0).fillna('-')

    aggregated_results.append([
        pop, 
        len(pep_df[pep_df['ID'] != '-']), 
        len(pep_df[(pep_df['ID'] != '-') & (pep_df['ID_sp'] == '-') & (pep_df['ID_isoform'] == '-')]), 
        len(pep_df[(pep_df['ID'] != '-') & (pep_df['ID_sp'] == '-') & (pep_df['ID_isoform'] != '-')]), 
        len(pep_df[(pep_df['ID'] != '-') & (pep_df['ID_sp'] != '-') & (pep_df['ID_isoform'] != '-')]), 
        len(pep_df[(pep_df['ID'] == '-') & (pep_df['ID_sp'] != '-') & (pep_df['ID_isoform'] != '-')]), 
        len(pep_df[(pep_df['ID'] == '-') & (pep_df['ID_sp'] == '-') & (pep_df['ID_isoform'] != '-')])
    ])

results_df = pd.DataFrame(
    data=aggregated_results, 
    columns= [
        'population', 
        'total_tryptic_peptides', 
        'only_prohap', 
        'overlap_prohap_isoform', 
        'overlap_prohap_isoform_swissprot', 
        'overlap_isoform_swissprot', 
        'only_isoform'
    ]
) 

results_df.to_csv(args.output_file, sep='\t', index=False)
