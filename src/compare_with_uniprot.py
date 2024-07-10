import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Computes the overlap of peptides between a given database and SwissProt and TrEMBL')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input CSV file, peptide database")

parser.add_argument("-sp", dest="input_file_sp", required=True,
                    help="input TSV file, SwissProt peptide database")

parser.add_argument("-iso", dest="input_file_iso", required=True,
                    help="input TSV file, SwissProt + TrEMBL + unreviewed isoforms peptide database")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file)
pep_df = pep_df[pep_df['enzyme'] == 'Trypsin'][['ID', 'sequence', 'pep_type1']]
pep_df.rename(columns={'sequence': 'Sequence'}, inplace=True)

print ('Reading', args.input_file_sp)
sp_pep_df = pd.read_table(args.input_file_sp)

print ('Reading', args.input_file_iso)
iso_pep_df = pd.read_table(args.input_file_iso)

print ('Computing overlap')
df_joined = pep_df.join(sp_pep_df[['ID', 'Sequence']].set_index('Sequence'), on='Sequence', how='outer', rsuffix='_sp')
df_joined = df_joined.join(iso_pep_df[['ID', 'Sequence']].set_index('Sequence'), on='Sequence', how='outer', rsuffix='_isoform')

print ('Done')
df_joined.to_csv(args.output_file, sep='\t', index=False)