import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Computes the overlap of peptides between a given database and SwissProt and TrEMBL')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input CSV file, peptide database")

parser.add_argument("-sp", dest="input_file_sp", required=True,
                    help="input TSV file, SwissProt peptide database")

#parser.add_argument("-iso", dest="input_file_iso", required=True,
#                    help="input TSV file, SwissProt + TrEMBL + unreviewed isoforms peptide database")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

def concat_string_series(series):
    return ';'.join(series.drop_duplicates().tolist())

# ProHap peptides
print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file, sep='\t')

# Select only peptides containing I or L
prohap_il_peptides = pep_df[pep_df['Sequence'].str.contains('I') | pep_df['Sequence'].str.contains('L')].copy()
prohap_il_peptides['Sequence_I_to_L'] = prohap_il_peptides['Sequence'].apply(lambda seq: seq.replace('I', 'L'))

# Group together the peptides if I == L, and drop those that don't have an ambiguity
prohap_pep_ambiguous = prohap_il_peptides.groupby('Sequence_I_to_L').agg({'Sequence': concat_string_series, 'Proteins': concat_string_series}).reset_index()
prohap_pep_ambiguous = prohap_pep_ambiguous[prohap_pep_ambiguous['Sequence'].str.contains(';')]

# SwissProt peptides
print ('Reading', args.input_file_sp)
sp_pep_df = pd.read_table(args.input_file_sp)
sp_il_peptides = sp_pep_df[sp_pep_df['Sequence'].str.contains('I') | sp_pep_df['Sequence'].str.contains('L')].copy()
sp_il_peptides['Sequence_I_to_L'] = sp_il_peptides['Sequence'].apply(lambda seq: seq.replace('I', 'L'))

sp_pep_ambiguous = sp_il_peptides.groupby('Sequence_I_to_L').agg({'Sequence': concat_string_series, 'Proteins': concat_string_series}).reset_index()
sp_pep_ambiguous = sp_pep_ambiguous[sp_pep_ambiguous['Sequence'].str.contains(';')]

'''
# UniProt Isoform peptides
print ('Reading', args.input_file_iso)
iso_pep_df = pd.read_table(args.input_file_iso)
iso_il_peptides = iso_pep_df[iso_pep_df['Sequence'].str.contains('I') | iso_pep_df['Sequence'].str.contains('L')].copy()
iso_il_peptides['Sequence_I_to_L'] = iso_il_peptides['Sequence'].apply(lambda seq: seq.replace('I', 'L'))

iso_pep_ambiguous = iso_il_peptides.groupby('Sequence_I_to_L').agg({'Sequence': concat_string_series, 'Proteins': concat_string_series}).reset_index()
iso_pep_ambiguous = iso_pep_ambiguous[iso_pep_ambiguous['Sequence'].str.contains(';')]
'''

df_joined = prohap_pep_ambiguous.join(sp_pep_ambiguous.set_index('Sequence_I_to_L'), on='Sequence_I_to_L', lsuffix='_ProHap', rsuffix='_SwissProt', how='outer')
df_joined.fillna('-')
df_joined.to_csv(args.output_file, sep='\t', header=True, index=False)