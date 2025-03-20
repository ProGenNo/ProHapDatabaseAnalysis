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

def get_conflicting_peptides(seq_corrected, seq_orig, df):
    return ';'.join(df[(df['Sequence_I_to_L'] == seq_corrected) & (df['Sequence'] != seq_orig)]['Sequence'].tolist())

print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file, sep='\t')
prohap_il_peptides = pep_df[pep_df['Sequence'].str.contains('I') | pep_df['Sequence'].str.contains('L')]
prohap_il_peptides['ID'] = prohap_il_peptides['ID'].apply(lambda x: 'prohap_'+x)
prohap_il_peptides['Database'] = 'ProHap_1kGP_ALL'
prohap_il_peptides['Sequence_I_to_L'] = prohap_il_peptides['Sequence'].apply(lambda seq: seq.replace('I', 'L'))
prohap_il_peptides['Ambiguous'] = prohap_il_peptides.apply(lambda row: get_conflicting_peptides(row['Sequence_I_to_L'], row['Sequence'], prohap_il_peptides))

print ('Reading', args.input_file_sp)
sp_pep_df = pd.read_table(args.input_file_sp)
sp_il_peptides = sp_pep_df[sp_pep_df['Sequence'].str.contains('I') | sp_pep_df['Sequence'].str.contains('L')]
sp_il_peptides['ID'] = sp_il_peptides['ID'].apply(lambda x: 'swissprot_'+x)
sp_il_peptides['Database'] = 'SwissProt'
sp_il_peptides['Sequence_I_to_L'] = sp_il_peptides['Sequence'].apply(lambda seq: seq.replace('I', 'L'))
sp_il_peptides['Ambiguous'] = sp_il_peptides.apply(lambda row: get_conflicting_peptides(row['Sequence_I_to_L'], row['Sequence'], sp_il_peptides))

print ('Reading', args.input_file_iso)
iso_pep_df = pd.read_table(args.input_file_iso)
iso_il_peptides = iso_pep_df[iso_pep_df['Sequence'].str.contains('I') | iso_pep_df['Sequence'].str.contains('L')]
iso_il_peptides['ID'] = iso_il_peptides['ID'].apply(lambda x: 'uniprot_'+x)
iso_il_peptides['Database'] = 'UniProt_isoform'
iso_il_peptides['Sequence_I_to_L'] = iso_il_peptides['Sequence'].apply(lambda seq: seq.replace('I', 'L'))
iso_il_peptides['Ambiguous'] = iso_il_peptides.apply(lambda row: get_conflicting_peptides(row['Sequence_I_to_L'], row['Sequence'], iso_il_peptides))

pd.concat([prohap_il_peptides, sp_il_peptides, iso_il_peptides]).to_csv(args.output_file, header=True, index=False)