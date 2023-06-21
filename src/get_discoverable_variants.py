import argparse
import re
import pandas as pd

parser = argparse.ArgumentParser(description="Lists the discoverable variants, and the frequency of the most frequent haplotype where this variant is present")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input annotated peptides")

parser.add_argument("-hap_tsv", dest="haplo_db", required=True,
                    help="haplotypes tab-separated file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file")

args = parser.parse_args()

print ("Reading", args.input_file)
pep_df = pd.read_csv(args.input_file, header=0, low_memory=False)

print ("Reading", args.haplo_db)
haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)
haplo_db.set_index('HaplotypeID', inplace=True)

variants_data = {}  # dictionary of variants accessed by their ID (chr:pos:ref>alt)

for index,row in pep_df.iterrows():
    if (row['covered_alleles_dna'] != row['covered_alleles_dna']):
        continue

    varIDs = [ varID for varID in re.split(r"[,;\|]", row['covered_alleles_dna']) if ('>' in varID) ]
    hapIDs = [ protID for protID in row['matching_proteins'].split(';') if protID.startswith('haplo') ]

    max_hap_freq = 0

    for hapID in hapIDs:
        max_hap_freq = max(max_hap_freq, haplo_db.loc[hapID]['frequency'])

    for var_idx, varID in enumerate(varIDs):

        if varID in variants_data:
            var_data = variants_data[varID]
            pep_type = []
            enzymes = []

            if row['pep_type1'] not in var_data[1]:
                pep_type = var_data[1] + [row['pep_type1']]
            else :
                pep_type = var_data[1]

            if row['enzyme'] not in var_data[2]:
                enzymes = var_data[2] + [row['enzyme']]
            else:
                enzymes = var_data[2]

            variants_data[varID] = [max(var_data[0], max_hap_freq), pep_type, enzymes]
        else:
            variants_data[varID] = [max_hap_freq, [row['pep_type1']], [row['enzyme']]]

result_df_data = []

for varID in variants_data:
    result_df_data.append([ varID, variants_data[varID][0], '|'.join(variants_data[varID][1]), '|'.join(variants_data[varID][2]) ])

result_df = pd.DataFrame(columns=['variantID', 'max_haplotype_frequency', 'pep_type', 'enzyme'], data=result_df_data)
result_df.to_csv(args.output_file, sep='\t', index=False, header=True)
