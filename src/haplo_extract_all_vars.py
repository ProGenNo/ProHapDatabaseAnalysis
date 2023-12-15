import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Annotate variant types inside the haplotype table")

parser.add_argument("-hap_tsv", dest="haplo_db", required=True,
                    help="haplotypes tab-separated file")

parser.add_argument("-o", dest="output_file_var", required=True,
                    help="output CSV file")

args = parser.parse_args()

print ("Reading", args.haplo_db)
haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)

variants_aggregated = {}

for index,row in haplo_db.iterrows():
    var_types = row['variant_types'].split(';')
    var_ids = [ str(row['chromosome']) +':' + ch for ch in row['DNA_changes'].split(';')]

    for i,varID in enumerate(var_ids):
        if (varID in variants_aggregated):
            if var_types[i] not in variants_aggregated[varID]:
                variants_aggregated[varID].append(var_types[i])
        else:
            variants_aggregated[varID] = [var_types[i]]

var_data = []

for varID in variants_aggregated:
    var_data.append([varID, '|'.join(variants_aggregated[varID])])

var_df = pd.DataFrame(columns=['variantID', 'possible_conseq'], data=var_data)
var_df.to_csv(args.output_file_var, index=False)
