import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Annotate variant types inside the haplotype table")

parser.add_argument("-hap_tsv", dest="haplo_db", required=True,
                    help="haplotypes tab-separated file")

parser.add_argument("-oh", dest="output_file_haplo", required=True,
                    help="output TSV file")

parser.add_argument("-ov", dest="output_file_var", required=True,
                    help="output CSV file")

args = parser.parse_args()

print ("Reading", args.haplo_db)
haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)

variants_aggregated = {}

def get_var_types(row):
    var_types = []
    prot_changes = [ [ch.split(':',1)[1].split('>',1)[0], ch.split(':',2)[2] ] for ch in row['all_protein_changes'].split(';') ]

    for i,varLoc in enumerate(row['DNA_changes'].split(';')):
        varID = row['chromosome'] + ':' + varLoc
        REF = varLoc.split(':')[1].split('>')[0]
        ALT = varLoc.split('>')[1]

        varType = ""
        if ('+fs' in prot_changes[i][1]):
            varType = 'frameshift'
        elif (prot_changes[i][0] == prot_changes[i][1].split('(')[0]):
            varType = 'synonymous'
        elif (len(REF) == len(ALT)):
            varType = 'SAV'
        else:
            varType = 'inframe_indel'

        if '(fs)' in prot_changes[i][1]:
            varType += '_after_fs'

        var_types.append(varType)

        if (varID in variants_aggregated):
            if varType not in variants_aggregated[varID]:
                variants_aggregated[varID].append(varType)
        else:
            variants_aggregated[varID] = [varType]

    return ';'.join(var_types)

haplo_db['variant_types'] = haplo_db.apply(get_var_types, axis=1)
haplo_db.to_csv(args.output_file_haplo, sep='\t', index=False)

var_data = []

for varID in variants_aggregated:
    var_data.append([varID, '|'.join(variants_aggregated[varID])])

var_df = pd.DataFrame(columns=['variantID', 'possible_conseq'], data=var_data)
var_df.to_csv(args.output_file_var, index=False)
