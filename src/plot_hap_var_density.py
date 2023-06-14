import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Creates a scatterplot of variant MAF (x axis) vs. the frequency of the most frequent haplotype where this variant is present")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input haplotype table")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output PNG file")

args = parser.parse_args()

print ("Reading", args.input_file)
haplo_df = pd.read_csv(args.input_file, sep='\t', header=0)

variants_data = {}  # dictionary of variants accessed by their ID (chr:pos:ref>alt)

variant_colors = {
    'SAP': '#C1A40D',
    'inframe-indel': '#CB00BE',
    'frameshift': '#980000',
    'synonym': '#AAAAAA'
}

def check_variant_class(ref, alt, prot_ref, prot_alt):
    if (prot_ref == prot_alt):
        return 'synonym'
    if (len(alt) == len(ref)):
        return 'SAP'
    if ((len(ref) - len(alt)) % 3) == 0:
        return 'inframe-indel'
    else:
        return 'frameshift'

for index,row in haplo_df.iterrows():

    varIDs = [ row['chromosome'] + ':' + varID for varID in row['DNA_changes'].split(';') ]
    varMAFs = [ int(MAF) for MAF in row['allele_frequencies'].split(';') ]
    prot_changes = [ [ch.split(':',1)[1].split('>',1)[0], ch.split(':',2)[2].split('(',1)[0]] for ch in row['all_protein_changes'].split(';') ]

    for var_idx, varID in enumerate(varIDs):

        var_class = check_variant_class(varID.split(':')[2].split('>')[0], varID.split('>')[1], prot_changes[var_idx][0], prot_changes[var_idx][1])

        if varID in variants_data:
            var_data = variants_data[varID]

            variants_data[varID][1] = max(var_data[1], row['frequency'])
            #if (var_data[2] == 'synonym') and (var_class != 'synonym'):
            #    variants_data[varID][2] == 'synonym/' + var_class
            #elif (var_class == 'synonym'):
            #    variants_data[varID][2] = 'synonym/' + var_data[2]
        else:
            variants_data[varID] = [varMAFs[var_idx], row['frequency'], var_class]

plot_x = [ var_data[0] for var_data in variants_data.values() ] 
plot_y = [ var_data[1] for var_data in variants_data.values() ]
plot_c = [ variant_colors[var_data[2]] for var_data in variants_data.values() ]

plt.scatter(plot_x, plot_y, c=plot_c, alpha=0.5)
plt.savefig(args.output_file, dpi=200)