import argparse
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser(description="Print out stats on variant types")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input: list of all peptides (tab-separated)")

parser.add_argument("-hap_tsv", dest="haplo_db", required=True,
                    help="haplotypes tab-separated file")

parser.add_argument("-pop", dest="population_code", required=True,
                    help="population code")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file", metavar="FILE")

args = parser.parse_args()

pep_df = pd.read_csv(args.input_file)
haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)
haplo_db.set_index('HaplotypeID', inplace=True)

def get_max_freq(prot_ids):
    if ('ENST' in prot_ids):
        return -1

    haplo_ids = [ protID for protID in prot_ids.split(';') if protID.startswith('haplo') ]
    if (len(haplo_ids) == 0):
        return -1
    
    all_frequencies = {}

    for haploID in haplo_ids:
        haplotype = haplo_db.loc[haploID]
        trID = haplotype['TranscriptID']
        if trID in all_frequencies:
            all_frequencies[trID].append(haplotype['frequency'])
        else:
            all_frequencies[trID] = [haplotype['frequency']]

    return max([ sum(freqs) for freqs in all_frequencies.values() ])

with Pool(args.threads) as p:
    frequencies = list(tqdm(p.imap(get_max_freq, pep_df['matching_proteins'].tolist()), total=len(pep_df)))

    pep_df['frequency'] = frequencies
    pep_df.to_csv(args.output_file, header=True, index=False)
