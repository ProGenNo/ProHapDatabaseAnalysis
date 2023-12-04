import argparse
import pandas as pd
from common import read_fasta
from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser(description="Remove haplotypes not present in a given population from the database")

parser.add_argument("-f", dest="input_fasta", required=True,
                    help="input FASTA file")

parser.add_argument("-hap_tsv", dest="haplo_db", required=True,
                    help="haplotypes tab-separated file")

parser.add_argument("-pop", dest="population_code", required=True,
                    help="population code")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")

parser.add_argument("-output_tsv", dest="output_file", required=True,
                    help="output TSV file")

parser.add_argument("-output_fasta", dest="output_fasta", required=True,
                    help="output FASTA file")

args = parser.parse_args()

print ("Reading", args.haplo_db)
haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)

all_proteins = read_fasta(args.input_fasta)

def check_population(freqs):
    if args.population_code in freqs:
        freq = float(freqs.split(args.population_code, 1)[1].split(';',1)[0][1:])
        if freq > 0.0:
            return True
    return False

def check_protein(prot):
    remaining_proteins = []
    for elem in prot['matching_proteins']:
        kept_IDs = [ idx for idx,protID in enumerate(elem) if ((not protID.startswith('haplo')) or (protID in retained_haplotypes)) ]
        remaining_proteins.append(kept_IDs)
    return remaining_proteins
    
haplo_db['keep'] = haplo_db['frequency_superpopulation'].apply(check_population)
haplo_db = haplo_db[haplo_db['keep']]
haplo_db.drop('keep', axis=1, inplace=True)
haplo_db.to_csv(args.output_file, sep='\t', index=False)

retained_haplotypes = haplo_db['HaplotypeID'].tolist()

with Pool(args.threads) as p:
    filtered_protIDs = list(tqdm(p.imap(check_protein, all_proteins.values()), total=len(all_proteins)))

    outfile = open(args.output_fasta, 'w')
    for i,prot in enumerate(all_proteins.values()):
        remaining_proteins = filtered_protIDs[i]
        if any([ len(kept_IDs) > 0 for kept_IDs in remaining_proteins ]):
            matching_prots = []
            rfs = []
            pos = []
            start = []

            for i,kept_IDs in enumerate(remaining_proteins):
                if (len(kept_IDs) > 0):
                    matching_prots.append(','.join([ prot['matching_proteins'][i][j] for j in kept_IDs]))
                    rfs.append(','.join([ str(prot['reading_frames'][i][j]) for j in kept_IDs]))
                    pos.append(str(prot['seq_positions'][i]))
                    start.append(prot['description'].split('start:',1)[1].split(';')[i].split(maxsplit=1)[0])

            description = 'position_within_protein:' + ';'.join(pos) + ' start:' + ';'.join(start) + ' matching_proteins:' + ';'.join(matching_prots) + ' reading_frame:' + ';'.join(rfs)

            outfile.write('>' + prot['tag'] + '|' + prot['accession'] + '|' + description + '\n')
            outfile.write(prot['sequence'] + '\n')
    outfile.close()
