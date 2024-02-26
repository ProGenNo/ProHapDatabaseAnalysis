from tqdm import tqdm
import argparse
from multiprocessing import Pool
import pandas as pd

def parse_snp_list(SNPs):
    result = [] 
    for haplotypeGroup in SNPs.split(';'):
        for haplotypeName in haplotypeGroup.split('|'):
            if haplotypeName == '-':
                result.append('-')
                continue
            proteinID, SNPs = haplotypeName.split(':')
            for SNP in SNPs.split(','):
                result.append(proteinID + ':' + SNP)
    
    return result

def get_protein_coverage(data):
    regions = []
    current_pos = 0
    region_start = 0
    current_freq = -1

    active_regions = []
    eventQ = [] if (data[0][0] == 0) else [{'pos': 0, 'event': 'start', 'freq': -1, 'pep_idx': -1}, {'pos': data[0][0], 'event': 'end', 'freq': -1, 'pep_idx': -1}]

    for idx, pep in enumerate(data):
        eventQ.append({ 'pos': pep[0], 'event': 'start', 'freq': pep[2], 'pep_idx': idx })
        eventQ.append({ 'pos': pep[0] + pep[1], 'event': 'end', 'freq': pep[2], 'pep_idx': idx })

    eventQ = sorted(eventQ, key=lambda x: x['pos'])

    for evt in eventQ:
        current_pos = evt['pos']

        if evt['event'] == 'start':
            if evt['freq'] > current_freq:
                if current_pos > region_start:
                    regions.append([region_start, current_pos, current_freq])
                current_freq = evt['freq']
                region_start = current_pos
            active_regions.append([evt['pep_idx'], evt['freq']])

        if evt['event'] == 'end':
            active_regions = list(filter(lambda x: x[0] != evt['pep_idx'], active_regions))
            new_freq = -1
            if (len(active_regions) > 0):
                new_freq = max(list(map(lambda x: x[1], active_regions)))
            if (new_freq < current_freq):
                if current_pos > region_start:
                    regions.append([region_start, current_pos, current_freq])
                region_start = current_pos
                current_freq = new_freq

    return regions

parser = argparse.ArgumentParser(
	description='Reads a peptide database, gives statistics on protein coverage by canonical/non-canonical peptides.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file, peptide database")

parser.add_argument("-g_id", dest="gene_ids", required=True,
                    help="csv file mapping transcript IDs to gene IDs")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="# threads to use")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file, header=0)
pep_df = pep_df[pep_df['pep_type1'].str.contains('variant')]

print ("Reading", args.gene_ids)
gene_id_df = pd.read_csv(args.gene_ids, header=0)

total_peptide_count = len(pep_df)

result_columns = ['TranscriptID', 'coverage']
result_data = []

# iterate through the peptide list - list is sorted by gene ID
current_gene = ''
 # accessed by protein ID, 'peptides' list of lists, 3 values: start, length, type

total_aa = [0, 0, 0, 0, 0]

all_peptides = {}

# split up the DF by gene ID
all_genes = gene_id_df['GeneID'].drop_duplicates().tolist()

def process_gene(geneIdx):
    geneID = all_genes[geneIdx]
    local_df = pep_df[(pep_df['matching_genes'] == pep_df['matching_genes']) & pep_df['matching_genes'].str.contains(geneID)]
    gene_transcripts = gene_id_df[gene_id_df['GeneID'] == geneID]['TranscriptID'].tolist()
    
    protein_peptides = {}       # dictionary of mapptings between peptides and proteins -> to be aggregated into coverage stats
    result = []

    # loop through peptides maping to this gene
    for index, row in local_df.iterrows():
        if geneID == '-' or row['pep_type1'] == 'contaminant':
            continue

        pep_length = len(row['sequence'])
        matching_transcripts = [ trID.split('.',1)[0] for trID in row['matching_transcripts'].split(';') ]
        relevant_transcript_idx = [ i for i,trID in enumerate(matching_transcripts) if trID in gene_transcripts ]

        # align this peptide to each matching protein (haplotype)
        for i in relevant_transcript_idx:
            trID = matching_transcripts[i]
            ref_stable_id = trID.split('.')[0]

            pep_type = row['pep_type1']
            
            if (row['covered_changes_peptide'] == row['covered_changes_peptide']):
                try:
                   prot_changes = row['covered_changes_protein'].split('|')[i]
                except:
                   # print(row['ID'], trID, row['pep_type1'], row['matching_transcripts'], row['covered_changes_protein'])
                   continue
                if (pep_type == 'single-variant') and (';' in prot_changes):
                    continue    # peptide was downgraded -> does not cover these variants reliably 

            pep_pos = int(row['positions_in_proteins'].split(';')[i]) + (int(row['preceding_indel_shift'].split(';')[i]) if row['preceding_indel_shift'] != '-' else 0)

            # add into a dictionary of peptide - ref. protein mappings
            if ref_stable_id in protein_peptides:
                protein_peptides[ref_stable_id]['peptides'].append([pep_pos, pep_length, row['frequency']])
            else:
                protein_peptides[ref_stable_id] = { 'peptides': [ [pep_pos, pep_length, row['frequency']] ] }

    # aggregate coverage by each protein
    for trID in protein_peptides:
        coverage = get_protein_coverage(sorted(protein_peptides[trID]['peptides'], key=lambda x: x[0]))
        coverage_str = ";".join(list(map(lambda x: "_".join(list(map(lambda y: str(y), x))), coverage)))

        result.append([trID, coverage_str])

    return result

print ('Annotating proteome coverage.')

with Pool(args.threads) as p:
    gene_results = list(tqdm(p.imap_unordered(process_gene, range(len(all_genes))), total=len(all_genes)))
    #gene_results = list(map(process_gene, range(len(all_genes))))
    p.close()
    p.join()

    for i, result in enumerate(gene_results):
        result_data.extend(result)

    print ('Done.')

    result_df = pd.DataFrame(columns=result_columns, data=result_data)
    result_df.to_csv(args.output_file, sep='\t', header=True, index=False)
