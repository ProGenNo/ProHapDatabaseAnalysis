from tqdm import tqdm
from common import get_protein_name_dict, read_fasta, get_protein_coverage
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

parser = argparse.ArgumentParser(
	description='Reads a peptide database, gives statistics on protein coverage by canonical/non-canonical peptides.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file, peptide database")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="input FASTA file")

parser.add_argument("-ref_fa", dest="ref_fasta", required=True,
                    help="reference proteome (Ensembl) fasta file") 

parser.add_argument("-g_id", dest="gene_ids", required=True,
                    help="csv file mapping transcript IDs to gene IDs")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="# threads to use")

parser.add_argument("-sep", dest="input_sep", required=False,
                    help="separator in input file (default: ',')", default=',')

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ('Reading', args.fasta_file)
name_dict, all_proteins = get_protein_name_dict(args.fasta_file)

print ("Reading", args.ref_fasta)
ref_proteins = read_fasta(args.ref_fasta)

print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file, header=0, sep=args.input_sep)

print ("Reading", args.gene_ids)
gene_id_df = pd.read_csv(args.gene_ids, header=0)

total_peptide_count = len(pep_df)

result_columns = ['TranscriptID', 'coverage']
result_data = []

# iterate through the peptide list - list is sorted by gene ID
current_gene = ''
 # accessed by protein ID, 'peptides' list of lists, 3 values: start, length, type

total_aa = [0, 0, 0, 0, 0]

# split up the DF by gene ID
all_genes = gene_id_df['GeneID'].drop_duplicates().tolist()

def process_gene(geneIdx):
    geneID = all_genes[geneIdx]
    local_df = pep_df[(pep_df['matching_genes'] == pep_df['matching_genes']) & pep_df['matching_genes'].str.contains(geneID)]
    gene_transcripts = gene_id_df[gene_id_df['GeneID'] == geneID]['TranscriptID'].tolist()
    
    protein_peptides = {}       # dictionary of mapptings between peptides and proteins -> to be aggregated into coverage stats
    local_aa = [0, 0, 0, 0, 0]     # length of covered regions by type: 0: should not be covered; 1: always canonical; 2: possibly single-variant; 3: possibly multi-variant
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

            encoded_type = -1
            if pep_type == 'canonical':
                encoded_type = 0
            elif pep_type == 'single-variant':
                encoded_type = 1
            elif pep_type == 'multi-variant':
                encoded_type = 2
            elif pep_type == 'frameshift':
                encoded_type = 3

            pep_pos = int(row['positions_in_proteins'].split(';')[i]) + (int(row['preceding_indel_shift'].split(';')[i]) if row['preceding_indel_shift'] != '-' else 0)

            # add into a dictionary of peptide - ref. protein mappings
            if ref_stable_id in protein_peptides:
                protein_peptides[ref_stable_id]['peptides'].append([pep_pos, pep_length, encoded_type])
            else:
                protein_peptides[ref_stable_id] = { 'peptides': [ [pep_pos, pep_length, encoded_type] ] }

    # aggregate coverage by each protein
    for trID in protein_peptides:
        coverage = get_protein_coverage(sorted(protein_peptides[trID]['peptides'], key=lambda x: x[0]))

        for region in coverage:
            region_len = region[1] - region[0]
            region_type = region[2] + 1
            local_aa[region_type] += region_len

        # add the remainder of the protein that's not covered, if any
        prot_len = len(ref_proteins[trID]['sequence'])
        if (len(coverage) > 0):
            local_aa[0] += (prot_len - coverage[-1][1])
            coverage.append([coverage[-1][1], prot_len, -1])
        else:
            coverage = [[0, prot_len, -1]]

        coverage_str = ";".join(list(map(lambda x: "_".join(list(map(lambda y: str(y), x))), coverage)))

        result.append([trID, coverage_str])

    return [ result, local_aa ]

print ('Annotating proteome coverage.')

with Pool(args.threads) as p:
    gene_results = list(tqdm(p.imap_unordered(process_gene, range(len(all_genes))), total=len(all_genes)))
#gene_results = list(map(process_gene, range(len(all_genes))))
    p.close()
    p.join()

for i, result in enumerate(gene_results):
    if (len(result[0]) == 0):
        continue 

    result_data.extend(result[0])

    total_aa[0] += result[1][0]
    total_aa[1] += result[1][1]
    total_aa[2] += result[1][2]
    total_aa[3] += result[1][3]
    total_aa[4] += result[1][4]

print ('Done.')

result_df = pd.DataFrame(columns=result_columns, data=result_data)
result_df.to_csv(args.output_file, sep='\t', header=True, index=False)

# check the length of the proteone (= sum of lengths of all canonical proteins in fasta)
total_aa_sum = 0

for protein in ref_proteins.values():
    total_aa_sum += len(protein['sequence'].replace('*', ''))

print ("Proteome length:", total_aa_sum)
print ("Canonical proteome: %d AAs - %.2f %%,\npossible single-variant peptides: %d AAs - %.2f %%,\npossible multi-variant peptides: %d AAs - %.2f %%,\npossible frameshift peptides: %d AAs - %.2f %%,\nsequences not matching to peptides: %d AAs - %.2f %%" % (total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100,total_aa[4], (total_aa[4] / total_aa_sum) * 100, (total_aa_sum - sum(total_aa[1:])), ((total_aa_sum - sum(total_aa[1:])) / total_aa_sum) * 100))
