from tqdm import tqdm
from common import get_protein_coverage
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

parser.add_argument("-g_id", dest="gene_ids", required=True,
                    help="csv file mapping transcript IDs to gene IDs")

parser.add_argument("-hap_prefix", dest="haplo_prefix", required=False,
                    help="prefix for haplotype protein ID (default: 'haplo_')", default='haplo_')

parser.add_argument("-hap_tsv", dest="haplo_db", required=False,
                    help="haplotypes tab-separated file (optional)", default=None)

parser.add_argument("-pop", dest="population_code", required=True,
                    help="population code")

parser.add_argument("-s", dest="samples_filename", required=True,
                    help="tab-separated file with sample information (must include 'Sample name' and 'Sex' columns)")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="# threads to use")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file, header=0)
pep_df = pep_df[pep_df['pep_type1'].str.contains('variant') | pep_df['pep_type1'].str.contains('frameshift')]

print ('Processing', len(pep_df), 'peptides.')

print ("Reading", args.haplo_db)
haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)
haplo_db.set_index('HaplotypeID', inplace=True)

print ("Reading", args.gene_ids)
gene_id_df = pd.read_csv(args.gene_ids, header=0)

samples_df = pd.read_csv(args.samples_filename, sep='\t')
sample_ids = [ sID for sID in samples_df[samples_df['Superpopulation code'] == args.population_code]['Sample name'].tolist() ]

total_peptide_count = len(pep_df)

result_columns = ['SampleID', 'TranscriptID', 'coverage_aa']
result_data = []

# iterate through the peptide list - list is sorted by gene ID

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

            pep_samples = []
            for protID in row['matching_proteins'].split(';'):
                if protID.startswith(args.haplo_prefix):
                    pep_samples.extend([ sample.split(':',1)[0] for sample in haplo_db.loc[protID]['samples'].split(';') ])
            pep_samples = list(dict.fromkeys(pep_samples))
            #print(row['matching_proteins'])

            # add into a dictionary of peptide - ref. protein mappings
            if ref_stable_id in protein_peptides:
                for sID in pep_samples:
                    if (sID in protein_peptides[ref_stable_id]):
                        protein_peptides[ref_stable_id][sID].append([pep_pos, pep_length, encoded_type])
                    else:
                        protein_peptides[ref_stable_id][sID] = [ [pep_pos, pep_length, encoded_type] ] 
            else:
                protein_peptides[ref_stable_id] = {}
                for sID in pep_samples:
                    protein_peptides[ref_stable_id][sID] = [ [pep_pos, pep_length, encoded_type] ] 

    # aggregate coverage by each protein
    for trID in protein_peptides:
        for sID in protein_peptides[trID]:
            coverage = get_protein_coverage(sorted(protein_peptides[trID][sID], key=lambda x: x[0]))
            coverage_aa = sum([ segment[1] - segment[0] for segment in coverage if (segment[2] > 0) ])
            #coverage_str = ";".join(list(map(lambda x: "_".join(list(map(lambda y: str(y), x))), coverage)))

            result.append([sID, trID, coverage_aa])

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
    result_df.to_csv('.'.join(args.output_file.split('.')[:-1]) + '_per_transcript.tsv', sep='\t', header=True, index=False)

    # sum up AAs across transcript for each sample
    summary_data = [ [sID, sum(result_df[result_df['SampleID'] == sID]['coverage_aa'].tolist())] for sID in sample_ids ]
    summary_df = pd.DataFrame(data=summary_data, columns=['sampleID', 'coverage_aa'])

    summary_df.to_csv(args.output_file, sep='\t', header=True, index=False)
