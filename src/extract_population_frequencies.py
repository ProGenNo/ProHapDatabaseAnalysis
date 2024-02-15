import argparse
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser(description="Creates a table of haplotype frequencies per transcript per population / superpopulation")

parser.add_argument("-i", dest="input_filename", required=True,
                    help="input files, replace the superpopulation abbreviation with \"{pop}\"")

parser.add_argument("-transcripts", dest="transcript_list", required=True,
                    help="list of transcript IDs, provided in a CSV file", metavar="FILE")

parser.add_argument("-s", dest="samples_filename", required=True,
                    help="tab-separated file with sample information (must include 'Sample name' and 'Population code' and 'Superpopulation code' columns)")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")

parser.add_argument("-o_superpop", dest="output_file_superpop", required=True,
                    help="output TSV file for superpopulations")

args = parser.parse_args()

dataframes = {}

#haplo_df = pd.read_table(args.input_filename)

print ('Reading', args.samples_filename)
samples_df = pd.read_csv(args.samples_filename, sep='\t')

all_transcripts = pd.read_csv(args.transcript_list)['TranscriptID'].apply(lambda trID: trID.split('.',1)[0]).tolist()

#pop_codes = samples_df['Population code'].drop_duplicates().tolist()
superpop_codes = samples_df['Superpopulation code'].drop_duplicates().tolist()

for pop_code in superpop_codes:
    filename = args.input_filename.replace('{pop}', pop_code)
    print ('Reading', filename)
    dataframes[pop_code] = pd.read_table(filename)

def process_transcript(trID):
    superpop_freqs = {}

    #df_tr = haplo_df[haplo_df['TranscriptID'] == trID]

    for pop_code in superpop_codes:
        df_pop = dataframes[pop_code]
        df_pop_tr = df_pop[df_pop['TranscriptID'] == trID]
        superpop_freqs[pop_code] = df_pop_tr['frequency'].tolist()
        # df_pop_tr['frequency'].apply(lambda x: [float(superpop.split(':',1)[1]) for superpop in x.split(';') if (pop_code in superpop)][0] if (pop_code in x) else 0.0).tolist()

    return [trID, superpop_freqs]

with Pool(args.threads) as p:
    transcripts_data = list(tqdm(p.imap_unordered(process_transcript, all_transcripts), total=len(all_transcripts)))
    superpopulations_data = []

    for row in transcripts_data:

        spop_row = [row[0]]
        for spop_code in superpop_codes:
            freqs = row[1][spop_code]
            freqs.sort()
            spop_row.append(';'.join([ str(f) for f in freqs]))

        superpopulations_data.append(spop_row)

#pop_df = pd.DataFrame(data=populations_data, columns=['TranscriptID'] + pop_codes)
#print ('Writing file:', args.output_file_pop)
#pop_df.to_csv(args.output_file_pop, sep='\t', index=False)

spop_df = pd.DataFrame(data=superpopulations_data, columns=['TranscriptID'] + superpop_codes)
print ('Writing file:', args.output_file_superpop)
spop_df.to_csv(args.output_file_superpop, sep='\t', index=False)
