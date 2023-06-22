import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Print out stats on variant types")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input: list of all included variants")

args = parser.parse_args()

var_df = pd.read_csv(args.input_file)

var_df['is_SAV'] = var_df['possible_conseq'].apply(lambda x: all([ 'SAV' in c for c in x.split('|') ]))
var_df['is_synon'] = var_df['possible_conseq'].apply(lambda x: all([ 'synonymous' in c for c in x.split('|') ]))
var_df['is_stop_lost'] = var_df['possible_conseq'].apply(lambda x: all([ 'stop_lost' in c for c in x.split('|') ]))
var_df['is_after_fs'] = var_df['possible_conseq'].apply(lambda x: all([ 'after_fs' in c for c in x.split('|') ]))
var_df['affects_splice_site'] = var_df['possible_conseq'].apply(lambda x: all([ 'splice_variant' in c for c in x.split('|') ]))

print ('Number of SAVs:', len(var_df[var_df['is_SAV']]))
print ('Number of synonymous variants:', len(var_df[var_df['is_synon']]))
print ('Number of stop-lost variants:', len(var_df[var_df['is_stop_lost']]))
print ('Number of variants certainly affecting a splice site:', len(var_df[var_df['affects_splice_site']]))
print ('Number of variants possibly affecting a splice site:', len(var_df[var_df['possible_conseq'].str.contains('splice_variant')]))
print ('Number of ambiguities (synonymous or SAV):', len(var_df[var_df['possible_conseq'].str.contains('SAV') & var_df['possible_conseq'].str.contains('synonymous')]))
print ('Number of ambiguities (stop-lost or SAV):', len(var_df[var_df['possible_conseq'].str.contains('SAV') & var_df['possible_conseq'].str.contains('stop_lost')]))
print ('Number of ambiguities (stop-lost or synonymous):', len(var_df[var_df['possible_conseq'].str.contains('synonymous') & var_df['possible_conseq'].str.contains('stop_lost')]))
print ('Number of in-frame indels:', len(var_df[var_df['possible_conseq'].str.contains('inframe_indel')]))
print ('Number of frameshifts:', len(var_df[var_df['possible_conseq'].str.contains('frameshift')]))
print ('Number of variants occurring only after a frameshift:', len(var_df[var_df['is_after_fs']]))
print ('Number of variants occurring possibly after a frameshift:', len(var_df[var_df['possible_conseq'].str.contains('after_fs')]))

