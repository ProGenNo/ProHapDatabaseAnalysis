import argparse
import bisect
import pandas as pd

parser = argparse.ArgumentParser(
	description='Prints variant discoverability and stats by consequence type.')

parser.add_argument("-i_all", dest="input_file_all", required=True,
                    help="input CSV file")

parser.add_argument("-i_disc", dest="input_file_discoverable", required=True,
                    help="input CSV file")

args = parser.parse_args()

df_all = pd.read_csv(args.input_file_all, header=0)
df_disc = pd.read_table(args.input_file_discoverable, header=0)
disc_ids = df_disc['variantID'].tolist()
disc_ids.sort()

df_all['is_discoverable'] = df_all['variantID'].apply(lambda varID: disc_ids[min(bisect.bisect_left(disc_ids, varID), len(disc_ids) - 1)] == varID)

df_all['is_synonym'] = df_all['possible_conseq'].apply(lambda x: all(['synonymous' in c for c in x.split('|') if c != 'splice_variant']) and (x != 'splice_variant'))
df_all['is_SAV'] = df_all['possible_conseq'].apply(lambda x: all(['SAV' in c for c in x.split('|') if c != 'splice_variant']) and (x != 'splice_variant'))
df_all['is_stop_lost'] = df_all['possible_conseq'].apply(lambda x: all(['stop_lost' in c for c in x.split('|') if c != 'splice_variant']) and (x != 'splice_variant'))
df_all['is_indel'] = df_all['possible_conseq'].apply(lambda x: all(['inframe_indel' in c for c in x.split('|') if c != 'splice_variant']) and (x != 'splice_variant'))
df_all['is_frameshift'] = df_all['possible_conseq'].apply(lambda x: all(['frameshift' in c for c in x.split('|') if c != 'splice_variant']) and (x != 'splice_variant'))
df_all['is_splice_site'] = df_all['possible_conseq'].apply(lambda x: all(['splice' in c for c in x.split('|')]))

df_all['is_SAV_synon'] = df_all['possible_conseq'].apply(lambda x: ('SAV' in x) and ('synonymous' in x))
df_all['is_SAV_stoplost'] = df_all['possible_conseq'].apply(lambda x: ('SAV' in x) and ('stop_lost' in x) and ('synonymous' not in x))
df_all['is_stoplost_synon'] = df_all['possible_conseq'].apply(lambda x: ('stop_lost' in x) and ('synonymous' in x) and ('SAV' not in x))

df_all['possible_after_fs'] = df_all['possible_conseq'].apply(lambda x: any([('after_fs' in c) for c in x.split('|') if ('synonymous' not in c) and ('frameshift' not in c)]) and ('after_fs' in x))
df_all['always_after_fs'] = df_all['possible_conseq'].apply(lambda x: all([('after_fs' in c) and ('synonymous' not in c) and ('frameshift' not in c) for c in x.split('|')]) and ('after_fs' in x))

print("------------- All included variants -------------")
print("Non-synonymous variants:", len(df_all[~df_all['is_synonym']]))
print("SAVs:", len(df_all[df_all['is_SAV']]))
print("Stop lost:", len(df_all[df_all['is_stop_lost']]))
print("Inframe indels:", len(df_all[df_all['is_indel']]))
print("Frameshifts:", len(df_all[df_all['is_frameshift']]))
print("Affecting a splice site", len(df_all[df_all['is_splice_site']]))
print("Synonymous", len(df_all[df_all['is_synonym']]))

print("Context-dependent, SAV or synonymous:", len(df_all[df_all['is_SAV_synon']]))
print("Context-dependent, SAV or stop lost:", len(df_all[df_all['is_SAV_stoplost']]))
print("Context-dependent, stop lost or synonymous:", len(df_all[df_all['is_stoplost_synon']]))

print("Possibly after frameshift:", len(df_all[df_all['possible_after_fs']]))
print("Always after frameshift:", len(df_all[df_all['always_after_fs']]))

print("------------- Discoverable variants -------------")
print("Total discoverable variants:", len(df_all[df_all['is_discoverable']]))
print("Discoverable SAV: %d (%.2f %%)" % (len(df_all[df_all['is_discoverable'] & (df_all['is_SAV'] | df_all['is_SAV_synon'])]), len(df_all[df_all['is_discoverable'] & (df_all['is_SAV'] | df_all['is_SAV_synon'])]) / len(df_all[(df_all['is_SAV'] | df_all['is_SAV_synon'])]) * 100))
print("Discoverable stop lost: %d (%.2f %%)" % (len(df_all[df_all['is_discoverable'] & (df_all['is_stop_lost'] | df_all['is_stoplost_synon'])]), len(df_all[df_all['is_discoverable'] & (df_all['is_stop_lost'] | df_all['is_stoplost_synon'])]) / len(df_all[(df_all['is_stop_lost'] | df_all['is_stoplost_synon'])]) * 100))
print("Discoverable indel: %d (%.2f %%)" % (len(df_all[df_all['is_discoverable'] & df_all['is_indel']]), len(df_all[df_all['is_discoverable'] & (df_all['is_indel'])]) / len(df_all[(df_all['is_indel'])]) * 100))
print("Discoverable frameshift: %d (%.2f %%)" % (len(df_all[df_all['is_discoverable'] & df_all['is_frameshift']]), len(df_all[df_all['is_discoverable'] & (df_all['is_frameshift'])]) / len(df_all[(df_all['is_frameshift'])]) * 100))

print("Discoverable in single-variant peptides only: %.2f %%" % (len(df_disc[df_disc['pep_type'].str.contains('single') & ~df_disc['pep_type'].str.contains('multi') & ~df_disc['pep_type'].str.contains('frameshift')]) / len(df_disc) * 100))
print("Discoverable in single-variant  and multi-variant peptides: %.2f %%" % (len(df_disc[df_disc['pep_type'].str.contains('single') & df_disc['pep_type'].str.contains('multi') & ~df_disc['pep_type'].str.contains('frameshift')]) / len(df_disc) * 100))
print("Discoverable in multi-variant peptides only: %.2f %%" % (len(df_disc[~df_disc['pep_type'].str.contains('single') & df_disc['pep_type'].str.contains('multi') & ~df_disc['pep_type'].str.contains('frameshift')]) / len(df_disc) * 100))
print("Discoverable in frameshift peptides only: %.2f %%" % (len(df_disc[~df_disc['pep_type'].str.contains('single') & ~df_disc['pep_type'].str.contains('multi') & df_disc['pep_type'].str.contains('frameshift')]) / len(df_disc) * 100))
print("Discoverable in frameshift peptides: %.2f %%" % (len(df_disc[df_disc['pep_type'].str.contains('frameshift')]) / len(df_disc) * 100))

print(df_disc[~df_disc['pep_type'].str.contains('single') & ~df_disc['pep_type'].str.contains('multi') & ~df_disc['pep_type'].str.contains('frameshift')][['variantID', 'pep_type']])
