import numpy as np
import pandas as pd
import os, glob, re, math, shutil
from Bio import Seq, SeqIO

########################## STEP 1: READ IN TABLE FROM READ TO SEARCH IDS IN SRA ##########################

# data table from Reed
df = pd.read_excel("./data/listeria_wwtp_mlst.xlsx", sheet_name=None)

sample_IDs = np.concatenate([df[sheet_name]['EHSID'] for sheet_name in df.keys()])
sample_IDs = np.unique(sample_IDs)
print(f"{len(sample_IDs)} unique IDs")

# search using search_SRA.py to get the SRA IDs and FASTQs
pd.Series(sample_IDs).to_csv("./data/search_SRA.csv", header=None, index=False)

########################## STEP 2: GET THE SRA METADATA AND PERFORM CHECKS ##########################

df_SRA = pd.read_csv("./data/search_SRA_metadata.csv")
print(len(df_SRA), len(df_SRA.dropna(subset='BioSample')))

assert len(df_SRA.query("LibrarySource != 'GENOMIC'")) == 0
assert len(df_SRA.query("LibraryLayout != 'PAIRED'")) == 0
assert len(df_SRA.query("Platform != 'ILLUMINA'")) == 0

print(df_SRA.CenterName.value_counts(dropna=False))

# check that all internal IDs were found in SRA
assert len(set(sample_IDs) - set(df_SRA.dropna(subset='BioSample').Query)) == 0

# this is for passing into the snakemake pipeline
df_SRA[['SampleName', 'Run']].to_csv("./samples.csv", index=False)

########################## STEP 3: COLLECTION METADATA FOR THE WASTEWATER SAMPLES ##########################

df = pd.read_excel("../data/from_Reed/listeria_wwtp_mlst.xlsx", sheet_name=None)
df_1 = df['all_listeria'].rename(columns={'MLST': 'MLST1', 'EHSID': 'acc_id'})
df_2 = df['dup_removed'].rename(columns={'Lineage ': 'Lineage', 'EHSID': 'acc_id'})

df_1[['CC', 'Lineage']] = df_1['CC, Lineage '].str.split(',', expand=True)
df_2[['CC', 'Lineage']] = df_2['CC, Lineage '].str.split(',', expand=True)
combine_cols = ['acc_id', 'ST', 'CC', 'Lineage', 'MLST1', 'MLST2', 'MLST3', 'MLST4', 'MLST5', 'MLST6', 'MLST7']

df_Reed = pd.concat([df_1, df_2])[combine_cols + ['Collection date', 'WWTP']].drop_duplicates().sort_values('acc_id')

df_Reed.loc[df_Reed['WWTP']=='M/L ', 'WWTP'] = 'Meadowbrook/Limestone'
df_Reed.loc[df_Reed['WWTP']=='Binghamton', 'WWTP'] = 'Binghamton - North'
df_Reed.loc[df_Reed['WWTP']=='OCSD - Harriman', 'WWTP'] = 'Harriman'
df_Reed.loc[df_Reed['WWTP']=='Town of Fallsburg', 'WWTP'] = 'South Fallsburg'
df_Reed.loc[df_Reed['WWTP']=='Whitehall WWTP', 'WWTP'] = 'Whitehall'
df_Reed.loc[df_Reed['WWTP']=='Rockland', 'WWTP'] = 'Rockland Co SD'
df_Reed.loc[df_Reed['WWTP']=='Rensselaer Co District 1 ', 'WWTP'] = 'Rensselaer'

wwtp_county_mapping = df_2[['WWTP', 'County']].drop_duplicates()
df_merged_data = pd.read_excel("data/from_Reed/merged_data_listeria_wwtp.xlsx", sheet_name='in')

df_location_metadata = df_Reed[['WWTP']].merge(wwtp_county_mapping).merge(df_merged_data[['sample_location_specify', 'wwtp_latitude', 'wwtp_longitude']].rename(columns={'sample_location_specify': 'WWTP'}))
df_location_metadata = df_location_metadata.sort_values(['WWTP', 'County', 'wwtp_latitude', 'wwtp_longitude'], na_position='last').drop_duplicates(subset=['WWTP'], keep='first')

# Read shapefile
ny_sewersheds = gpd.read_file("../data/NYS_sewersheds/New York State Sewersheds.shp")

# Reproject to WGS84 longitude/latitude
ny_sewersheds = ny_sewersheds.to_crs("EPSG:4326")

# replacements because of slight differences in namking
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='Binghamton (North)', 'Sewershed'] = 'Binghamton - North'
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='OCSD No 1 - Harriman WWTP', 'Sewershed'] = 'Harriman' # there are 2 with the name Harriman in it, but only one is in Orange County
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='Herkimer County Sewer District', 'Sewershed'] = 'Herkimer County SD'
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='Northwest quadrant', 'Sewershed'] = 'Northwest Quadrant'
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='Rockland Co Sewer District', 'Sewershed'] = 'Rockland Co SD'
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='Saratoga Co. Sewer District No. 1', 'Sewershed'] = 'Saratoga Co SD No 1'
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='Southtowns AWTF', 'Sewershed'] = 'Southtowns'
ny_sewersheds.loc[ny_sewersheds['Sewershed']=='Town of Victor and Town of Farmington', 'Sewershed'] = 'Town of Victor and Town Of Farmington'

df_location_metadata = df_location_metadata.merge(ny_sewersheds[['Sewershed', 'geometry']].rename(columns={'Sewershed': 'WWTP'}), on='WWTP')
df_location_metadata.to_csv("./data/WWTP_location_metadata.csv", index=False)

df_collection_metadata = df_Reed[['acc_id', 'Collection date', 'WWTP']].drop_duplicates().merge(df_location_metadata)
df_collection_metadata.to_csv("../data/collection_metadata.csv", index=False)

########################## STEP 4: BLAST RESULTS AGAINST THE CGMLST1748 SCHEME ##########################

df_blast_results = []

for sample in os.listdir("results"):

    df_blast = pd.read_csv(f"../results/{sample}/strain_typing/cgMLST1748_blast_results.tsv.gz", sep='\t', 
                           names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
                           )

    df_blast['Sample'] = sample
    df_blast_results.append(df_blast)

df_blast_results = pd.concat(df_blast_results)
print(len(df_blast_results))

# add information about the locus lengths so that we can keep only hits that cover at least 90% of the average locus length
# df_blast_results['s_aln_len'] = np.abs(df_blast_results['send'] - df_blast_results['sstart']) + 1
df_blast_results['locus'] = df_blast_results['sseqid'].str.split('_').str[0]
df_blast_results['allele'] = df_blast_results['sseqid'].str.split('_').str[1].astype(int)

# compute proportion of the locus' min length that is covered by the hit. min length is the smalles length of any of the alleles for that locus
# df_blast_results = df_blast_results.merge(df_cgmlst1748[['locus', 'min length']], how='left')
# df_blast_results['locus_prop_length'] = df_blast_results['s_aln_len'] / df_blast_results['min length']

# df_blast_results = df_blast_results.query("locus_prop_length >= 0.9")
# print(len(df_blast_results))

# df_blast_results.set_index(['Sample']).to_csv("./data/cgmlst1758_blast_results.csv.gz")

df_allele_typing = df_blast_results.sort_values(['Sample', 'locus', 'bitscore'], ascending=[True, True, False]).drop_duplicates(['Sample', 'locus'], keep='first')
df_allele_typing['allele'] = df_allele_typing['allele'].astype(str)
df_allele_typing.set_index(['Sample']).to_csv("../data/cgmlst1758_allele_typing.csv")