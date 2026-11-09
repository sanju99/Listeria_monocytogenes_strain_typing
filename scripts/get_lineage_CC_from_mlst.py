#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", dest="mlst_results", type=str, required=True, help='Output of mlst on this sample')
parser.add_argument("-I", dest='mlst_reference_profiles', type=str, required=True, help="MLST profiles file from L'Institut Pasteur")
parser.add_argument("-o", dest="lineage_CC_file", type=str, required=True, help='Output CSV file to write the lineage and CC type to')


def infer_lineage_and_CC(df_mlst_profiles, df_mlst_allele_nums):
    '''
    Use this function for novel MLST profiles
    '''
    ref_alleles = df_mlst_profiles[df_mlst_allele_nums.columns].to_numpy(dtype=int)
    query_alleles = df_mlst_allele_nums.to_numpy(dtype=int)
    
    # Calculate allelic distances: rows = query isolates, columns = reference profiles
    distances = (query_alleles[:, None, :] != ref_alleles[None, :, :]).sum(axis=2)
    
    # Find the closest reference profile(s) for each isolate
    df_closest_profiles = []
    
    for i, isolate in enumerate(df_mlst_allele_nums.index):
        min_distance = distances[i].min()
        closest = df_mlst_profiles.loc[distances[i] == min_distance]
        closest['Distance'] = min_distance
        df_closest_profiles.append(closest)
    
    df_closest_profiles = pd.concat(df_closest_profiles)

    # there really shouldn't be multiple because they would be differing by so many alleles
    assert df_closest_profiles.Lineage.nunique() == 1

    return df_closest_profiles



def main():

    args = parser.parse_args()

    mlst_results = args.mlst_results
    mlst_reference_profiles = args.mlst_reference_profiles
    lineage_CC_file = args.lineage_CC_file

    # 7-gene MLST results for this sample
    df_mlst = pd.read_csv(mlst_results)

    # reference 7-gene MLST profiles mapped to different lineages and clonal complexes
    df_mlst_profiles = pd.read_csv(mlst_reference_profiles, sep='\t')

    # extract the allele numbers. Allow for tildes in front if it's an approximate match
    # this will handle both cases where the values are 6; 6,6; and 6,5
    # 6 returns 6. 6,6 returns 6, and 6,5 returns 6,5
    df_mlst_allele_nums = df_mlst['ALLELES'].str.findall(r'\(~?(\d+)(?:,(\d+))?\)').apply(lambda x: pd.Series([','.join(dict.fromkeys(filter(None, pair))) for pair in x])).astype(str)
    df_mlst_allele_nums.columns = df_mlst['ALLELES'].str.findall(r'([a-zA-Z]+)\(').apply(pd.Series).iloc[0, :].values

    # because of the possibility of values like 6,5, the numbers in df_mlst_allele_nums are strings, but the numbers in df_mlst_profiles are integers. So need to convert
    df_mlst_profiles[df_mlst_allele_nums.columns] = df_mlst_profiles[df_mlst_allele_nums.columns].astype(str)

    # merge with the profiles file.
    df_merged = df_mlst_allele_nums.merge(df_mlst_profiles, on=list(df_mlst_allele_nums.columns))

    # when the MLST profile is novel, then get the reference profiles that are the closest distance away
    if len(df_merged) == 0:
        df_merged = infer_lineage_and_CC(df_mlst_profiles, df_mlst_allele_nums)
    
    df_merged.to_csv(lineage_CC_file, index=False)

if __name__ == "__main__":
    main()