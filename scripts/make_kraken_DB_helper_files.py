import argparse
import pandas as pd


def load_nodes_dmp(kraken_DB):
    parent_map = {}
    rank_map = {}
    name_map = {}

    with open(f"{kraken_DB}/nodes.dmp", 'r') as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            taxid = int(parts[0])
            parent_map[taxid] = int(parts[1])
            rank_map[taxid] = parts[2]

    # names.dmp
    with open(f"{kraken_DB}/names.dmp", 'r') as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            taxid = int(parts[0])
            name = parts[1]
            name_class = parts[3]

            if name_class == "scientific name":
                name_map[taxid] = name
                
    return parent_map, rank_map, name_map



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--kraken_DB", dest='kraken_DB', type=str, required=True, help='Full path to a kraken database. nodes.dmp MUST be present in this directory')

    args = parser.parse_args()
    kraken_DB = args.kraken_DB

    # get dictionaries mapping each taxid to its parent and each taxid to its rank (species, genus, family, etc.)
    parent_map, rank_map, name_map = load_nodes_dmp(kraken_DB)
    
    # save to CSV files in the same directory
    pd.DataFrame(parent_map.items(), columns=['TaxID', 'Parent_TaxID']).to_csv(f"{kraken_DB}/parent_taxa_mapping.csv.gz", index=False, compression='gzip')
    pd.DataFrame(rank_map.items(), columns=['TaxID', 'Rank']).to_csv(f"{kraken_DB}/taxid_rank_mapping.csv.gz", index=False, compression='gzip')
    pd.DataFrame(name_map.items(), columns=['TaxID', 'TaxName']).to_csv(f"{kraken_DB}/taxid_name_mapping.csv.gz", index=False, compression='gzip')


if __name__ == "__main__":
    main()