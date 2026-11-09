import pandas as pd
import sys, subprocess, os, glob, argparse
from io import StringIO

parser = argparse.ArgumentParser()

# use edirect_env environment
parser.add_argument("-i", "--input", dest='input_file', type=str, required=True, help='Input CSV file of IDs to search. Must be a single column with no header')
parser.add_argument("-o", "--output", dest='output_file', type=str, required=True, help='Output CSV file for the returned metadata')
parser.add_argument("-g", "--genus", dest='genus_to_search', type=str, required=True, help='Genus of the organism to search for')

cmd_line_args = parser.parse_args()

input_file = cmd_line_args.input_file
output_file = cmd_line_args.output_file
genus_to_search = cmd_line_args.genus_to_search

if input_file[-4:] in ['.txt', '.tsv']:
    search_lst = pd.read_csv(input_file, header=None, sep='\t')[0].values
else:
    search_lst = pd.read_csv(input_file, header=None)[0].values

if os.path.isfile(output_file):
    results_df = pd.read_csv(output_file)
    print(f"Finished querying {results_df['Query'].nunique()}/{len(search_lst)} samples")

    if results_df['Query'].nunique() == len(search_lst):
        exit()
else:
    results_df = pd.DataFrame(columns=['Query'])

print(f"Total samples to query: {len(search_lst)}")
all_data = []

for i, query_name in enumerate(search_lst):
    command = f'esearch -db sra -query "{genus_to_search} AND {query_name}" | efetch -format runinfo'

    if query_name not in results_df['Query'].values:
        try:
            res = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=30
            )
        
            output = res.stdout
        
        except subprocess.TimeoutExpired:
            print(f"Did not find {query_name} in SRA")
            output = ""
        
        if output.strip() != "":
            res_df = pd.read_csv(StringIO(output))
            res_df["Query"] = query_name
            all_data.append(res_df)
        else:
            all_data.append(pd.DataFrame({"Query": [query_name]}))
    
    # if query_name not in results_df['Query'].values:
        
    #     res = subprocess.Popen(command, shell=True, encoding='utf8', stdout=subprocess.PIPE)
    #     res = res.communicate()[0]

        # if output != '':
        #     res_df = pd.read_csv(StringIO(res))
        #     res_df['Query'] = [query_name]*len(res_df)
        #     all_data.append(res_df)
        # else:
        #     all_data.append(pd.DataFrame({"Query": query_name}, index=[i]))
        
    if i % 10 == 0:
        print(f"Finished {i}")

        if len(all_data) > 0:
            pd.concat([results_df, pd.concat(all_data, axis=0)], axis=0).to_csv(output_file, index=False)
        else:
            results_df.to_csv(output_file, index=False)
    
pd.concat([results_df, pd.concat(all_data, axis=0)], axis=0).to_csv(output_file, index=False)