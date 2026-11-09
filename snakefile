import numpy as np
import pandas as pd

configfile: "config.yaml"

# Define the list of samples, either directly or through config
output_dir = config["output_dir"]

# dataframe of isolates to run. First column is the sample ID, second column is a comma-separated string of the sequencing run IDs
df_samples_runs = pd.read_csv(config['samples_file'])
sample_SRA_dict = dict(zip(df_samples_runs['SampleName'], df_samples_runs['Run']))

include: "rules/main_rules.smk"

# Define rules in the Snakefile

rule all:
    input:
        [f"{output_dir}/{sample_ID}/strain_typing/lineage_CC.csv" for sample_ID in sample_SRA_dict.keys()],
        [f"{output_dir}/{sample_ID}/assembly/quast/report.tsv" for sample_ID in sample_SRA_dict.keys()],
        # [f"{output_dir}/{sample_ID}/genome_annotation/bakta/{sample_ID}.gff3" for sample_ID in sample_SRA_dict.keys()],
        # [f"{output_dir}/{sample_ID}/strain_typing/cgMLST1748_blast_results.tsv.gz" for sample_ID in sample_SRA_dict.keys()],