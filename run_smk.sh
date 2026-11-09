#!/bin/bash 

# source activate snakemake

# snakemake --use-conda --conda-create-envs-only --conda-frontend mamba --configfile config.yaml --cores 18

snakemake --snakefile snakefile \
          --configfile config.yaml \
          --use-conda --conda-frontend mamba \
          --cores 16 --resources mem_mb=36000 \
          --rerun-incomplete --keep-going \
          --unlock


snakemake --snakefile snakefile \
          --configfile config.yaml \
          --use-conda --conda-frontend mamba \
          --cores 16 --resources mem_mb=36000 \
          --rerun-incomplete --keep-going #--allowed-rules run_mlst lineage_CC_typing compute_downloaded_assembly_stats