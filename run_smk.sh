#!/bin/bash 

# snakemake --use-conda --conda-create-envs-only --conda-frontend mamba --configfile config.yaml --cores 18

# unlock the directory in case there are any locks
snakemake --snakefile snakefile \
          --configfile config.yaml \
          --use-conda --conda-frontend mamba \
          --cores 16 --resources mem_mb=36000 \
          --rerun-incomplete --keep-going \
          --unlock

# run the pipeline
snakemake --snakefile snakefile \
          --configfile config.yaml \
          --use-conda --conda-frontend mamba \
          --cores 16 --resources mem_mb=36000 \
          --rerun-incomplete --keep-going --rerun-triggers mtime