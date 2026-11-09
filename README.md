# README

This pipeline has been tested on an Apple Silicon Mac using a base snakemake environment with `snakemake v7.32.4` and `python v3.11.15`.

The conda/mamba environment creation is what will change if you have a different machine (Mac Apple Silion, Mac Intel, Linux, etc.). This pipeline will not work on Windows without substantial changes.

Create the base environment with the following command:

```bash
mamba env create -f envs/snakemake.yaml -n <INSERT YOUR DESIRED ENV NAME HERE>
```

Then activate this environment, then run the pipeline with:

```bash
bash run_smk.sh
```