# <i>Listeria monocytogenes</i> strain typing workflow

This pipeline has been tested on an Apple Silicon Mac using a base snakemake environment with `snakemake v7.32.4` and `python v3.11.15`.

The conda/mamba environment creation is what will change if you have a different machine (Mac Apple Silion, Mac Intel, Linux, etc.). This pipeline will not work on Windows without substantial changes.

Download the `mamba` package manager, then create the base environment with the following command:

```bash
mamba env create -f envs/snakemake.yaml
```

Then activate this environment, then run the pipeline:

```bash
mamba activate snakemake
bash run_smk.sh
```

## Workflow

The full workflow for a single sample is shown in `dag.png`. The steps are:

1. Download FASTQ files using `fasterq-dump` and trim adapters and merge overlapping paired reads using `fastp`.
2. Assemble reads using `spades`. It uses the merged reads and the remaining paired reads that did not overlap enough to be merged.
3. Align the 3 read files to the draft assembly using `bwa` and compute coverage statistics using `samtools`.
4. Perform taxonomic classification on the assembled contigs using `kraken2`.
5. Filter out low-quality contigs using the depth and taxonomic classification information from steps 3 and 4. The default parameters are a minimum length of 500 base pairs and a median coverage of 10x to keep a contig. After that, only contigs that are assigned by `kraken2` to any of the following are retained: the <i>L. monocytogenes</i> taxid (1639), all parent and children taxids of 1639, or the Listeria genus. Contigs assigned to the same genus are retained because `kraken2` sometimes misclassifies them as different species (e.g., <i>Listeria innocua</i> instead of <i>Listeria monocytogenes</i>.)
6. Compute assembly statistics using `quast` on both the pre- and post-filtered assemblies.
7. Run multi-locus sequence typing (MLST) using `mlst`.
8. Compare the MLST results with the list of 7-gene MLST profiles downloaded from L'Institut Pasteur to assign lineages, sequence types (ST), and clonal complexes (CC). The allele profiles were downloaded from <a href="https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/2/profiles_csv" target="_blank">here</a> in November 2026.

## Inputs

The workflow requires an input CSV file with two columns: `SampleName` is the sample name (it can be whatever you want) and `Run` is the SRA ID to fetch the FASTQ files using `fasterq-dump`. `Run` must be a valid SRA ID, otherwise it will fail to download the FASTQ files.