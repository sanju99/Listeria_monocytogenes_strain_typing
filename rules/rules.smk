import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"

primary_directory = "/Users/szk017/Documents/Projects/Listeria"

def get_sra_id(wildcards):
    return sample_SRA_dict[wildcards.sample_ID]



rule download_FASTQ:
    output:
        fastq1 = f"{sample_out_dir}/{{sample_ID}}_1.fastq.gz",
        fastq2 = f"{sample_out_dir}/{{sample_ID}}_2.fastq.gz",
    params:
        sample_out_dir = sample_out_dir,
        SRA_ID = get_sra_id,
    conda:
        f"{primary_directory}/envs/read_processing.yaml"
    threads:
        8
    shell:
        """
        fastq-dump --gzip --split-files --outdir {params.sample_out_dir} {params.SRA_ID}

        mv {params.sample_out_dir}/{params.SRA_ID}_1.fastq.gz {output.fastq1}
        mv {params.sample_out_dir}/{params.SRA_ID}_2.fastq.gz {output.fastq2}
        """



rule download_FASTQ_trim_adapters:
    output:
        fastq1 = temp(f"{sample_out_dir}/{{sample_ID}}_1.fastq"),
        fastq2 = temp(f"{sample_out_dir}/{{sample_ID}}_2.fastq"),
        fastq1_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R1.fastq.gz",
        fastq2_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R2.fastq.gz",
        merged_reads = f"{sample_out_dir}/fastp/{{sample_ID}}_mergedReads.fastq.gz",
        fastp_html = f"{sample_out_dir}/fastp/{{sample_ID}}.html",
        fastp_json = f"{sample_out_dir}/fastp/{{sample_ID}}.json"
    params:
        sample_out_dir = sample_out_dir,
        SRA_ID = get_sra_id,
        min_read_length = config["min_read_length"],
        poly_G_length = config["poly_G_length"],
        adapters_file = "references/NexteraPE-PE.fa",
    conda:
        f"{primary_directory}/envs/read_processing.yaml"
    threads:
        8
    shell:
        """
        fasterq-dump --split-files --outdir {params.sample_out_dir} {params.SRA_ID}

        # change the file names so that we can use them as declared in the input block
        mv {params.sample_out_dir}/{params.SRA_ID}_1.fastq {output.fastq1}
        mv {params.sample_out_dir}/{params.SRA_ID}_2.fastq {output.fastq2}

        fastp -i {output.fastq1} -I {output.fastq2} \
              -o {output.fastq1_trimmed} -O {output.fastq2_trimmed} \
              -h {output.fastp_html} -j {output.fastp_json} \
              --length_required {params.min_read_length} \
              --adapter_fasta {params.adapters_file} \
              --trim_poly_g --poly_g_min_len {params.poly_G_length} \
              --merge --merged_out {output.merged_reads} \
              --thread {threads}
        """




rule spades_assembly:
    input:
        fastq1_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R1.fastq.gz",
        fastq2_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R2.fastq.gz",
        merged_reads = f"{sample_out_dir}/fastp/{{sample_ID}}_mergedReads.fastq.gz",
    output:
        f"{sample_out_dir}/spades/contigs.fasta",
    params:
        spades_mode = config['spades_mode'],
        output_dir = f"{sample_out_dir}/spades",
    conda:
        f"{primary_directory}/envs/genome_assembly.yaml"
    threads:
        8
    shell:
        """
        spades.py \
            --{params.spades_mode} \
            -1 {input.fastq1_trimmed} \
            -2 {input.fastq2_trimmed} \
            --merged {input.merged_reads} \
            -o {params.output_dir} \
            --threads {threads}
        """


rule compute_draft_assembly_stats:
    input:
        contigs_file = f"{sample_out_dir}/spades/contigs.fasta",
    output:
        stats_file = f"{sample_out_dir}/spades/seqkit_stats.tsv",
    conda:
        f"{primary_directory}/envs/genome_assembly.yaml"
    threads:
        1
    shell:
        """
        seqkit stats --all {input.contigs_file} --basename --tabular > {output.stats_file}
        """



# finally, create CRAM file to cut down on disk space
# samtools view -T {input.contigs_file} -C -o {output.combined_cram} {output.combined_bam}
# samtools index -@ {threads} {output.combined_cram}


rule map_reads_to_ASM_compute_metrics:
    input:
        contigs_file = f"{sample_out_dir}/spades/contigs.fasta",
        fastq1_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R1.fastq.gz",
        fastq2_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R2.fastq.gz",
        merged_reads = f"{sample_out_dir}/fastp/{{sample_ID}}_mergedReads.fastq.gz",
    output:
        paired_reads_bam = temp(f"{sample_out_dir}/alignment/paired_reads.bam"),
        merged_reads_bam = temp(f"{sample_out_dir}/alignment/merged_reads.bam"),
        combined_bam = temp(f"{sample_out_dir}/alignment/bwa_mapped_sorted.bam"),
        combined_bam_index = temp(f"{sample_out_dir}/alignment/bwa_mapped_sorted.bam.bai"),

        depth_file = f"{sample_out_dir}/alignment/{{sample_ID}}.depth.tsv.gz",
        flagstat_file = f"{sample_out_dir}/alignment/{{sample_ID}}.flagstat.tsv",
        idxstats_file = f"{sample_out_dir}/alignment/{{sample_ID}}.idxstats.tsv",
    conda:
        f"{primary_directory}/envs/genome_assembly.yaml"
    threads:
        8
    shell:
        """
        bwa index {input.contigs_file}

        bwa mem -t {threads} {input.contigs_file} {input.fastq1_trimmed} {input.fastq2_trimmed} | samtools view -b -@ {threads} | samtools sort -@ {threads} -o {output.paired_reads_bam}

        bwa mem -t {threads} {input.contigs_file} {input.merged_reads} | samtools view -b -@ {threads} | samtools sort -@ {threads} -o {output.merged_reads_bam}

        # merge BAMs and index
        samtools merge {output.paired_reads_bam} {output.merged_reads_bam} -@ {threads} -o {output.combined_bam}
        samtools index -@ {threads} {output.combined_bam}

        # compute coverage at all sites
        samtools depth -aa {output.combined_bam} | gzip -c > {output.depth_file}
    
        # this is to get global assembly statistics. 
        # -O tsv makes it easier to read into pandas for analysis
        samtools flagstat -@ {threads} -O tsv {output.combined_bam} > {output.flagstat_file}
    
        # this is to get per-contig information.
        samtools idxstats {output.combined_bam} > {output.idxstats_file}
        """



rule kraken2:
    input:
        contigs_file = f"{sample_out_dir}/spades/contigs.fasta",
    output:
        kraken_report = f"{sample_out_dir}/kraken/report.txt",
        kraken_classifications = f"{sample_out_dir}/kraken/classifications.tsv",
    params:
        output_dir = f"{sample_out_dir}/spades",
        kraken_db = config['kraken_db'],
    conda:
        f"{primary_directory}/envs/genome_assembly.yaml"
    threads:
        8
    shell:
        """
        kraken2 --db {params.kraken_db} \
                --use-names \
                --report {output.kraken_report} \
                --output {output.kraken_classifications} \
                {input.contigs_file} \
                --threads {threads}
        """


rule filter_assembly_contigs:
    input:
        contigs_file = f"{sample_out_dir}/spades/contigs.fasta",
        depth_file = f"{sample_out_dir}/alignment/{{sample_ID}}.depth.tsv.gz",
        kraken_classifications = f"{sample_out_dir}/kraken/classifications.tsv",
    output:
        filtered_contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    params:
        taxid = config['taxid'],
        min_contig_length = 500,
        min_contig_median_cov = 10,
        kraken_db = config['kraken_db'],
    shell:
        """
        python scripts/filter_contigs.py \
                --fasta {input.contigs_file} \
                --depth {input.depth_file} \
                --kraken {input.kraken_classifications} \
                --taxid {params.taxid} \
                --output {output.filtered_contigs_file} \
                --min_contig_length {params.min_contig_length} \
                --min_median_coverage {params.min_contig_median_cov} \
                --kraken_db {params.kraken_db}
        """



# seqkit stats --all {input.contigs_file} --basename --tabular > {output.stats_file}

rule compute_filtered_assembly_stats:
    input:
        contigs_file = f"{sample_out_dir}/spades/contigs.fasta",
        filtered_contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        stats_file = f"{sample_out_dir}/assembly/quast/report.tsv",
    params:
        quast_dir = f"{sample_out_dir}/assembly/quast",
    conda:
        f"{primary_directory}/envs/genome_annotation.yaml"
    threads:
        1
    shell:
        """
        quast.py --labels contigs,filtered_contigs --threads {threads} -o {params.quast_dir} {input.contigs_file} {input.filtered_contigs_file}
        """



rule compute_downloaded_assembly_stats:
    input:
        filtered_contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        stats_file = f"{sample_out_dir}/assembly/quast/report.tsv",
    params:
        quast_dir = f"{sample_out_dir}/assembly/quast",
    conda:
        f"{primary_directory}/envs/genome_annotation.yaml"
    threads:
        1
    shell:
        """
        quast.py --labels assembly --threads {threads} -o {params.quast_dir} {input.filtered_contigs_file}
        """



rule run_mlst:
    input:
        contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        mlst_file = f"{sample_out_dir}/strain_typing/mlst.csv"
    threads:
        1
    conda:
        f"{primary_directory}/envs/strain_typing.yaml"
    shell:
        """
        # --full includes the header. -q is quiet mode
        mlst {input.contigs_file} --full --csv --outfile {output.mlst_file} -q
        """



rule lineage_CC_typing:
    input:
        mlst_file = f"{sample_out_dir}/strain_typing/mlst.csv"
    output:
        lineage_CC_file = f"{sample_out_dir}/strain_typing/lineage_CC.csv"
    params:
        mlst_profiles = config['mlst_profiles'],
    shell:
        """
        python scripts/get_lineage_CC_from_mlst.py -i {input.mlst_file} -I {params.mlst_profiles} -o {output.lineage_CC_file}
        """



rule make_blast_cgMLST1748_DB:
    output:
        blast_DB_files = [f"{config['cgMLST1748_blast_DB']}.{suffix}" for suffix in ['ndb', 'nhr', 'nin', 'njs', 'not', 'nsq', 'ntf', 'nto']],
    params:
        cgMLST1748_dir = directory(config['cgMLST1748_dir']),
        cgMLST1748_blast_DB = config['cgMLST1748_blast_DB']
    threads:
        1
    conda:
        f"{primary_directory}/envs/strain_typing.yaml"
    shell:
        """
        # cat all the FASTAs into 
        cat {params.cgMLST1748_dir}/lmo*.fasta > {params.cgMLST1748_dir}/full_DB.fasta

        makeblastdb -in {params.cgMLST1748_dir}/full_DB.fasta -dbtype nucl -out {params.cgMLST1748_blast_DB}

        rm {params.cgMLST1748_dir}/full_DB.fasta
        """


rule blast_cgMLST1748_DB:
    input:
        contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        cgMLST1748_file = f"{sample_out_dir}/strain_typing/cgMLST1748_blast_results.tsv.gz",
    params:
        cgMLST1748_blast_DB = config['cgMLST1748_blast_DB'],
    threads:
        8
    conda:
        f"{primary_directory}/envs/strain_typing.yaml"
    shell:
        """
        blastn -query {input.contigs_file} \
               -db {params.cgMLST1748_blast_DB} \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
               -num_threads {threads} \
               -max_target_seqs 1000000 \
               -perc_identity 70 | gzip -c > {output.cgMLST1748_file}
        """


rule run_snippy:
    input:
        contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        snippy_dir = directory(f"{sample_out_dir}/{{sample_ID}}"),
    params:
        ref_genome = config['ref_genome'],
    threads:
        8
    shell:
        """
        source ~/miniforge3/etc/profile.d/conda.sh
        set +u
        conda activate snippy
        set -u

        snippy --outdir {output.snippy_dir} --ctgs {input.contigs_file} --ref {params.ref_genome} --basequal 20 --mapqual 30 --mincov 10 --minfrac 0.95 --cleanup --cpus {threads}
        """



rule snippy_core_genome_aln:
    input:
        snippy_dirs=expand(f"{sample_out_dir}/{{sample_ID}}", sample_ID=df_samples_runs['SampleName'].values),
    output:
        core_genome_aln = f"phylogeny/core_genome.aln"
    params:
        ref_genome = config['ref_genome'],
    shell:
        """
        source ~/miniforge3/etc/profile.d/conda.sh
        set +u
        conda activate snippy
        set -u

        snippy-core --prefix phylogeny/core_genome --ref {params.ref_genome} {input.snippy_dirs}
        """


rule download_bakta_db:
    params:
        bakta_db = config['bakta_db'],
        bakta_db_type = config['bakta_db_type'],
    output:
        bakta_db = directory(f"{config['bakta_db']}/db-{config['bakta_db_type']}"),
    conda:
        f"{primary_directory}/envs/genome_annotation.yaml"
    shell:
        """
        bakta_db download --output {params.bakta_db} --type {params.bakta_db_type}
        """


rule run_bakta:
    input:
        contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        genome_annotation_dir = directory(f"{sample_out_dir}/genome_annotation/bakta"),
        gff_file = f"{sample_out_dir}/genome_annotation/bakta/{{sample_ID}}.gff3",
    params:
        sample_ID = f"{{sample_ID}}",
        bakta_db = f"{config['bakta_db']}/db-{config['bakta_db_type']}",
        genus = config['genus'],
        species = config['species'],
    threads:
        8
    conda:
        f"{primary_directory}/envs/genome_annotation.yaml"
    shell:
        """
        bakta --db {params.bakta_db} \
              --prefix {params.sample_ID} \
              --output {output.genome_annotation_dir} \
              --genus {params.genus} --species {params.species} \
              --threads {threads} --force \
              {input.contigs_file}
        """