import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"

# function to get the SRA ID for a given sample from the dictionary
def get_sra_id(wildcards):
    return sample_SRA_dict[wildcards.sample_ID]



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
        "../envs/read_processing.yaml"
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
        "../envs/genome_assembly.yaml"
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




rule map_reads_to_ASM_compute_metrics:
    input:
        contigs_file = f"{sample_out_dir}/spades/contigs.fasta",
        fastq1_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R1.fastq.gz",
        fastq2_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}_R2.fastq.gz",
        merged_reads = f"{sample_out_dir}/fastp/{{sample_ID}}_mergedReads.fastq.gz",
    output:
        paired_reads_bam = temp(f"{sample_out_dir}/alignment/paired_reads.bam"),
        merged_reads_bam = temp(f"{sample_out_dir}/alignment/merged_reads.bam"),
        
        combined_bam = temp(f"{sample_out_dir}/alignment/{{sample_ID}}.bam"),
        combined_bam_index = temp(f"{sample_out_dir}/alignment/{{sample_ID}}.bam.bai"),
        
        combined_cram = f"{sample_out_dir}/alignment/{{sample_ID}}.cram",
        combined_cram_index = f"{sample_out_dir}/alignment/{{sample_ID}}.cram.crai",

        depth_file = f"{sample_out_dir}/alignment/{{sample_ID}}.depth.tsv.gz",
        flagstat_file = f"{sample_out_dir}/alignment/{{sample_ID}}.flagstat.tsv",
        idxstats_file = f"{sample_out_dir}/alignment/{{sample_ID}}.idxstats.tsv",
    conda:
        "../envs/genome_assembly.yaml"
    threads:
        8
    shell:
        """
        # index reference genome (which is the spades assembly)
        bwa index {input.contigs_file}

        # align the paired reads to the reference
        bwa mem -t {threads} {input.contigs_file} {input.fastq1_trimmed} {input.fastq2_trimmed} | samtools view -b -@ {threads} | samtools sort -@ {threads} -o {output.paired_reads_bam}

        # align the merged reads to the reference
        bwa mem -t {threads} {input.contigs_file} {input.merged_reads} | samtools view -b -@ {threads} | samtools sort -@ {threads} -o {output.merged_reads_bam}

        # merge the two BAMs and index
        samtools merge {output.paired_reads_bam} {output.merged_reads_bam} -@ {threads} -o {output.combined_bam}
        samtools index -@ {threads} {output.combined_bam}

        # compute coverage at all sites
        samtools depth -aa {output.combined_bam} | gzip -c > {output.depth_file}
    
        # this is to get global assembly statistics. 
        # -O tsv makes it easier to read into pandas for analysis
        samtools flagstat -@ {threads} -O tsv {output.combined_bam} > {output.flagstat_file}
    
        # this is to get per-contig information.
        samtools idxstats {output.combined_bam} > {output.idxstats_file}

        # finally, create CRAM file to cut down on disk space
        samtools view -T {input.contigs_file} -C -o {output.combined_cram} {output.combined_bam}
        samtools index -@ {threads} {output.combined_cram}
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
        "../envs/genome_assembly.yaml"
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
        min_contig_length = config['min_contig_length'],
        min_contig_median_cov = config['min_contig_median_cov'],
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


rule compute_assembly_stats:
    input:
        contigs_file = f"{sample_out_dir}/spades/contigs.fasta",
        filtered_contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        stats_file = f"{sample_out_dir}/assembly/quast/report.tsv",
    params:
        quast_dir = f"{sample_out_dir}/assembly/quast",
    conda:
        "../envs/genome_annotation.yaml"
    threads:
        1
    shell:
        """
        quast.py --labels contigs,filtered_contigs --threads {threads} -o {params.quast_dir} {input.contigs_file} {input.filtered_contigs_file}
        """



rule run_mlst:
    input:
        contigs_file = f"{sample_out_dir}/assembly/contigs.fasta",
    output:
        mlst_file = f"{sample_out_dir}/strain_typing/mlst.csv"
    threads:
        1
    conda:
        "../envs/strain_typing.yaml"
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