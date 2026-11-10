rule compute_downloaded_assembly_stats:
    input:
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
        quast.py --labels assembly --threads {threads} -o {params.quast_dir} {input.filtered_contigs_file}
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
        "../envs/strain_typing.yaml"
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
        "../envs/strain_typing.yaml"
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
        "../envs/genome_annotation.yaml"
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
        "../envs/genome_annotation.yaml"
    shell:
        """
        bakta --db {params.bakta_db} \
              --prefix {params.sample_ID} \
              --output {output.genome_annotation_dir} \
              --genus {params.genus} --species {params.species} \
              --threads {threads} --force \
              {input.contigs_file}
        """