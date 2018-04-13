# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using MetaPhlAn2

rule download_metaphlan2_database:
    """Download MetaPhlAn2 db_v20_m200"""
    output:
        config["dbdir"]+"/metaphlan2/mpa_v20_m200.fna",
        config["dbdir"]+"/metaphlan2/mpa_v20_m200.pkl",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan2.yaml"
    params:
        dbdir=config["dbdir"]+"/metaphlan2"
    shell:
        """
        cd {params.dbdir}
        wget https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.tar \
        && \
        tar -xf mpa_v20_m200.tar \
        && \
        bunzip2 mpa_v20_m200.fna.bz2 \
        && \
        rm -v mpa_v20_m200.tar
        """


rule build_metaphlan2_index:
    """Build MetaPhlAn2 bowtie2 index."""
    input:
        config["dbdir"]+"/metaphlan2/mpa_v20_m200.fna"
    output:
        [config["dbdir"]+"/metaphlan2/mpa_v20_m200.{n}.bt2".format(n=num) for num in (1,2,3,4)],
        [config["dbdir"]+"/metaphlan2/mpa_v20_m200.rev.{n}.bt2".format(n=num) for num in (1,2)],
    log:
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan2.yaml"
    threads:
        4
    params:
        dbdir=config["dbdir"]+"/metaphlan2"
    shell:
        """
        cd {params.dbdir}
        bowtie2-build \
            mpa_v20_m200.fna \
            mpa_v20_m200 \
            --threads {threads}
        """


mpa_config = config["metaphlan2"]
rule metaphlan2:
    """Taxonomic profiling using MetaPhlAn2."""
    input:
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        bt2_out=config["outdir"]+"/metaphlan2/{sample}.bowtie2.bz2",
        mpa_out=config["outdir"]+"/metaphlan2/{sample}.metaphlan2.txt",
    log:
        stdout=config["outdir"]+"/logs/metaphlan2/{sample}.metaphlan2.stdout.log",
        stderr=config["outdir"]+"/logs/metaphlan2/{sample}.metaphlan2.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan2.yaml"
    threads:
        4
    params:
        mpa_pkl=mpa_config["mpa_pkl"],
        bt2_db_prefix=mpa_config["bt2_db_prefix"],
        extra=mpa_config["extra"],
    shell:
        """
        metaphlan2.py \
            --input_type fastq \
            --nproc {threads} \
            --sample_id {wildcards.sample} \
            --bowtie2out {output.bt2_out} \
            --mpa_pkl {params.mpa_pkl} \
            --bowtie2db {params.bt2_db_prefix} \
            {input.read1},{input.read2} \
            {output.mpa_out} \
            {params.extra} \
            > {log.stdout} \
            2> {log.stderr}
        """


rule combine_metaphlan2_outputs:
    """Combine metaphlan2 outputs into a large table."""
    input:
        expand(config["outdir"]+"/metaphlan2/{sample}.metaphlan2.txt", sample=SAMPLES)
    output:
        config["outdir"]+"/metaphlan2/all_samples.metaphlan2.txt"
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan2.yaml"
    threads:
        1
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """

