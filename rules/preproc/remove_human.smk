# vim: syntax=python expandtab
# Rules to filter human sequences from metagenomic reads

rule download_hg19:
    """Download masked hg19 from: 
    https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk"""
    output:
        config["dbdir"]+"/hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    conda:
        "../../envs/bbmap.yaml"
    params:
        dbdir=config["dbdir"]+"/hg19/"
    shell:
        """
        scripts/download_from_gdrive.py \
            -o {output}.gz \
            0B3llHR93L14wd0pSSnFULUlhcUk \
        && \
        gunzip {output}.gz
    """


rule index_hg19:
    """Create BBMap index of hg19 fasta file."""
    input:
        config["dbdir"]+"/hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    output:
        config["dbdir"]+"/hg19/ref/genome/1/chr1.chrom.gz",
        config["dbdir"]+"/hg19/ref/genome/1/chr2.chrom.gz",
        config["dbdir"]+"/hg19/ref/genome/1/chr3.chrom.gz",
        config["dbdir"]+"/hg19/ref/genome/1/chr4.chrom.gz",
        config["dbdir"]+"/hg19/ref/genome/1/chr5.chrom.gz",
        config["dbdir"]+"/hg19/ref/genome/1/chr6.chrom.gz",
        config["dbdir"]+"/hg19/ref/genome/1/chr7.chrom.gz",
        config["dbdir"]+"/hg19/ref/genome/1/info.txt",
        config["dbdir"]+"/hg19/ref/genome/1/scaffolds.txt.gz",
        config["dbdir"]+"/hg19/ref/genome/1/summary.txt",
        config["dbdir"]+"/hg19/ref/index/1/chr1-3_index_k13_c2_b1.block",
        config["dbdir"]+"/hg19/ref/index/1/chr1-3_index_k13_c2_b1.block2.gz",
        config["dbdir"]+"/hg19/ref/index/1/chr4-7_index_k13_c2_b1.block",
        config["dbdir"]+"/hg19/ref/index/1/chr4-7_index_k13_c2_b1.block2.gz",
    conda:
        "../../envs/bbmap.yaml"
    params:
        dbdir=config["dbdir"]+"/hg19/"
    shell:
        """
        bbmap.sh \
            ref={input} \
            path={params.dbdir}
        """


rh_config = config["remove_human"]
rule remove_human:
    input:
        read1=config["outdir"]+"/trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
        read2=config["outdir"]+"/trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
    output:
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
        human=config["outdir"]+"/filtered_human/{sample}_human.fq.gz",
    log:
        statsfile=config["outdir"]+"/logs/remove_human/{sample}.statsfile.txt",
    shadow:
        "shallow"
    conda:
        "../../envs/bbmap.yaml"
    threads:
        2
    shell:
        """
        bbmap.sh \
            threads={threads} \
            in1={input.read1} \
            in2={input.read2} \
            path={rh_config[hg19_path]} \
            outu1={output.read1} \
            outu2={output.read2} \
            outm={output.human} \
            statsfile={log.statsfile} \
            minid={rh_config[minid]} \
            maxindel={rh_config[maxindel]} \
            minhits={rh_config[minhits]} \
            bandwidthratio={rh_config[bandwidthratio]} \
            bandwidth={rh_config[bandwidth]} \
            qtrim={rh_config[qtrim]} \
            trimq={rh_config[trimq]} \
            {rh_config[quickmatch]} \
            {rh_config[fast]} \
            {rh_config[untrim]} \
        """

