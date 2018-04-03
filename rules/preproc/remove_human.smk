# vim: syntax=python expandtab
# Rules to filter human sequences from metagenomic reads

rule download_hg19:
    """Downloads masked hg19 from: 
    https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk"""
    output:
        config["dbdir"]+"/hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    conda:
        "../../envs/bbmap.yaml"
    shell:
        """
        scripts/download_from_gdrive.py \
            -o {output}.gz \
            0B3llHR93L14wd0pSSnFULUlhcUk \
        && \
        gunzip {output}.gz \
        && \
        bbmap.sh \
            ref={output} \
            path={config[dbdir]}/hg19/
        """


rule remove_human:
    input:
        read1=config["outdir"]+"/trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
        read2=config["outdir"]+"/trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
        hg19=config["dbdir"]+"/hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
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
            path={config[dbdir]}/hg19 \
            outu1={output.read1} \
            outu2={output.read2} \
            outm={output.human} \
            statsfile={log.statsfile} \
            minid={config[bbmap_removehuman_minid]} \
            maxindel={config[bbmap_removehuman_maxindel]} \
            minhits={config[bbmap_removehuman_minhits]} \
            bandwidthratio={config[bbmap_removehuman_bandwidthratio]} \
            bandwidth={config[bbmap_removehuman_bandwidth]} \
            qtrim={config[bbmap_removehuman_qtrim]} \
            trimq={config[bbmap_removehuman_trimq]} \
            {config[bbmap_removehuman_quickmatch]} \
            {config[bbmap_removehuman_fast]} \
            {config[bbmap_removehuman_untrim]} \
        """

