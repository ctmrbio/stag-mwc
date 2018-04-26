# vim: syntax=python expandtab
# Rules to filter human sequences from metagenomic reads
from snakemake.exceptions import WorkflowError
import os.path

localrules:
    download_hg19

if not os.path.isdir(os.path.join(config["remove_human"]["hg19_path"], "ref")):
    err_message = "Cannot find hg19 database for human sequence removal at: '{}'!\n".format(config["remove_human"]["hg19_path"])
    err_message += "Specify path to folder containing BBMap index of hg19 files in config.yaml.\n"
    err_message += "Run 'snakemake index_hg19' to download and create a BBMap index in '{dbdir}/hg19'".format(dbdir=config["dbdir"])
    raise WorkflowError(err_message)

# Add final output files from this module to 'all_outputs' from the main
# Snakefile scope. SAMPLES is also from the main Snakefile scope.
filtered_human = expand("{outdir}/filtered_human/{sample}_R{readpair}.filtered_human.fq.gz",
        outdir=config["outdir"],
        sample=SAMPLES,
        readpair=[1,2])
all_outputs.extend(filtered_human)


rule download_hg19:
    """Download masked hg19 from: 
    https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk"""
    output:
        config["dbdir"]+"/hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    conda:
        "../../envs/stag-mwc.yaml"
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
        "../../envs/stag-mwc.yaml"
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
    """Filter reads matching hg19. NB: requires about 16GB of memory."""
    input:
        read1=config["outdir"]+"/trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
        read2=config["outdir"]+"/trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
    output:
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
        human=config["outdir"]+"/filtered_human/{sample}_human.fq.gz",
    log:
        statsfile=config["outdir"]+"/logs/remove_human/{sample}.statsfile.txt",
        stderr=config["outdir"]+"/logs/remove_human/{sample}.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        4
    params:
        minid=rh_config["minid"],
        maxindel=rh_config["maxindel"],
        minhits=rh_config["minhits"],
        bandwidthratio=rh_config["bandwidthratio"],
        bandwidth=rh_config["bandwidth"],
        qtrim=rh_config["qtrim"],
        trimq=rh_config["trimq"],
        quickmatch=rh_config["quickmatch"],
        fast=rh_config["fast"],
        untrim=rh_config["untrim"],
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
            minid={params.minid} \
            maxindel={params.maxindel} \
            minhits={params.minhits} \
            bandwidthratio={params.bandwidthratio} \
            bandwidth={params.bandwidth} \
            qtrim={params.qtrim} \
            trimq={params.trimq} \
            {params.quickmatch} \
            {params.fast} \
            {params.untrim} \
            2> {log.stderr}
        """

