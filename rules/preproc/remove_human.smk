# vim: syntax=python expandtab
# Rules to filter human sequences from metagenomic reads
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path 

from snakemake.exceptions import WorkflowError

localrules:
    download_hg19,
    plot_proportion_human

if config["host_removal"]:
    hg19_path = Path(config["remove_human"]["hg19_path"])
    if not Path(hg19_path/"ref").exists():
        err_message = "Cannot find hg19 database for human sequence removal at: '{}'!\n".format(hg19_path)
        err_message += "Specify path to folder containing BBMap index of hg19 files in config.yaml.\n"
        err_message += "Run 'snakemake index_hg19' to download and create a BBMap index in '{dbdir}'".format(dbdir=DBDIR/"hg19")
        raise WorkflowError(err_message)

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    filtered_human = expand(str(OUTDIR/"filtered_human/{sample}_R{readpair}.filtered_human.fq.gz"),
            sample=SAMPLES,
            readpair=[1,2])
    proportion_human = str(OUTDIR/"filtered_human/proportion_human.pdf")
    all_outputs.extend(filtered_human)
    all_outputs.append(proportion_human)

    citations.add((
        "Bushnell, B. (2016).",
        "BBMap short read aligner.",
        "University of California, Berkeley, California.",
        "Available online at: http://sourceforge.net/projects/bbmap.",
    ))
else:
    filtered_human = expand(str(OUTDIR/"filtered_human/{sample}_R{readpair}.filtered_human.fq.gz"),
            sample=SAMPLES,
            readpair=[1,2])
    all_outputs.extend(filtered_human)
    


rule download_hg19:
    """Download masked hg19 from: 
    https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk"""
    output:
        OUTDIR/"hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        dbdir=DBDIR/"hg19"
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
        DBDIR/"hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    output:
        DBDIR/"hg19/ref/genome/1/chr1.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr2.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr3.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr4.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr5.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr6.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr7.chrom.gz",
        DBDIR/"hg19/ref/genome/1/info.txt",
        DBDIR/"hg19/ref/genome/1/scaffolds.txt.gz",
        DBDIR/"hg19/ref/genome/1/summary.txt",
        DBDIR/"hg19/ref/index/1/chr1-3_index_k13_c2_b1.block",
        DBDIR/"hg19/ref/index/1/chr1-3_index_k13_c2_b1.block2.gz",
        DBDIR/"hg19/ref/index/1/chr4-7_index_k13_c2_b1.block",
        DBDIR/"hg19/ref/index/1/chr4-7_index_k13_c2_b1.block2.gz",
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        dbdir=DBDIR/"hg19"
    shell:
        """
        bbmap.sh \
            ref={input} \
            path={params.dbdir}
        """


if config["host_removal"]:
    rh_config = config["remove_human"]
    rule remove_human:
        """Filter reads matching hg19. NB: requires about 16GB of memory."""
        input:
            read1=OUTDIR/"trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
            read2=OUTDIR/"trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
        output:
            read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
            read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
            human=OUTDIR/"filtered_human/{sample}_human.fq.gz",
        log:
            statsfile=str(LOGDIR/"remove_human/{sample}.statsfile.txt"),
            stderr=str(LOGDIR/"remove_human/{sample}.stderr.log"),
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            16
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


    rule plot_proportion_human:
        """Plot proportion of reads that matched hg19."""
        input:
            expand(str(LOGDIR/"remove_human/{sample}.statsfile.txt"), sample=SAMPLES)
        output:
            pdf=report(OUTDIR/"filtered_human/proportion_human.pdf",
                       category="Proportion human reads",
                       caption="../../report/proportion_human_reads.rst"),
            tsv=report(OUTDIR/"filtered_human/proportion_human.tsv",
                       category="Proportion human reads",
                       caption="../../report/tsv_proportion_human_reads.rst"),
        log:
            str(LOGDIR/"remove_human/proportion_human.log")
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            1
        params:
            unambigous=lambda _: "--unambigous" if rh_config["plot_unambigous"] else ""
        shell:
            """
            scripts/plot_proportion_human.py \
                {input} \
                {params.unambigous} \
                --outfile {output.pdf} \
                --table {output.tsv}
            """
else:
    rh_config = config["remove_human"]
    rule skip_remove_human:
        """Do not filter human sequences"""
        input:
            read1=OUTDIR/"trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
            read2=OUTDIR/"trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
        output:
            read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
            read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
        log:
            stderr=str(LOGDIR/"remove_human/{sample}.stderr.log"),
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            16
        params:
        shell:
            """
            ln -sv $(readlink -f {input.read1}) {output.read1} >> {log.stderr}
            ln -sv $(readlink -f {input.read2}) {output.read2} >> {log.stderr}
            """
