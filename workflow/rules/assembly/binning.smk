# vim: syntax=python expandtab
# Assembly and binning of metagenomic assemblies
from pathlib import Path

from snakemake.exceptions import WorkflowError

a_conf = config["assembly_"]
b_conf = config["binning_"]

localrules:


if config["assembly"]["binning"]:
    concatenated_contigs = f"{OUTDIR}/assembly/{a_conf['assembler']}/concatenated_contigs.fa.gz"
    bam_files = expand(f"{OUTDIR}/binning/{{sample}}.bam", sample=SAMPLES)
    bins = expand(f"{OUTDIR}/binning/{{sample}}/bins", sample=SAMPLES)

    all_outputs.append(concatenated_contigs)
    all_outputs.extend(bam_files)
    all_outputs.extend(bins)

    citations.add(publications["BBMap"])
    citations.add(publications["Vamb"])


rule concatenate_index_contigs:
    """Concatenate and index contigs before mapping and binning."""
    input:
        contigs=expand(f"{OUTDIR}/assembly/{a_conf['assembler']}/{{sample}}/{{sample}}.contigs.fa", sample=SAMPLES),
    output:
        concatenated=f"{OUTDIR}/assembly/{a_conf['assembler']}/concatenated_contigs.fa.gz",
    log:
        stdout=f"{LOGDIR}/binning/concatenate_and_index.stdout",
        stderr=f"{LOGDIR}/binning/concatenate_and_index.stderr",
    conda:
        "../../envs/assembly.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:assembly"+singularity_branch_tag
    threads: 8
    params:
        min_contig_length=b_conf["min_contig_length"],
        index_dir=f"{OUTDIR}/assembly/{a_conf['assembler']}/",
    shell:
        """
        concatenate.py \
            -m {params.min_contig_length} \
            {output.concatenated} \
            {input.contigs} \
            > {log.stdout} \
            2> {log.stderr}

        bbmap.sh \
            ref={output.concatenated} \
            path={params.index_dir} \
            >> {log.stdout} \
            2>> {log.stderr}
        """


rule map_to_contigs:
    """Map reads to concatenated contigs."""
    input:
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
        contigs=rules.concatenate_index_contigs.output.concatenated,
    output:
        sam=f"{OUTDIR}/binning/{{sample}}.sam.gz",
        bamscript=f"{OUTDIR}/binning/{{sample}}.bamscript.sh",
        bam=f"{OUTDIR}/binning/{{sample}}.bam",
    log:
        stdout=f"{LOGDIR}/binning/{{sample}}.bbmap.stdout",
        stderr=f"{LOGDIR}/binning/{{sample}}.bbmap.stderr",
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:assembly"+singularity_branch_tag
    threads: 20
    shell:
        """
        bbmap.sh \
            in1={input.read1} \
            in2={input.read2} \
            path=$(dirname {input.contigs}) \
            out={output.sam} \
            bamscript={output.bamscript} \
            > {log.stdout} \
            2> {log.stdout} 

        sed -i 's/_sorted//g' {output.bamscript}

        ./{output.bamscript} >> {log.stdout} 2>> {log.stderr}
        """


rule vamb:
    """Metagenomic binning.""" 
    input:
        contigs=rules.concatenate_index_contigs.output.concatenated,
        bam_files=expand(f"{OUTDIR}/binning/{{sample}}.bam", sample=SAMPLES),
    output:
        bins=directory(f"{OUTDIR}/binning/{{sample}}/bins"),
    log:
        stdout=f"{LOGDIR}/binning/{{sample}}.vamb.stdout",
        stderr=f"{LOGDIR}/binning/{{sample}}.vamb.stderr",
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:assembly"+singularity_branch_tag
    threads: 20
    params:
        minfasta=b_conf["minfasta"],
    shell:
        """
        vamb \
            -o C \
            --outdir {output.bins} \
            --fasta {input.contigs} \
            --bamfiles {input.bam_files} \
            > {log.stdout} \
            2> {log.stderr}
        """

