# vim: syntax=python expandtab
# Assembly of metagenomic samples
from pathlib import Path

from snakemake.exceptions import WorkflowError

a_conf = config["assembly_"]
b_conf = config["binning_"]

localrules:


if config["assembly"]["assemble"]:
    assemblies = expand(f"{OUTDIR}/assembly/{a_conf['assembler']}/{{sample}}/{{sample}}.contigs.fa", sample=SAMPLES)
    all_outputs.extend(assemblies)

    citations.add(publications["MEGAHIT"])


rule megahit:
    f"""Metagenomic assembly with MEGAHIT."""
    input:
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
    output:
        contigs=f"{OUTDIR}/assembly/megahit/{{sample}}/{{sample}}.contigs.fa",
    log:
        stdout=f"{LOGDIR}/megahit/{{sample}}.stdout.log",
        stderr=f"{LOGDIR}/megahit/{{sample}}.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    container:
        "docker://vout/megahit:release-v1.2.9"
    threads: 20
    params:
        outdir=lambda w: f"{OUTDIR}/assembly/megahit/{w.sample}",
        min_contig_length=a_conf["min_contig_length"],
    shell:
        """
        rm -rfv {params.outdir}
        megahit \
            -1 {input.read1} \
            -2 {input.read2} \
            --out-dir {params.outdir} \
            --out-prefix {wildcards.sample} \
            --num-cpu-threads {threads} \
            --memory {resources.mem_mb} \
            --min-contig-len {params.min_contig_length} \
            2> {log.stderr} \
            > {log.stdout}
        """

