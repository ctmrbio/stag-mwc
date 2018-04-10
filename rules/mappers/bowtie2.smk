# Generic rules for alignment of reads to a reference database using Bowtie2
from snakemake.exceptions import WorkflowError
import os.path

bt2_db_extensions = (".1.bt2", ".1.bt2l")
if not any([os.path.isfile(config["bowtie2"]["db_prefix"]+ext) for ext in bt2_db_extensions]):
    err_message = "Bowtie2 index not found at: '{}'\n".format(config["bowtie2"]["db_prefix"])
    err_message += "Check path in config setting 'bowtie2:db_prefix'."
    raise WorkflowError(err_message)
bt2_db_name = os.path.basename(config["bowtie2"]["db_prefix"])


rule bowtie2:
    """Align reads using Bowtie2."""
    input:
        sample=[config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
                config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz"]
    output:
        config["outdir"]+"/bowtie2/{db_name}/{{sample}}.bam".format(db_name=bt2_db_name)
    log:
        config["outdir"]+"/logs/bowtie2/{db_name}/{{sample}}.log".format(db_name=bt2_db_name)
    params:
        index=config["bowtie2"]["db_prefix"],
        extra=config["bowtie2"]["extra"],
    threads:
        8
    wrapper:
        "0.23.1/bio/bowtie2/align"


rule bowtie2_mapping_stats:
    """Summarize bowtie2 mapping statistics."""
    input:
        bam=config["outdir"]+"/bowtie2/{dbname}/{{sample}}.bam".format(dbname=bt2_db_name)
    output:
        covstats=config["outdir"]+"/bowtie2/{dbname}/{{sample}}.covstats.txt".format(dbname=bt2_db_name),
        rpkm=config["outdir"]+"/bowtie2/{dbname}/{{sample}}.rpkm.txt".format(dbname=bt2_db_name)
    log:
        config["outdir"]+"/logs/bowtie2/{dbname}/{{sample}}.pileup.log".format(dbname=bt2_db_name)
    shadow:
        "shallow"
    conda:
        "../../envs/bbmap.yaml"
    shell:
        """
        pileup.sh \
            in={input.bam} \
            out={output.covstats} \
            rpkm={output.rpkm} \
            > {log}
        """
