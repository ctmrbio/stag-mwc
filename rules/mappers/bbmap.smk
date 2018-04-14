# Rules for generic read mapping using BBMap
from snakemake.exceptions import WorkflowError
import os.path

if not os.path.isdir(config["bbmap"]["db_path"]):
    err_message = "BBMap index not found at: '{}'\n".format(config["bbmap"]["db_path"])
    err_message += "Check path in config setting 'bbmap:db_path'.\n"
    err_message += "If you want to skip mapping with BBMap, set mappers:bbmap:False in config.yaml."
    raise WorkflowError(err_message)

bbmap_config = config["bbmap"]
bbmap_output_folder = config["outdir"]+"/bbmap/{db_name}/".format(db_name=bbmap_config["db_name"])
rule bbmap:
    """BBMap"""
    input:
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        sam=bbmap_output_folder+"{sample}.sam.gz",
        covstats=bbmap_output_folder+"{sample}.covstats.txt",
        rpkm=bbmap_output_folder+"{sample}.rpkm.txt",
    log:
        stdout=config["outdir"]+"/logs/bbmap/{sample}.bbmap.stdout.log",
        stderr=config["outdir"]+"/logs/bbmap/{sample}.bbmap.statsfile.txt"
    shadow:
        "shallow"
    conda:
        "../../envs/bbmap.yaml"
    threads:
        8
    params:
        db_path=bbmap_config["db_path"],
        min_id=bbmap_config["min_id"],
        extra=bbmap_config["extra"],
    shell:
        """
        bbmap.sh \
            threads={threads} \
            minid={params.min_id} \
            path={params.db_path} \
            in1={input.read1} \
            in2={input.read2} \
            out={output.sam} \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            {params.extra} \
            > {log.stdout} \
            2> {log.stderr}
        """
