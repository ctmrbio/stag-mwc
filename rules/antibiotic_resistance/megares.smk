# Generic rules for detection of antibiotic resistance genes using MEGARes
from snakemake.exceptions import WorkflowError
import os.path

if not os.path.isdir(os.path.join(config["megares"]["db_path"], "ref")):
    err_message = "No MEGARes database found at: {}!\n".format(config["megares"]["db_path"])
    err_message += "Specify the DB path in the megares section of config.yaml.\n"
    err_message += "Run 'snakemake create_megares_index' to download and build a BBMap index in '{dbdir}/megares'\n".format(dbdir=config["dbdir"])
    err_message += "If you do not want to map reads against MEGARes for antibiotic resistance gene detection, set antibiotic_resistance: False in config.yaml"
    raise WorkflowError(err_message)

megares_outputs = expand("{outdir}/megares/{sample}.{output_type}",
        outdir=outdir,
        sample=SAMPLES,
        output_type=("sam.gz", "mapped_reads.fq.gz", "mhist.txt", "covstats.txt", "rpkm.txt"))
all_outputs.extend(megares_outputs)

rule download_megares:
    """Download MEGARes database."""
    output:
        config["dbdir"]+"/megares/megares_annotations_v1.01.csv",
        config["dbdir"]+"/megares/megares_database_v1.01.fasta",
        config["dbdir"]+"/megares/megares_to_external_header_mappings_v1.01.tsv",
    shadow:
        "shallow"
    params:
        dbdir=config["dbdir"]+"/megares"
    shell:
        """
        cd {params.dbdir}
        wget http://megares.meglab.org/download/megares_v1.01.zip \
        && \
        unzip megares_v1.01.zip \
        && \
        mv megares_v1.01/* . \
        && \
        rm -rfv megares_v1.01 megares_v1.01.zip
        """


rule create_megares_index:
    """Create BBMap index for MEGARes."""
    input:
        fasta=config["dbdir"]+"/megares/megares_database_v1.01.fasta"
    output:
        config["dbdir"]+"/megares/ref/genome/1/chr1.chrom.gz",
        config["dbdir"]+"/megares/ref/genome/1/info.txt",
        config["dbdir"]+"/megares/ref/genome/1/scaffolds.txt.gz",
        config["dbdir"]+"/megares/ref/genome/1/summary.txt",
        config["dbdir"]+"/megares/ref/index/1/chr1_index_k13_c8_b1.block",
        config["dbdir"]+"/megares/ref/index/1/chr1_index_k13_c8_b1.block2.gz",
    shadow:
        "shallow"
    conda:
        "../../envs/bbmap.yaml"
    params:
        dbdir=config["dbdir"]+"/megares"
    shell:
        """
        bbmap.sh ref={input} path={params.dbdir}
        """


megares_config = config["megares"]
rule bbmap_to_megares:
    """BBMap to MEGARes."""
    input:
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        sam=config["outdir"]+"/megares/{sample}.sam.gz",
        mapped_reads=config["outdir"]+"/megares/{sample}.mapped_reads.fq.gz",
        covstats=config["outdir"]+"/megares/{sample}.covstats.txt",
        rpkm=config["outdir"]+"/megares/{sample}.rpkm.txt",
        mhist=config["outdir"]+"/megares/{sample}.mhist.txt",
    log:
        stdout=config["outdir"]+"/logs/megares/{sample}.bbmap.stdout.log",
        stderr=config["outdir"]+"/logs/megares/{sample}.bbmap.statsfile.txt"
    shadow:
        "shallow"
    conda:
        "../../envs/bbmap.yaml"
    threads:
        8
    params:
        db_path=megares_config["db_path"],
        min_id=megares_config["min_id"],
        max_indel=megares_config["max_indel"],
        extra=megares_config["extra"],
    shell:
        """
        bbmap.sh \
            threads={threads} \
            minid={params.min_id} \
            maxindel={params.max_indel} \
            path={params.db_path} \
            in1={input.read1} \
            in2={input.read2} \
            out={output.sam} \
            outm={output.mapped_reads} \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            mhist={output.mhist} \
            {params.extra} \
            > {log.stdout} \
            2> {log.stderr}
        """
