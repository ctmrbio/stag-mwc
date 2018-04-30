# Generic rules for detection of antibiotic resistance genes using MEGARes
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
    download_megares

megares_db_path = Path(config["megares"]["db_path"])
if not Path(megares_db_path/"ref").exists():
    err_message = "No MEGARes database found at: '{}'!\n".format(megares_db_path)
    err_message += "Specify the DB path in the megares section of config.yaml.\n"
    err_message += "Run 'snakemake create_megares_index' to download and build a BBMap index in '{dbdir}'\n".format(dbdir=DBDIR/"megares")
    err_message += "If you do not want to map reads against MEGARes for antibiotic resistance gene detection, set antibiotic_resistance: False in config.yaml"
    raise WorkflowError(err_message)

megares_outputs = expand(str(OUTDIR/"megares/{sample}.{output_type}"),
        sample=SAMPLES,
        output_type=("sam.gz", "mapped_reads.fq.gz", "mhist.txt", "covstats.txt", "rpkm.txt"))
all_outputs.extend(megares_outputs)

rule download_megares:
    """Download MEGARes database."""
    output:
        DBDIR/"megares/megares_annotations_v1.01.csv",
        DBDIR/"megares/megares_database_v1.01.fasta",
        DBDIR/"megares/megares_to_external_header_mappings_v1.01.tsv",
    log:
        str(LOGDIR/"megares/megares.download.log")
    shadow:
        "shallow"
    params:
        dbdir=DBDIR/"megares"
    shell:
        """
        cd {params.dbdir}
        wget http://megares.meglab.org/download/megares_v1.01.zip \
            > {log} \
        && \
        unzip megares_v1.01.zip \
            >> {log} \
        && \
        mv megares_v1.01/* . \
        && \
        rm -rfv megares_v1.01 megares_v1.01.zip \
            >> {log}
        """


rule create_megares_index:
    """Create BBMap index for MEGARes."""
    input:
        fasta=DBDIR/"megares/megares_database_v1.01.fasta"
    output:
        DBDIR/"megares/ref/genome/1/chr1.chrom.gz",
        DBDIR/"megares/ref/genome/1/info.txt",
        DBDIR/"megares/ref/genome/1/scaffolds.txt.gz",
        DBDIR/"megares/ref/genome/1/summary.txt",
        DBDIR/"megares/ref/index/1/chr1_index_k13_c8_b1.block",
        DBDIR/"megares/ref/index/1/chr1_index_k13_c8_b1.block2.gz",
    log:
        str(LOGDIR/"megares/megares.bbmap_index.log")
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        dbdir=DBDIR/"megares"
    shell:
        """
        bbmap.sh ref={input} path={params.dbdir} > {log}
        """


megares_config = config["megares"]
rule bbmap_to_megares:
    """BBMap to MEGARes."""
    input:
        read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        sam=OUTDIR/"megares/{sample}.sam.gz",
        mapped_reads=OUTDIR/"megares/{sample}.mapped_reads.fq.gz",
        covstats=OUTDIR/"megares/{sample}.covstats.txt",
        rpkm=OUTDIR/"megares/{sample}.rpkm.txt",
        mhist=OUTDIR/"megares/{sample}.mhist.txt",
    log:
        stdout=str(LOGDIR/"megares/{sample}.bbmap.stdout.log"),
        stderr=str(LOGDIR/"megares/{sample}.bbmap.statsfile.txt"),
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
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
