# vim: syntax=python expandtab
# Read pre-processing 
# TODO: Remove superfluous str conversions of paths in expand and log statements
#      when Snakemake is pathlib compatible.

fastp_config = config["fastp"]
if config["qc_reads"]:
    # Add final output files from this module to 'all_outputs' from
    # the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
    trimmed_qc = expand(str(OUTDIR/"fastp/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1, 2])

    if fastp_config["keep_output"]:
        all_outputs.extend(trimmed_qc)

    citations.add(publications["fastp"])

    rule fastp:
        input:
            read1=INPUT_read1,
            read2=INPUT_read2,
        output:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz" if fastp_config["keep_output"] else temp(OUTDIR/"fastp/{sample}_1.fq.gz"),
            read2=OUTDIR/"fastp/{sample}_2.fq.gz" if fastp_config["keep_output"] else temp(OUTDIR/"fastp/{sample}_2.fq.gz"),
            json=LOGDIR/"fastp/{sample}.fastp.json",
            html=LOGDIR/"fastp/{sample}.fastp.html",
        log:
            stdout=str(LOGDIR/"fastp/{sample}.stdout.log"),
            stderr=str(LOGDIR/"fastp/{sample}.stderr.log"),
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads:
            cluster_config["fastp"]["n"] if "fastp" in cluster_config else 4
        params:
            extra=fastp_config["extra"],
        shell:
            """
            fastp \
                --in1 {input.read1} \
                --in2 {input.read2} \
                --out1 {output.read1} \
                --out2 {output.read2} \
                --json {output.json} \
                --html {output.html} \
                --thread {threads} \
                {params.extra} \
                > {log.stdout} \
                2> {log.stderr}
            """

else:
    trimmed_qc = expand(str(OUTDIR/"fastp/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1, 2])

    if fastp_config["keep_output"]:
        all_outputs.extend(trimmed_qc)

    localrules:
        skip_fastp,

    rule skip_fastp:
        input:
            read1=INPUT_read1,
            read2=INPUT_read2,
        output:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        log:
            stderr=str(LOGDIR/"fastp/{sample}.stderr.log"),
        shell:
            """
            ln -sv $(readlink -f {input.read1}) {output.read1} >> {log.stderr}
            ln -sv $(readlink -f {input.read2}) {output.read2} >> {log.stderr}
            """
