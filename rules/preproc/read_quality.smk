# vim: syntax=python expandtab
# Read pre-processing 
#TODO: Remove superfluous str conversions of paths in expand and log statements
#      when Snakemake is pathlib compatible.

if config["qc_reads"]:
    # Add final output files from this module to 'all_outputs' from
    # the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
    trimmed_qc = expand(str(OUTDIR/"fastp/{sample}_R{readpair}.qc.fq.gz"),
            sample=SAMPLES,
            readpair=[1, 2])
    all_outputs.extend(trimmed_qc)

    citations.add((
        "Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu (2018).",
        "fastp: an ultra-fast all-in-one FASTQ preprocessor.",
        "Bioinformatics, Volume 34, Issue 17, 01 September 2018, Pages i884–i890.",
        "https://doi.org/10.1093/bioinformatics/bty560"
    ))


    fastp_config = config["fastp"]
    rule trim_adapters_quality:
        input:
            read1=INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="1"),
            read2=INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="2")
        output:
            read1=OUTDIR/"fastp/{sample}_R1.qc.fq.gz",
            read2=OUTDIR/"fastp/{sample}_R2.qc.fq.gz",
            json=LOGDIR/"fastp/{sample}.fastp.json",
            html=LOGDIR/"fastp/{sample}.fastp.html",
        log:
            stdout=str(LOGDIR/"fastp/{sample}.stdout.log"),
            stderr=str(LOGDIR/"fastp/{sample}.stderr.log"),
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            4
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
