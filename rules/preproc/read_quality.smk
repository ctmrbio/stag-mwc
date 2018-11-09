# vim: syntax=python expandtab
# MWC read pre-processing rules
#TODO: Remove superfluous str conversions of paths in expand and log statements
#      when Snakemake is pathlib compatible.

if config["qc_reads"]:
    # Add final output files from this module to 'all_outputs' from
    # the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
    fastqc_output = expand(str(OUTDIR/"fastqc/{sample}_R{readpair}_fastqc.{ext}"),
            sample=SAMPLES,
            readpair=[1, 2],
            ext=["zip", "html"])
    trimmed_qa = expand(str(OUTDIR/"trimmed_qa/{sample}_R{readpair}.trimmed_qa.fq.gz"),
            sample=SAMPLES,
            readpair=[1, 2])
    all_outputs.extend(fastqc_output)
    all_outputs.extend(trimmed_qa)

rule fastqc:
    input:
        INPUTDIR/config["input_fn_pattern"]
    output:
        html=OUTDIR/"fastqc/{sample}_R{readpair}_fastqc.html",
        zip=OUTDIR/"fastqc/{sample}_R{readpair}_fastqc.zip",
    log:
        str(LOGDIR/"fastqc/{sample}_R{readpair}_fastq.log")
    shadow: 
        "shallow"
    wrapper:
        "0.22.0/bio/fastqc"


bbduk_config = config["bbduk"]
rule trim_adapters_quality:
    input:
        read1=INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="1"),
        read2=INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="2")
    output:
        read1=OUTDIR/"trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
        read2=OUTDIR/"trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
    log:
        stdout=str(LOGDIR/"bbduk/{sample}.stdout.log"),
        stderr=str(LOGDIR/"bbduk/{sample}.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        4
    params:
        minlen=bbduk_config["minlen"],
        qtrim=bbduk_config["qtrim"],
        trimq=bbduk_config["trimq"],
        ktrim=bbduk_config["ktrim"],
        k=bbduk_config["k"],
        mink=bbduk_config["mink"],
        hdist=bbduk_config["hdist"],
        trimbyoverlap=bbduk_config["trimbyoverlap"],
        trimpairsevenly=bbduk_config["trimpairsevenly"],
    shell:
        """
        bbduk.sh \
            in1={input.read1} \
            in2={input.read2} \
            out1={output.read1} \
            out2={output.read2} \
            ref=adapters \
            threads={threads} \
            minlen={params.minlen} \
            qtrim={params.qtrim} \
            trimq={params.trimq} \
            ktrim={params.ktrim} \
            k={params.k} \
            mink={params.mink} \
            hdist={params.hdist} \
            {params.trimbyoverlap} \
            {params.trimpairsevenly} \
            > {log.stdout} \
            2> {log.stderr}
        """
