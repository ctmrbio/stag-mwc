# vim: syntax=python expandtab
# MWC read pre-processing rules

# Add final output files from this module to 'all_outputs' from
# the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
fastqc_output = expand(outdir/"fastqc/{sample}_R{readpair}_fastqc.{ext}",
        sample=SAMPLES,
        readpair=[1, 2],
        ext=["zip", "html"])
trimmed_qa = expand(outdir/"trimmed_qa/{sample}_R{readpair}.trimmed_qa.fq.gz",
        sample=SAMPLES,
        readpair=[1, 2])
all_outputs.extend(fastqc_output)
all_outputs.extend(trimmed_qa)

rule fastqc:
    input:
        inputdir/config["input_fn_pattern"]
    output:
        html=outdir/"fastqc/{sample}_R{readpair}_fastqc.html",
        zip=outdir/"fastqc/{sample}_R{readpair}_fastqc.zip",
    shadow: 
        "shallow"
    wrapper:
        "0.22.0/bio/fastqc"


bbduk_config = config["bbduk"]
rule trim_adapters_quality:
    input:
        read1=inputdir/config["input_fn_pattern"].format(sample="{sample}", readpair="1"),
        read2=inputdir/config["input_fn_pattern"].format(sample="{sample}", readpair="2")
    output:
        read1=outdir/"trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
        read2=outdir/"trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
    log:
        stdout=logdir/"bbduk/{sample}.stdout.log",
        stderr=logdir/"bbduk/{sample}.stderr.log",
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
