# vim: syntax=python expandtab
# MWC read pre-processing rules

rule fastqc:
    input:
        "input/{sample}_R{readpair}.fastq.gz"
    output:
        html=config["outdir"]+"/fastqc/{sample}_R{readpair}.html",
        zip=config["outdir"]+"/fastqc/{sample}_R{readpair}.zip",
    shadow: 
        "shallow"
    wrapper:
        "0.22.0/bio/fastqc"


rule trim_adapters_quality:
    input:
        read1="input/{sample}_R1.fastq.gz",
        read2="input/{sample}_R2.fastq.gz",
    output:
        read1=config["outdir"]+"/trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
        read2=config["outdir"]+"/trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
    conda: "../../envs/bbmap.yaml"
    shadow: "shallow"
    shell:
        """
        bbduk.sh \
            in1={input.read1} \
            in2={input.read2} \
            out1={output.read1} \
            out2={output.read2} \
            ref=adapters \
            threads={threads} \
            minlen={config[bbduk_minlen]} \
            qtrim={config[bbduk_qtrim]} \
            trimq={config[bbduk_trimq]} \
            ktrim={config[bbduk_ktrim]} \
            k={config[bbduk_k]} \
            mink={config[bbduk_mink]} \
            hdist={config[bbduk_hdist]} \
            {config[bbduk_trimbyoverlap]} \
            {config[bbduk_trimpairsevenly]} 
        """
