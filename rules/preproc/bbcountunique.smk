# vim: syntax=python expandtab
# Assess sequencing depth of sample using BBCountUnique from the BBMap suite.
import os.path

# Add final output files from this module to 'all_outputs' from
# the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
bcu_output = expand("{outdir}/bbcountunique/{sample}.{output_type}",
        outdir=config["outdir"],
        sample=SAMPLES,
        output_type=["bbcountunique.txt", "bbcountunique.pdf"])
all_outputs.extend(bcu_output)

rule bbcountunique:
    """Assess sequencing depth using BBCountUnique."""
    input:
        os.path.join(config["inputdir"], config["input_fn_pattern"]).format(sample="{sample}", readpair="1")
    output:
        txt=config["outdir"]+"/bbcountunique/{sample}.bbcountunique.txt",
        pdf=config["outdir"]+"/bbcountunique/{sample}.bbcountunique.pdf",
    log:
        stdout=config["outdir"]+"/logs/bbcountunique/{sample}.bbcountunique.stdout.log",
        stderr=config["outdir"]+"/logs/bbcountunique/{sample}.bbcountunique.stderr.log",
    shadow: 
        "shallow"
    threads:
        2
    conda:
        "../../envs/stag-mwc.yaml",
    params:
        interval=config["bbcountunique"]["interval"]
    shell:
        """
        bbcountunique.sh in={input} out={output.txt} interval={params.interval} > {log.stdout} 2> {log.stderr} \
        && \
        scripts/plot_bbcountunique.py {output.txt} {output.pdf} >> {log.stdout} 2>> {log.stderr}
        """
