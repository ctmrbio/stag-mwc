# vim: syntax=python expandtab
# Assess sequencing depth of sample using BBCountUnique from the BBMap suite.

if config["naive"]["assess_depth"]:
    # Add final output files from this module to 'all_outputs' from
    # the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
    bcu_output = expand(str(OUTDIR/"bbcountunique/{sample}.{output_type}"),
            sample=SAMPLES,
            output_type=["bbcountunique.txt", "bbcountunique.pdf"])
    all_outputs.extend(bcu_output)

    citations.add(publications["BBMap"])

    rule bbcountunique:
        """Assess sequencing depth using BBCountUnique."""
        input:
            INPUT_read1,
        output:
            txt=OUTDIR/"bbcountunique/{sample}.bbcountunique.txt",
        log:
            stdout=LOGDIR/"bbcountunique/{sample}.bbcountunique.stdout.log",
            stderr=LOGDIR/"bbcountunique/{sample}.bbcountunique.stderr.log",
        shadow: 
            "shallow"
        threads: 2
        conda:
            "../../envs/stag-mwc.yaml",
        container:
            config["containers"]["bbmap"]
        params:
            interval=config["bbcountunique"]["interval"]
        shell:
            """
            bbcountunique.sh \
              in={input} \
              out={output.txt} \
              interval={params.interval} \
              > {log.stdout} \
              2> {log.stderr}
            """

    rule plot_bbcountunique:
        """Assess sequencing depth using BBCountUnique."""
        input:
            txt=OUTDIR/"bbcountunique/{sample}.bbcountunique.txt",
        output:
            pdf=report(OUTDIR/"bbcountunique/{sample}.bbcountunique.pdf",
                       caption="../../report/bbcountunique.rst",
                       category="Sequencing depth")
        log:
            stdout=LOGDIR/"bbcountunique/{sample}.plot_bbcountunique.stdout.log",
            stderr=LOGDIR/"bbcountunique/{sample}.plot_bbcountunique.stderr.log",
        shadow: 
            "shallow"
        threads: 1
        conda:
            "../../envs/stag-mwc.yaml",
        container:
            config["containers"]["stag"]
        shell:
            """
            workflow/scripts/plot_bbcountunique.py \
              {input.txt} \
              {output.pdf} \
              >> {log.stdout} \
              2>> {log.stderr}
            """
