# vim: syntax=python expandtab
# Assess sequencing depth of sample using BBCountUnique from the BBMap suite.
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.

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
            pdf=report(OUTDIR/"bbcountunique/{sample}.bbcountunique.pdf",
                       caption="../../report/bbcountunique.rst",
                       category="Sequencing depth")
        log:
            stdout=str(LOGDIR/"bbcountunique/{sample}.bbcountunique.stdout.log"),
            stderr=str(LOGDIR/"bbcountunique/{sample}.bbcountunique.stderr.log"),
        shadow: 
            "shallow"
        threads:
            cluster_config["bbcountunique"]["n"] if "bbcountunique" in cluster_config else 2
        conda:
            "../../envs/stag-mwc.yaml",
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
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

            workflow/scripts/plot_bbcountunique.py \
              {output.txt} \
              {output.pdf} \
              >> {log.stdout} \
              2>> {log.stderr}
            """
