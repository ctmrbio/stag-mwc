# vim: syntax=python expandtab
# MultiQC

if config["multiqc_report"]:
    citations.add((
        "Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller.",
        "MultiQC: Summarize analysis results for multiple tools and samples in a single report.",
        "Bioinformatics (2016).",
        "doi: 10.1093/bioinformatics/btw354.",
    ))

    mqc_config = config["multiqc"]
    rule multiqc:
        input:
            all_outputs
        output:
            report=report(f"{OUTDIR}/multiqc/multiqc_report.html",
                category="Sequencing data quality",
                caption="../../report/multiqc.rst"),
        log:
           f"{LOGDIR}/multiqc/multiqc.log"
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            1
        params:
            extra=mqc_config["extra"],
        shell:
            """
            multiqc {OUTDIR} \
                --filename {output.report} \
                --force \
                2> {log}
            """

    # Appended after the rule definition to avoid circular dependency
    all_outputs.append(f"{OUTDIR}/multiqc/multiqc_report.html")
