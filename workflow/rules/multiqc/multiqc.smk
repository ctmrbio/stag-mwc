# vim: syntax=python expandtab
# MultiQC

if config["multiqc_report"]:
    citations.add(publications["MultiQC"])

    mqc_config = config["multiqc"]
    rule multiqc:
        input:
            all_outputs
        output:
            report=report(OUTDIR/"multiqc/multiqc_report.html",
                category="Sequencing data quality",
                caption="../../report/multiqc.rst"),
        log:
           LOGDIR/"multiqc/multiqc.log"
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            config["containers"]["multiqc"]
        threads: 1
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
    all_outputs.append(OUTDIR/"multiqc/multiqc_report.html")

