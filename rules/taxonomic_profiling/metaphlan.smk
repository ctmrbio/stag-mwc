# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using MetaPhlAn
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    plot_metaphlan_heatmap,
    create_metaphlan_krona_plots,
    metaphlan_krona,

mpa_config = config["metaphlan"]
if config["taxonomic_profile"]["metaphlan"]:
    if not mpa_config["bt2_db_dir"] or not Path(mpa_config["bt2_db_dir"]).exists():
        err_message = "No MetaPhlAn database dir found at: '{}'!\n".format(mpa_config["bt2_db_dir"])
        err_message += "Specify relevant paths in the metaphlan section of config.yaml.\n"
        err_message += "If you do not want to run MetaPhlAn for taxonomic profiling, set metaphlan: False in config.yaml"
        raise WorkflowError(err_message)

    # Add MetaPhlAn output files to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    mpa_outputs = expand(f"{OUTDIR}/metaphlan/{{sample}}.{{output_type}}",
            sample=SAMPLES,
            output_type=("bowtie2.bz2", "metaphlan.txt", "metaphlan.krona"))
    mpa_combined = expand(f"{OUTDIR}/metaphlan/all_samples.metaphlan.{{ext}}",
            ext=("txt", "krona.html"))
    mpa_plot = f"{OUTDIR}/metaphlan/all_samples.{mpa_config['heatmap']['level']}_top{mpa_config['heatmap']['topN']}.pdf"
    mpa_area_plot = f"{OUTDIR}/metaphlan/area_plot.metaphlan.pdf"
    all_outputs.extend(mpa_outputs)
    all_outputs.extend(mpa_combined)
    all_outputs.append(mpa_plot)
    all_outputs.append(mpa_area_plot)

    citations.add(publications["MetaPhlAn"])
    citations.add(publications["Krona"])


rule metaphlan:
    """Taxonomic profiling using MetaPhlAn."""
    input:
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
    output:
        bt2_out=f"{OUTDIR}/metaphlan/{{sample}}.bowtie2.bz2",
        mpa_out=f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.txt",
    log:
        stdout=f"{LOGDIR}/metaphlan/{{sample}}.metaphlan.stdout.log",
        stderr=f"{LOGDIR}/metaphlan/{{sample}}.metaphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc-biobakery"
    threads:
        cluster_config["metaphlan"]["n"] if "metaphlan" in cluster_config else 5
    params:
        bt2_db_dir=mpa_config["bt2_db_dir"],
        bt2_index=mpa_config["bt2_index"],
        extra=mpa_config["extra"],
    shell:
        """
        metaphlan \
            --input_type fastq \
            --nproc {threads} \
            --sample_id {wildcards.sample} \
            --bowtie2out {output.bt2_out} \
            --bowtie2db {params.bt2_db_dir} \
            --index {params.bt2_index} \
            {input.read1},{input.read2} \
            {output.mpa_out} \
            {params.extra} \
            > {log.stdout} \
            2> {log.stderr}
        """


rule metaphlan_krona:
    """Convert MPA2 output to Krona input"""
    input:
        mpa_out=f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.txt",
    output:
        krona=f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.krona",
    log:
        f"{LOGDIR}/metaphlan/{{sample}}.metaphlan_krona.log",
    shadow:
        "shallow"
    shell:
        """
        # This command is broken in MetaPhlAn v2.96.1 when running with old DB
        #metaphlankrona.py \
        #    --profile {input.mpa_out} \
        #    --krona {output.krona} \
        #    2>&1 > {log}
       
        set +o pipefail  # Small samples can produce empty output files failing the pipeline
        gsed '/#/d' {input.mpa_out} \
            | grep -E "s__|unclassified" \
            | cut -f1,3 \
            | awk '{{print $2,"\t",$1}}' \
            | gsed 's/|\w__/\t/g' \
            | gsed 's/k__//' \
            > {output.krona} \
            2> {log}
        """

rule combine_metaphlan_tables:
    """Combine metaphlan outputs into a large table."""
    input:
        expand(f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.txt", sample=SAMPLES)
    output:
        txt=report(f"{OUTDIR}/metaphlan/all_samples.metaphlan.txt",
                   category="Taxonomic profiling",
                   caption="../../report/metaphlan_table.rst"),
    log:
        f"{LOGDIR}/metaphlan/combine_metaphlan_tables.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc-biobakery"
    threads:
        1
    shell:
        """
        merge_metaphlan_tables.py {input} > {output.txt}
        gsed --in-place 's/\.metaphlan//g' {output.txt} 
        """


rule metaphlan_area_plot:
    input:
        f"{OUTDIR}/metaphlan/all_samples.metaphlan.txt"
    output:
        report(f"{OUTDIR}/metaphlan/area_plot.metaphlan.pdf",
            category="Taxonomic profiling",
            caption="../../report/area_plot.rst")
    log:
        f"{LOGDIR}/metaphlan/area_plot.log"
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/area_plot.py \
            --table {input} \
            --output {output} \
            --mode metaphlan \
            2>&1 > {log}
        """


rule plot_metaphlan_heatmap:
    """Plot MetaPhlAn heatmap."""
    input:
        f"{OUTDIR}/metaphlan/all_samples.metaphlan.txt",
    output:
        pdf=report(f"{OUTDIR}/metaphlan/all_samples.{{level}}_top{{topN}}.pdf",
                   category="Taxonomic profiling",
                   caption="../../report/metaphlan.rst"),
    log:
        f"{LOGDIR}/metaphlan/plot_metaphlan_heatmap.{{level}}_top{{topN}}.log",
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    threads:
        1
    params:
        outfile_prefix=lambda w: f"{OUTDIR}/metaphlan/all_samples",
        pseudocount=mpa_config["heatmap"]["pseudocount"],
        colormap=mpa_config["heatmap"]["colormap"],
        method=mpa_config["heatmap"]["method"],
        metric=mpa_config["heatmap"]["metric"],
        extra=mpa_config["heatmap"]["extra"],
    shell:
        """
        scripts/plot_metaphlan_heatmap.py \
            --outfile-prefix {params.outfile_prefix} \
            --level {wildcards.level} \
            --topN {wildcards.topN} \
            --pseudocount {params.pseudocount} \
            --colormap {params.colormap} \
            --method {params.method} \
            --metric {params.metric} \
            --force \
            {input} \
            2> {log}
        """

rule create_metaphlan_krona_plots:
    input:
        expand(f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.krona", sample=SAMPLES)
    output:
        html=report(OUTDIR/"metaphlan/all_samples.metaphlan.krona.html",
                    category="Taxonomic profiling",
                    caption="../../report/metaphlan_krona.rst"),
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc-biobakery"
    threads:
        1
    shell:
        """
        ktImportText \
            -o {output.html} \
            {input}
        """
