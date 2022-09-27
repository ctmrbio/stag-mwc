# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using MetaPhlAn

from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
    combine_metaphlan_tables,
    metaphlan_krona,
    metaphlan_area_plot,
    plot_metaphlan_heatmap,
    create_metaphlan_krona_plots,
    metaphlan_outputs,

mpa_config = config["metaphlan"]
if config["taxonomic_profile"]["metaphlan"] or config["functional_profile"]["humann"]:
    if not mpa_config["bt2_db_dir"] or not Path(mpa_config["bt2_db_dir"]).exists():
        err_message = "No MetaPhlAn database dir found at: '{}'!\n".format(mpa_config["bt2_db_dir"])
        err_message += "Specify relevant paths in the metaphlan section of config.yaml.\n"
        err_message += "If you do not want to run MetaPhlAn for taxonomic profiling, set metaphlan: False in config.yaml"
        raise WorkflowError(err_message)

    heatmap = f"{OUTDIR}/metaphlan/all_samples.{mpa_config['heatmap']['level']}_top{mpa_config['heatmap']['topN']}.pdf",
    mpa_area_plot = f"{OUTDIR}/metaphlan/area_plot.metaphlan.pdf",
    krona_plots = expand(f"{OUTDIR}/metaphlan/{{krona}}.metaphlan.krona.html",
        krona=("all_samples","combined_samples"))
    mpa_outputs = expand(f"{OUTDIR}/metaphlan/levels/{{taxlvl}}.tsv",
        taxlvl=("species", "genus", "family", "order"))

    all_outputs.append(heatmap)
    all_outputs.append(mpa_area_plot)
    all_outputs.append(mpa_outputs)

    citations.add(publications["MetaPhlAn"])

    if mpa_config["run_krona"]:
        all_outputs.append(krona_plots)
        citations.add(publications["Krona"])


rule metaphlan:
    """Taxonomic profiling using MetaPhlAn."""
    input:
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
    output:
        bt2_out=f"{OUTDIR}/metaphlan/{{sample}}.bowtie2.bz2" if mpa_config["keep_bt2"] else temp(f"{OUTDIR}/metaphlan/{{sample}}.bowtie2.bz2"),
        mpa_out=f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.txt",
        sam_out=f"{OUTDIR}/metaphlan/{{sample}}.sam.bz2" if mpa_config["keep_sam"] else temp(f"{OUTDIR}/metaphlan/{{sample}}.sam.bz2"),
    log:
        stdout=f"{LOGDIR}/metaphlan/{{sample}}.metaphlan.stdout.log",
        stderr=f"{LOGDIR}/metaphlan/{{sample}}.metaphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    container:
        "docker://quay.io/biocontainers/metaphlan:4.0.2--pyhca03a8a_0"
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
            --nproc 10 \
            --sample_id {wildcards.sample} \
            --samout {output.sam_out} \
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
    """Convert metaphlan output to Krona input"""
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
        set +o pipefail  # Small samples can produce empty output files failing the pipeline
        sed '/#/d' {input.mpa_out} \
            | grep -E "s__|unclassified" \
            | cut -f1,3 \
            | awk '{{print $2,"\t",$1}}' \
            | sed 's/|\w__/\t/g' \
            | sed 's/k__//' \
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
    container:
        "docker://quay.io/biocontainers/metaphlan:4.0.2--pyhca03a8a_0"
    threads:
        1
    shell:
        """
        merge_metaphlan_tables.py {input} > {output.txt} 2> {log}
        sed --in-place 's/\.metaphlan//g' {output.txt} 
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
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    shell:
        """
        scripts/area_plot.py \
            --table {input} \
            --output {output} \
            --mode metaphlan4 \
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
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
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
        html_samples=report(f"{OUTDIR}/metaphlan/all_samples.metaphlan.krona.html",
            category="Taxonomic profiling",
            caption="../../report/metaphlan_krona.rst"),
        html_all=report(f"{OUTDIR}/metaphlan/combined_samples.metaphlan.krona.html",
            category="Taxonomic profiling",
            caption="../../report/metaphlan_krona.rst"),
    log:
        f"{LOGDIR}/metaphlan/create_metaphlan_krona_plots.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    threads:
        1
    shell:
        """
        ktImportText \
            -o {output.html_samples} \
            {input} \
            > {log}

        ktImportText \
            -o {output.html_all} \
            -c \
            {input} \
            >> {log}
        """

rule metaphlan_outputs:
    """Separate the metaphlan abundance table into species, genus, family and order levels"""
    input:
        mpa_combined=f"{OUTDIR}/metaphlan/all_samples.metaphlan.txt",
    output:
        species=f"{OUTDIR}/metaphlan/levels/species.tsv",
        genus=f"{OUTDIR}/metaphlan/levels/genus.tsv",
        family=f"{OUTDIR}/metaphlan/levels/family.tsv",
        order=f"{OUTDIR}/metaphlan/levels/order.tsv",
    shell:
        """
        set +o pipefail
        sed '/#.*/d' {input.mpa_combined} | cut -f 1- | head -n1 | tee {output.species} {output.genus} {output.family} {output.order} > /dev/null

        sed '/#.*/d' {input.mpa_combined} | cut -f 1- | grep s__ | sed 's/^.*s__/s__/g' >> {output.species}
        sed '/#.*/d' {input.mpa_combined} | cut -f 1- | grep g__ | sed 's/^.*s__.*//g' | grep g__ | sed 's/^.*g__/g__/g' >> {output.genus}
        sed '/#.*/d' {input.mpa_combined} | cut -f 1- | grep f__ | sed 's/^.*g__.*//g' | grep f__ | sed 's/^.*f__/f__/g' >> {output.family}
        sed '/#.*/d' {input.mpa_combined} | cut -f 1- | grep o__ | sed 's/^.*f__.*//g' | grep o__ | sed 's/^.*o__/o__/g' >> {output.order}
        """

