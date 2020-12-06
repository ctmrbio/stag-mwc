# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using MetaPhlAn

from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
    plot_metaphlan_heatmap,
    create_metaphlan_krona_plots,
    metaphlan_krona,
    gather_visualizations,

mpa_config = config["metaphlan"]
if config["taxonomic_profile"]["metaphlan"]:
    if not mpa_config["bt2_db_dir"] or not Path(mpa_config["bt2_db_dir"]).exists():
        err_message = "No MetaPhlAn database dir found at: '{}'!\n".format(mpa_config["bt2_db_dir"])
        err_message += "Specify relevant paths in the metaphlan section of config.yaml.\n"
        err_message += "If you do not want to run MetaPhlAn for taxonomic profiling, set metaphlan: False in config.yaml"
        raise WorkflowError(err_message)

    mpa_outputs=f"{OUTDIR}/metaphlan/visualizations_success.txt"
    all_outputs.append(mpa_outputs)

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
        sam_out=f"{OUTDIR}/metaphlan/{{sample}}.sam.bz2",
    log:
        stdout=f"{LOGDIR}/metaphlan/{{sample}}.metaphlan.stdout.log",
        stderr=f"{LOGDIR}/metaphlan/{{sample}}.metaphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
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
        "shub://AroArz/singularity_playground:biobakery"
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
        html=report(f"{OUTDIR}/metaphlan/all_samples.metaphlan.krona.html",
                    category="Taxonomic profiling",
                    caption="../../report/metaphlan_krona.rst"),
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads:
        1
    shell:
        """
        ktImportText \
            -o {output.html} \
            {input}
        """

rule metaphlan_hclust2:
    """Create hclust2 species heatmap"""
    input:
        f"{OUTDIR}/metaphlan/all_samples.metaphlan.txt",
    output:
        hclust_heatmap=f"{OUTDIR}/metaphlan/abundance_heatmap_species.png",
        species=f"{OUTDIR}/metaphlan/merged_abundance_table_species.txt",
    log:
        f"{LOGDIR}/metaphlan/merged_abundance_table_species.log",
    shadow:
        "shallow"
    conda:
        "../../envs/graphlan.yaml"
    singularity:
        "shub://AroArz/singularity_playground:graphlan"
    params:
        scale=mpa_config["hclust_heatmap"]["scale"],
        feature_distance=mpa_config["hclust_heatmap"]["feature_distance"],
        sample_distance=mpa_config["hclust_heatmap"]["sample_distance"],
    shell:
        """
        grep -E "s__|clade" {input} \
            | sed 's/^.*s__//g' \
            | cut -f1,3-8 \
            | sed -e 's/clade_name/species/g' \
            > {output.species}

        hclust2.py \
            -i {output.species} \
            -o {output.hclust_heatmap} \
            --ftop 25 \
            --f_dist_f {params.feature_distance} \
            --s_dist_f {params.sample_distance} \
            --cell_aspect_ratio 0.5 \
            {params.scale} \
            --flabel_size 6 \
            --slabel_size 6 \
            --max_flabel_len 100 \
            --max_slabel_len 100 \
            --minv 0.1 \
            --dpi 300 \
            2> {log}
        """

rule metaphlan_cladogram:
    input:
        txt=f"{OUTDIR}/metaphlan/all_samples.metaphlan.txt",
        hclust_heatmap=f"{OUTDIR}/metaphlan/abundance_heatmap_species.png",
    output:
        clado_format=f"{OUTDIR}/metaphlan/all_samples.metaphlan.clado.txt",
        abundance_tree=f"{OUTDIR}/metaphlan/all_samples_abundance.tree.txt",
        abundance_annot=f"{OUTDIR}/metaphlan/all_samples_abundance.annot.txt",
        abundance_xml=f"{OUTDIR}/metaphlan/all_samples_abundance.xml",
        cladogram=f"{OUTDIR}/metaphlan/all_samples_cladogram.png",
    log:
        f"{LOGDIR}/metaphlan/merged_abundance_table_species.log",
    shadow:
        "shallow"
    conda:
        "../../envs/graphlan.yaml"
    singularity:
        "shub://AroArz/singularity_playground:graphlan"
    params:
        skip_rows=mpa_config["cladogram"]["skip_rows"],
        most_abundance=mpa_config["cladogram"]["most_abundance"],
        threshold=mpa_config["cladogram"]["minimum_abundance"],
        least_biomarkers=mpa_config["cladogram"]["least_biomarkers"],
        annotation=mpa_config["cladogram"]["annotation"],
        external_annotation=mpa_config["cladogram"]["external_annotation"],
        min_clade_size=mpa_config["cladogram"]["min_clade_size"],
    shell:
        """
        tail \
            -n +2 {input.txt} \
            | cut -f1,3- \
            > {output.clado_format}
        
        export2graphlan.py \
            --skip_rows {params.skip_rows} \
            -i {output.clado_format} \
            --tree {output.abundance_tree} \
            --annotation {output.abundance_annot} \
            --most_abundant {params.most_abundance} \
            --abundance_threshold {params.threshold} \
            --least_biomarkers {params.least_biomarkers} \
            --annotations {params.annotation} \
            --external_annotations {params.external_annotation} \
            --min_clade_size {params.min_clade_size} \
            2> {log}

        graphlan_annotate.py \
            --annot \
                {output.abundance_annot} \
                {output.abundance_tree} \
                {output.abundance_xml} \
                2>> {log}

        graphlan.py \
            --dpi 300 \
                {output.abundance_xml} \
                {output.cladogram} \
            --external_legends \
            2>> {log}
        """



rule gather_visualizations:
    input:
        mpa_combined = expand(f"{OUTDIR}/metaphlan/all_samples.metaphlan.{{ext}}", ext=("txt", "krona.html")),
        hclust_heatmap=f"{OUTDIR}/metaphlan/abundance_heatmap_species.png",
        species=f"{OUTDIR}/metaphlan/merged_abundance_table_species.txt",
        mpa_plot = f"{OUTDIR}/metaphlan/all_samples.{mpa_config['heatmap']['level']}_top{mpa_config['heatmap']['topN']}.pdf",
        mpa_area_plot = f"{OUTDIR}/metaphlan/area_plot.metaphlan.pdf",
        cladogram=f"{OUTDIR}/metaphlan/all_samples_cladogram.png",
    output:
        visualuzations_success=f"{OUTDIR}/metaphlan/visualizations_success.txt"
    shell:
        """
        echo "area_plot done" > {output}
        echo "heatmap done" >> {output}
        echo "kron_plot done" >> {output}
        echo "hclust heatmap done" >> {output}
        echo "cladogram done" >> {output}
        """
