# vim: syntax=python expandtab
# Compare all samples against all samples using MinHash sketches
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.

localrules:
    compare_sketches,
    plot_sample_similarity,


if config["naive"]["sketch_compare"]:
    # Add final output files from this module to 'all_outputs' from the
    # main Snakefile scope.
    sample_similarity_plot = str(OUTDIR/"sketch_compare/sample_similarity.pdf")
    all_outputs.append(sample_similarity_plot)

    citations.add(publications["BBMap"])


rule sketch:
    """Create MinHash sketches of samples using BBMap's sketch.sh.
    Uses only the first readpair of each sample."""
    input:
        INPUT_read1,
    output:
        sketch=OUTDIR/"sketch_compare/{sample}.sketch.gz",
    log:
        str(LOGDIR/"sketch_compare/{sample}.sketch.log")
    shadow: 
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    threads: 4
    shell:
        """
        sketch.sh \
            in={input} \
            out={output} \
            name0={wildcards.sample} \
            2> {log}
        """


rule compare_sketches:
    """Compare all samples using BBMap's comparesketch.sh"""
    input:
        samples=expand(str(OUTDIR/"sketch_compare/{sample}.sketch.gz"), 
                        sample=SAMPLES)
    output:
        alltoall=OUTDIR/"sketch_compare/alltoall.txt",
    log:
        str(LOGDIR/"sketch_compare/comparesketch.log")
    shadow: 
        "shallow"
    conda: 
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    shell:
        """
        comparesketch.sh \
            format=3 \
            out={output} \
            alltoall \
            {input} \
            2> {log}
        """


rule plot_sample_similarity:
    """Plot sample sketch similarity matrix"""
    input:
        OUTDIR/"sketch_compare/alltoall.txt"
    output:
        heatmap=OUTDIR/"sketch_compare/sample_similarity.pdf",
        clustered=report(OUTDIR/"sketch_compare/sample_similarity.clustered.pdf", 
                         category="Sketch comparison", 
                         caption="../../report/sketch_compare.rst")
    log:
        stdout=str(LOGDIR/"sketch_compare/sample_similarity_plot.stdout.log"),
        stderr=str(LOGDIR/"sketch_compare/sample_similarity_plot.stderr.log"),
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    shell:
        """
        workflow/scripts/plot_sketch_comparison_heatmap.py \
            --outfile {output.heatmap} \
            --clustered {output.clustered} \
            {input} \
            > {log.stdout} \
            2> {log.stderr}
        """
