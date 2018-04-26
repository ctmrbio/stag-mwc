# vim: syntax=python expandtab
# Compare all samples against all samples using MinHash sketches
import os.path

localrules:
    compare_sketches
    plot_sample_similarity


# Add final output files from this module to 'all_outputs' from the
# main Snakefile scope.
sample_similarity_plot = expand("{outdir}/sketch_compare/sample_similarity.pdf",
        outdir=outdir)
all_outputs.extend(sample_similarity_plot)


rule sketch:
    """Create MinHash sketches of samples using BBMap's sketch.sh.
    Uses only the first readpair of each sample."""
    input:
        os.path.join(config["inputdir"], config["input_fn_pattern"]).format(sample="{sample}", readpair="1")
    output:
        sketch=config["outdir"]+"/sketch_compare/{sample}.sketch.gz",
    log:
        config["outdir"]+"/logs/sketch_compare/{sample}.sketch.log"
    shadow: 
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        2
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
        samples=expand(config["outdir"]+"/sketch_compare/{sample}.sketch.gz",
                sample=SAMPLES)
    output:
        alltoall=config["outdir"]+"/sketch_compare/alltoall.txt",
    log:
        config["outdir"]+"/logs/sketch_compare/comparesketch.log"
    shadow: 
        "shallow"
    conda: 
        "../../envs/stag-mwc.yaml"
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
        config["outdir"]+"/sketch_compare/alltoall.txt"
    output:
        config["outdir"]+"/sketch_compare/sample_similarity.pdf"
    log:
        stdout=config["outdir"]+"/logs/sketch_compare/sample_similarity_plot.stdout.log",
        stderr=config["outdir"]+"/logs/sketch_compare/sample_similarity_plot.stderr.log",
    conda:
        "../../envs/stag-mwc.yaml"
    shell:
        """
        scripts/plot_sketch_comparison_heatmap.py \
            --outfile {output} \
            {input} \
            > {log.stdout} \
            2> {log.stderr}
        """
