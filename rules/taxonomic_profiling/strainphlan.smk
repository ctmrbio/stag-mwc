# vim: syntax=python expandtab
# Metagenomic strain-level population genomics with StrainPhlAn.
 
from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
#     extract_markers,

mpa_config = config["metaphlan"]
spa_config = config["strainphlan"]
if config["strain_level_profiling"]["strainphlan"]:
    if not mpa_config["bt2_db_dir"] or not Path(mpa_config["bt2_db_dir"]).exists():
        err_message = "No MetaPhlAn database dir found at: '{}'!\n".format(mpa_config["bt2_db_dir"])
        err_message += "bt2_db_dir and bt2_index are required to run StrainPhlAn as it uses the output from MetaPhlAn as input.\n"
        err_message += "Specify relevant paths in the metaphlan section of config.yaml.\n"
        err_message += "If you do not want to run MetaPhlAn or StrainPhlAn, set \"metaphlan: False\" and \"strainphlan: false\" in config.yaml"
        raise WorkflowError(err_message)

    spa_outputs=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.tre"
    all_outputs.append(spa_outputs)

    citations.add(publications["MetaPhlAn"])
    citations.add(publications["Krona"])

rule consensus_markers:
    """Generate consensus markers"""
    input:
        sam=f"{OUTDIR}/metaphlan/{{sample}}.sam.bz2", 
    output:
        consensus_markers=f"{OUTDIR}/strainphlan/{{sample}}.pkl",
    log:
        stdout=f"{LOGDIR}/strainphlan/sample2markers.{{sample}}.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/sample2markers.{{sample}}.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc-biobakery"
    threads:
        cluster_config["strainphlan"]["n"] if "strainphlan" in cluster_config else 8
    params:
        output_dir=f"{OUTDIR}/strainphlan/"
    shell:
        """
        sample2markers.py \
             -i {input.sam} \
             -o {params.output_dir} \
             -n 8 \
             > {log.stdout} \
             2> {log.stderr}
        """

rule extract_markers:
    """extract marker sequences"""
    input:
        consensus_markers=expand(f"{OUTDIR}/strainphlan/{{sample}}.pkl", sample=SAMPLES),
    output:
        reference_markers=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.fna",
    log:
        stdout=f"{LOGDIR}/strainphlan/extract_markers.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/extract_markers.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc-biobakery"
    threads:
        cluster_config["strainphlan"]["n"] if "strainphlan" in cluster_config else 8
    params:
        clade=spa_config["clade_of_interest"],
        out_dir=f"{OUTDIR}/strainphlan/",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
    shell:
        """
        extract_markers.py \
             -c {params.clade} \
             -o {params.out_dir} \
             -d {params.database} \
             > {log.stdout} \
             2> {log.stderr}
        """

rule strainphlan:
    """generate tree and alignment"""
    input:
        consensus_markers=expand(f"{OUTDIR}/strainphlan/{{sample}}.pkl", sample=SAMPLES),
        reference_markers=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.fna",
    output:
        alignment=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}_concatenated.aln",
        tree=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.tre",
    log:
        stdout=f"{LOGDIR}/strainphlan/alignment.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/alignment.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc-biobakery"
    threads:
        cluster_config["strainphlan"]["n"] if "strainphlan" in cluster_config else 8
    params:
        clade=spa_config["clade_of_interest"],
        out_dir=f"{OUTDIR}/strainphlan/",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
        mode=spa_config["mode"],
        extra=spa_config["extra"],  # This is extremely useful if you want to include a reference genome
    shell:
        """
        strainphlan \
             -s {input.consensus_markers} \
             -m {input.reference_markers} \
             {params.extra} \
             -d {params.database} \
             -o {params.out_dir} \
             -n {threads} \
             -c {params.clade} \
             --phylophlan_mode accurate \
             --mutation_rates \
             > {log.stdout} \
             2> {log.stderr}
        """
# 
# rule visualize_tree:
#     """visualize"""
#     input:
#         consensus_markers=expand(f"{OUTDIR}/strainphlan/{{sample}}.pkl", sample=SAMPLES),
#      reference_markers=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.fna",
#     output:
#         alignment=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}_concatenated.aln",
#      tree=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.tre"
#     log:
#         stdout=f"{LOGDIR}/strainphlan/alignment.strainphlan.stdout.log",
#         stderr=f"{LOGDIR}/strainphlan/alignment.strainphlan.stderr.log",
#     shadow:
#         "shallow"
#     conda:
#         "../../envs/metaphlan.yaml"
#     singularity:
#         "shub://ctmrbio/stag-mwc:stag-mwc-biobakery"
#     threads:
#         cluster_config["strainphlan"]["n"] if "strainphlan" in cluster_config else 8
#     params:
#         clade=spa_config["clade_of_interest"],
#      out_dir=f"{OUTDIR}/strainphlan/",
#      database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_idex']}.pkl",
#      mode=spa_config["phylophlan_mode"],
#      extra=spa_config["extra"],  # This is extremely useful if you want to include a reference genome
#     shell:
#         """
#      strainphlan \
#           -s {input.consensus_markers} \
#           -m {input.reference_markers} \
#           {params.extra} \
#           -o {output.alignment} \
#           -n {threads} \
#           -c {params.clade} \
#           --phylophlan_mode accurate \
#           --mutation_rates \
#           > {log.stdout} \
#           2> {log.stderr}
#      """


