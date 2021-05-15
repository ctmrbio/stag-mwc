# vim: syntax=python expandtab
# Metagenomic strain-level population genomics with StrainPhlAn.
 
from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:

mpa_config = config["metaphlan"]
spa_config = config["strainphlan"]
if config["strain_level_profiling"]["strainphlan"]:
    if not mpa_config["bt2_db_dir"] or not Path(mpa_config["bt2_db_dir"]).exists():
        err_message = "No MetaPhlAn database dir found at: '{}'!\n".format(mpa_config["bt2_db_dir"])
        err_message += "bt2_db_dir and bt2_index are required to run StrainPhlAn as it uses the output from MetaPhlAn as input.\n"
        err_message += "Specify relevant paths in the metaphlan section of config.yaml.\n"
        err_message += "If you do not want to run MetaPhlAn or StrainPhlAn, set \"metaphlan: False\" and \"strainphlan: false\" in config.yaml"
        raise WorkflowError(err_message)
    if not spa_config["clade_of_interest"]:
        available_clades=f"{LOGDIR}/strainphlan/available_clades.txt",
        all_outputs.append(available_clades)
        print("Clade of interest not specified in strainphlan section of config.yaml.")
        print("Based on your samples strainphlan will create a list of available clades in output/strainphlan/available_clades.txt")
        print("If you still want to run strainphlan, please update config.yaml e.g. \"clade_of_interest: s__Bifidobacterium_longum\".")
    if spa_config["clade_of_interest"]:
       spa_alignment=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.StrainPhlAn3_concatenated.aln",
       spa_tree=f"{OUTDIR}/strainphlan/RAxML_bestTree.{spa_config['clade_of_interest']}.StrainPhlAn3.tre",
       all_outputs.append(spa_alignment)
       all_outputs.append(spa_tree)
       print("If strainphlan crashes, take a look through output/strainphlan/available_clades.txt to ensure the clade_of_interest you specified is available.")
       citations.add(publications["MetaPhlAn"])
       citations.add(publications["StrainPhlAn"])

rule consensus_markers:
    """Generate consensus markers"""
    input:
        sam=f"{OUTDIR}/metaphlan/{{sample}}.sam.bz2", 
    output:
        consensus_markers=f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl",
    log:
        stdout=f"{LOGDIR}/strainphlan/sample2markers.{{sample}}.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/sample2markers.{{sample}}.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads:
        cluster_config["strainphlan"]["n"] if "strainphlan" in cluster_config else 8
    params:
        output_dir=f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}"
    shell:
        """
        sample2markers.py \
             -i {input.sam} \
             -o {params.output_dir} \
             -n 8 \
             > {log.stdout} \
             2> {log.stderr}
        """

rule print_clades:
    """print a list of available clades"""
    input:
        consensus_markers=expand(f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl", sample=SAMPLES),
    output:
       symlink_to_available_clades=f"{OUTDIR}/strainphlan/available_clades.txt",
    log:
        available_clades=f"{LOGDIR}/strainphlan/available_clades.txt",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads:
        cluster_config["strainphlan"]["n"] if "strainphlan" in cluster_config else 8
    params:
        out_dir=f"{OUTDIR}/strainphlan",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
    shell:
        """
        strainphlan \
             -s {input.consensus_markers} \
             --print_clades_only \
             -d {params.database} \
             -o {params.out_dir} \
             -n {threads} \
             2>&1 > {log.available_clades}

        cd {params.out_dir} && ln -s ../logs/strainphlan/available_clades.txt
        """

rule extract_markers:
    """extract marker sequences for clade of interest"""
    input:
        available_clades=f"{OUTDIR}/strainphlan/available_clades.txt",
        consensus_markers=expand(f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl", sample=SAMPLES),
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
        "shub://AroArz/singularity_playground:biobakery"
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
        consensus_markers=expand(f"{OUTDIR}/strainphlan/consensus_markers/{{sample}}/{{sample}}.pkl", sample=SAMPLES),
        reference_markers=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.fna",
    output:
        alignment=f"{OUTDIR}/strainphlan/{spa_config['clade_of_interest']}.StrainPhlAn3_concatenated.aln",
        tree=f"{OUTDIR}/strainphlan/RAxML_bestTree.{spa_config['clade_of_interest']}.StrainPhlAn3.tre",
    log:
        stdout=f"{LOGDIR}/strainphlan/alignment.strainphlan.stdout.log",
        stderr=f"{LOGDIR}/strainphlan/alignment.strainphlan.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads:
        cluster_config["strainphlan"]["n"] if "strainphlan" in cluster_config else 8
    params:
        clade=spa_config["clade_of_interest"],
        out_dir=f"{OUTDIR}/strainphlan",
        database=f"{mpa_config['bt2_db_dir']}/{mpa_config['bt2_index']}.pkl",
        mode=spa_config["mode"],
        extra=spa_config["extra"],  # This is extremely useful if you want to include a reference genome
    shell:
        """
        echo "please compare your clade_of_interest to list of available clades in available_clades.txt" > {log.stderr}

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
             2>> {log.stderr}
        """

