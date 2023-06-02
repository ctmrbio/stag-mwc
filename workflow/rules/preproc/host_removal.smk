# vim: syntax=python expandtab
# Rules to filter host sequences from metagenomic reads
from pathlib import Path 

from snakemake.exceptions import WorkflowError


if config["host_removal"]["kraken2"] and config["host_removal"]["bowtie2"]:
    err_message = "Running both Kraken2 and Bowtie2 for host removal is not supported"
    raise WorkflowError(err_message)

#########################################
#              kraken2
#########################################
rh_kraken2 = config["remove_host"]["kraken2"]
if config["host_removal"]["kraken2"]:
    k2_db_path = Path(rh_kraken2["db_path"])
    if not Path(k2_db_path/"taxo.k2d").is_file():
        err_message = "Cannot find Kraken2 database for host sequence removal at: '{}/*.k2d'!\n".format(k2_db_path)
        err_message += "Specify path to folder containing Kraken2 database for host removal in config.yaml.\n"
        raise WorkflowError(err_message)

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    k2_filtered = expand(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz",
            sample=SAMPLES,
            readpair=[1,2])
    host_proportions = OUTDIR/"host_removal/host_proportions.txt"
    if rh_kraken2["keep_fastq"]:
        all_outputs.extend(k2_filtered)
    all_outputs.append(host_proportions)

    citations.add(publications["Kraken2"])

    localrules:
        plot_proportion_host_kraken2

    rule kraken2_host_removal:
        """Filter reads matching host database using Kraken2."""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz" if rh_kraken2["keep_fastq"] else temp(OUTDIR/"host_removal/{sample}_1.fq.gz"),
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz" if rh_kraken2["keep_fastq"] else temp(OUTDIR/"host_removal/{sample}_2.fq.gz"),
            host1=OUTDIR/"host_removal/{sample}.host_1.fq.gz" if rh_kraken2["keep_host_fastq"] else temp(OUTDIR/"host_removal/{sample}.host_1.fq.gz"),
            host2=OUTDIR/"host_removal/{sample}.host_2.fq.gz" if rh_kraken2["keep_host_fastq"] else temp(OUTDIR/"host_removal/{sample}.host_2.fq.gz"),
            kraken=OUTDIR/"host_removal/{sample}.kraken" if rh_kraken2["keep_kraken"] else temp(OUTDIR/"host_removal/{sample}.kraken"),
            kreport=OUTDIR/"host_removal/{sample}.kreport" if rh_kraken2["keep_kreport"] else temp(OUTDIR/"host_removal/{sample}.kreport"),
        log:
            stderr=LOGDIR/"host_removal/{sample}.kraken2.log",
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 8
        params:
            db=rh_kraken2["db_path"],
            confidence=rh_kraken2["confidence"],
            extra=rh_kraken2["extra"],
            classified=lambda w: f"{OUTDIR}/host_removal/{w.sample}.host#.fq",
            unclassified=lambda w: f"{OUTDIR}/host_removal/{w.sample}#.fq",
            fq_to_compress=lambda w: f"{OUTDIR}/host_removal/{w.sample}*.fq",
        shell:
            """
            kraken2 \
                --db {params.db} \
                --threads {threads} \
                --output {output.kraken} \
                --classified-out {params.classified} \
                --unclassified-out {params.unclassified} \
                --report  {output.kreport} \
                --paired \
                --confidence {params.confidence} \
                {params.extra} \
                {input.read1} {input.read2} \
                2> {log.stderr}
            pigz \
                --processes {threads} \
                --verbose \
                --force \
                {params.fq_to_compress} \
                2>> {log.stderr}
            """


    rule plot_proportion_host_kraken2:
        """Plot proportion of reads that matched the host DB."""
        input:
            expand(LOGDIR/"host_removal/{sample}.kraken2.log", sample=SAMPLES)
        output:
            histogram=report(OUTDIR/"host_removal/host_histogram.pdf",
                       category="Preprocessing",
                       caption="../../report/host_histogram.rst"),
            barplot=report(OUTDIR/"host_removal/host_barplot.pdf",
                       category="Preprocessing",
                       caption="../../report/host_barplot.rst"),
            txt=report(OUTDIR/"host_removal/host_proportions.txt",
                       category="Preprocessing",
                       caption="../../report/host_proportions.rst"),
        log:
            LOGDIR/"host_removal/proportion_host.log",
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 1
        shell:
            """
            workflow/scripts/plot_proportion_kraken2.py \
                {input} \
                --histogram {output.histogram} \
                --barplot {output.barplot} \
                --table {output.txt} \
                2>&1 > {log}
            """

#########################################
#              bowtie2
#########################################
rh_bowtie2 = config["remove_host"]["bowtie2"]
if config["host_removal"]["bowtie2"]:
    bt2_db_path = Path(rh_bowtie2["db_path"])
    if not list(bt2_db_path.parent.glob(bt2_db_path.name+".1.bt2*")):
        err_message = "Cannot find Bowtie2 database for host sequence removal at: '{}/'!\n".format(bt2_db_path)
        err_message += "Specify path to folder containing Bowtie2 database for host removal in config.yaml.\n"
        err_message += "Note that the path should not contain the '.1.bt2' suffix."
        raise WorkflowError(err_message)

    bt2_filtered = expand(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz",
            sample=SAMPLES,
            readpair=[1,2])
    all_outputs.extend(bt2_filtered)

    citations.add(publications["Bowtie2"])


    rule bowtie2_host_removal:
        """Map reads against host sequence database."""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            sam=temp(OUTDIR/"host_removal/{sample}.sam"),
        log:
            stderr=LOGDIR/"host_removal/{sample}.bowtie2.stderr",
        threads:
            12
        conda:
            "../../envs/metaphlan.yaml"
        container:
            "docker://quay.io/biocontainers/metaphlan:4.0.3--pyhca03a8a_0"
        params:
            db_path=rh_bowtie2["db_path"],
            extra=rh_bowtie2["extra"],
        shell:
            """
            bowtie2 \
                --threads {threads} \
                -x {params.db_path} \
                -1 {input.read1} \
                -2 {input.read2} \
                -S {output.sam} \
                2> {log.stderr}
            """

    rule bt2_sam2bam:
        """Convert SAM to BAM."""
        input:
            sam=rules.bowtie2_host_removal.output.sam,
        output:
            bam=temp(OUTDIR/"host_removal/{sample}.bam"),
        log:
            stderr=LOGDIR/"host_removal/{sample}.sam2bam.stderr",
        threads:
            4
        conda:
            "../../envs/metaphlan.yaml"
        container:
            "docker://quay.io/biocontainers/metaphlan:4.0.3--pyhca03a8a_0"
        shell:
            """
            samtools view \
                -b \
                --threads {threads} \
                {input.sam} \
                -o {output.bam} \
                2> {log.stderr}
            """

    rule bt2_get_unmapped_pairs:
        """SAM-flag filter: get unmapped pairs."""
        input:
            bam2=rules.bt2_sam2bam.output.bam,
        output:
            unmapped=temp(OUTDIR/"host_removal/{sample}_unmapped.bam"),
        log:
            stderr=LOGDIR/"host_removal/{sample}.unmapped.stderr",
        threads:
            4
        conda:
            "../../envs/metaphlan.yaml"
        container:
            "docker://quay.io/biocontainers/metaphlan:4.0.3--pyhca03a8a_0"
        shell:
            """
            samtools view \
                -b \
                -f 13 \
                -F 256 \
                --threads {threads} \
                {input.bam2} \
                -o {output.unmapped} \
                2> {log.stderr}
            """

    rule bt2_sort_bam_files:
        """Sort BAM file by read name ( -n ) to get paired reads."""
        input:
            pairs=rules.bt2_get_unmapped_pairs.output.unmapped,
        output:
            sorted=temp(OUTDIR/"host_removal/{sample}_unmapped.sorted.bam"),
        log:
            stderr=LOGDIR/"host_removal/{sample}.sort.stderr",
        threads:
            4
        conda:
            "../../envs/metaphlan.yaml"
        container:
            "docker://quay.io/biocontainers/metaphlan:4.0.3--pyhca03a8a_0"
        shell: 
            """
            samtools sort \
                -n \
                -m 5G \
                --threads {threads} \
                {input.pairs} \
                -o {output.sorted} \
                2> {log.stderr}
            """

    rule bt2_get_read_pairs:
        """Output paired reads in FASTQ format."""
        input:
            sorted_pairs=rules.bt2_sort_bam_files.output.sorted,
        output:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz" if rh_bowtie2["keep_fastq"] else temp(OUTDIR/"host_removal/{sample}_1.fq.gz"),
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz" if rh_bowtie2["keep_fastq"] else temp(OUTDIR/"host_removal/{sample}_2.fq.gz"),
        log:
            stderr=LOGDIR/"host_removal/{sample}.samtools.fastq.log",
        threads:
            4
        conda:
            "../../envs/metaphlan.yaml"
        container:
            "docker://quay.io/biocontainers/metaphlan:4.0.3--pyhca03a8a_0"
        shell: 
            """
            samtools fastq \
                --threads {threads} \
                -1 {output.read1} \
                -2 {output.read2} \
                -0 /dev/null \
                -s /dev/null \
                -n \
                {input.sorted_pairs} \
                2> {log.stderr}
            """

#########################################
#           skip host removal
#########################################
if not any(config["host_removal"].values()):
    if not config["fastp"]["keep_output"]:
        err_message = "Set fastp keep_output in config.yaml to True in order to skip host removal.\n"
        err_message += "If you want to run host removal set remove_host in config.yaml to True"
        raise WorkflowError(err_message)
    
    filtered_host = expand(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz",
            sample=SAMPLES,
            readpair=[1,2])
    all_outputs.extend(filtered_host)

    localrules:
        skip_remove_host,

    rule skip_remove_host:
        """Do not remove host sequences"""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
        log:
            stderr=str(LOGDIR/"host_removal/{sample}.stderr.log"),
        shell:
            """
            ln -sv $(readlink -f {input.read1}) {output.read1} >> {log.stderr}
            ln -sv $(readlink -f {input.read2}) {output.read2} >> {log.stderr}
            """
