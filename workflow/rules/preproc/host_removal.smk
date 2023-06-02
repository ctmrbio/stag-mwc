# vim: syntax=python expandtab
# Rules to filter host sequences from metagenomic reads
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
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
    db_path = Path(rh_kraken2["db_path"])
    if not Path(db_path/"taxo.k2d").is_file():
        err_message = "Cannot find Kraken2 database for host sequence removal at: '{}/*.k2d'!\n".format(db_path)
        err_message += "Specify path to folder containing Kraken2 database for host removal in config.yaml.\n"
        raise WorkflowError(err_message)

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    filtered_host = expand(str(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1,2])
    host_proportions = str(OUTDIR/"host_removal/host_proportions.txt")
    if rh_kraken2["keep_fastq"]:
        all_outputs.extend(filtered_host)
    all_outputs.append(host_proportions)

    citations.add(publications["Kraken2"])

    localrules:
        plot_proportion_host

    rule remove_host:
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
            stderr=str(LOGDIR/"host_removal/{sample}.kraken2.log"),
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


    rule plot_proportion_host:
        """Plot proportion of reads that matched the host DB."""
        input:
            expand(str(LOGDIR/"host_removal/{sample}.kraken2.log"), sample=SAMPLES)
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
            str(LOGDIR/"host_removal/proportion_host.log")
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
    rule bowtie2_host_removal:
        """Map reads against host sequence database."""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            sam=temp(OUTDIR/"host_removal/{sample}.sam"),
        log:
            stdout=OUTDIR/"logs/host_removal/{sample}.bowtie2.stdout",
            stderr=OUTDIR/"logs/host_removal/{sample}.bowtie2.stderr",
        threads:
            10
        conda:
            "envs/conda.yaml"
        params:
            db_path="/ceph/db/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as",
        shell:
            """
            bowtie2 -p 8 -x {params.db_path} \
                -1 {input.read1} \
                -2 {input.read2} \
                -S {output.file} \
                > {log.stdout} \
                2> {log.stderr}
            """

    rule bt2_sam2bam:
        """Convert SAM to BAM."""
        input:
            sam=rules.bowtie2_host_removal.output.sam,
        output:
            bam=temp("output/host_removal/{sample}.bam"),
        threads:
            10
        conda:
            "envs/conda.yaml"
        shell:
            """
            samtools view \
                -bS {input.sam} \
                -o {output.bam}
            """

    rule bt2_get_unmapped_pairs:
        """SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)."""
        input:
            bam2=rules.bt2_sam2bam.output.bam,
        output:
            unmapped=temp("output/host_removal/{sample}_unmapped.bam"),
        threads:
            10
        conda:
            "envs/conda.yaml"
        shell:
            """
            samtools view \
                -b \
                -f 12 \
                -F 256 {input.bam2} \
                -o {output.unmapped}
            """

    rule bt2_sort_bam_files:
        """Sort BAM file by read name ( -n ) to get paired reads."""
        input:
            pairs=rules.bt2_get_unmapped_pairs.output.unmapped,
        output:
            sorted=temp("output/sam_files/{sample}_unmapped.sorted.bam"),
        threads:
            10
        conda:
            "envs/conda.yaml"
        shell: 
            """
            samtools sort \
                -n \
                -m 5G \
                -@ 2 {input.pairs} \
                -o {output.sorted}
            """

    rule bt2_get_read_pairs:
        """Output paired reads in FASTQ format."""
        input:
            sorted_pairs=rules.bt2_sort_bam_files.output.sorted,
        output:
            pair1="output/{sample}_host_removed_1.fq.gz",
            pair2="output/{sample}_host_removed_2.fq.gz",
        log:
            stderr="output/logs/{sample}.log",
        threads:
            10
        conda:
            "envs/conda.yaml"
        shell: 
            """
            samtools fastq \
                -@ 8 {input.sorted_pairs} \
                -1 {output.pair1} \
                -2 {output.pair2} \
                -0 /dev/null -s /dev/null -n \
                2> {log.stderr}
            """

#########################################
#           skip host removal
#########################################
if not config["host_removal"]["kraken2"] or config["host_removal"]["bowtie2"]:
    if not config["fastp"]["keep_output"]:
        err_message = "Set fastp keep_output in config.yaml to True in order to skip host removal.\n"
        err_message += "If you want to run host removal set remove_host in config.yaml to True"
        raise WorkflowError(err_message)
    
    filtered_host = expand(str(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz"),
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
