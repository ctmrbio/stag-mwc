.. _BBCountUnique: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/calcuniqueness-guide/
.. _FastP:  https://github.com/OpenGene/fastp
.. _Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _BBMap: https://sourceforge.net/projects/bbmap/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Kaiju: http://kaiju.binf.ku.dk/
.. _Kraken2: https://ccb.jhu.edu/software/kraken2/
.. _KrakenUniq: https://github.com/fbreitwieser/krakenuniq
.. _Bracken: https://ccb.jhu.edu/software/bracken/
.. _MetaPhlAn: https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4
.. _StrainpPlAn: https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4
.. _featureCounts: https://subread.sourceforge.net/featureCounts.html
.. _HUMAnN: https://github.com/biobakery/biobakery/wiki/humann3
.. _GTF format: https://genome.ucsc.edu/FAQ/FAQformat.html#format4
.. _SAF format: https://subread.sourceforge.net/featureCounts.html
.. _MEGAHIT: https://github.com/voutcn/megahit
.. _MultiQC: https://multiqc.info/
.. _megares: https://megares.meglab.org/

Modules
=======
|full_name| is a workflow framework that connects several other tools. The
basic assumption is that all analyses start with a quality control of the
sequencing reads (using `FastP`_), followed by host sequence removal (using
`Kraken2`_ or `Bowtie2`_). This section of the documentation aims to describe
useful details about the separate tools that are used in |full_name|.

The following subsections describe the function of each module included in
|full_name|.  Each tool produces output in a separate subfolder inside the
output folder specified in the configuration file.  The output folders
mentioned underneath each module heading refers to the subfolder inside the
main output folder specified in the configuration file.

Pre-processing
**************

qc_reads
--------------
:Tools: `FastP`_
:Output folder: ``fastp``

The quality control module uses `FastP`_ to produce HTML reports of the quality
of the input and output reads. The quality control module also trims adapter
sequences and performs quality trimming of the input reads. The quality assured
reads are output into ``fastp``. Output filenames are::

    <sample>_{1,2}.fq.gz

.. note:: 

    It is possible to skip fastp processing. StaG then replaces the
    output files with symlinks to the input files.

remove_host
--------------
The ``remove_host`` module can use either `Kraken2`_ or `Bowtie2`_ to classify
reads against a database of host sequences to remove reads matching to
non-desired host genomes.

.. note::

    It is possible to skip host removal. StaG then replaces the output files
    with symlinks to the fastp output files.


:Tool: `Kraken2`_
:Output folder: ``host_removal``

The output from Kraken2 are two sets of pairs of paired-end FASTQ files, and
optionally one Kraken2 classification file and one Kraken2 summary report.  In
addition, two PDF files with 1) a basic histogram plot of the proportion of
host reads detected in each sample, and 2) a barplot of the same. A TSV table
with the raw proportion data is also provided::

    <sample>_{1,2}.fq.gz
    <sample>.host_{1,2}.fq.gz
    host_barplot.pdf
    host_histogram.pdf
    host_proportions.txt


:Tool: `Bowtie2`_
:Output folder: ``host_removal``

The output from Bowtie2 is a set of paired-end FASTQ files::

    <sample>_{1,2}.fq.gz



preprocessing_summary
---------------------
This module summarize the number of reads passing through each preprocessing
step and produces a summary table showing the number of reads after each step.
For more detailed information about read QC please refer to the MulitQC report.


multiqc
-------
:Tool: `MultiQC`_
:Output folder: ``multiqc``

`MultiQC`_ summarizes information about several steps in StaG in an easy-to-use
HTML report. Refer to this report for details about e.g. read QC.


Naive sample analysis
***********************

sketch_compare
--------------
:Tools: ``sketch.sh``, ``comparesketch.sh`` from `BBMap`_
:Output folder: ``sketch_compare``

The ``sketch_compare`` module uses the very fast MinHash implementation in
`BBMap`_ to compute MinHash sketches of all samples to do an all-vs-all
comparison of all samples based on their kmer content. The module outputs
gzip-compressed sketches for each sample, as well as two heatmap plots showing
the overall similarity of all samples (one with hierarchical clustering).


assess_depth
--------------
:Tool: `BBCountUnique`_
:Output folder: ``bbcountunique``

The ``assess_depth`` module uses `BBMap`_'s very fast kmer counting algorithms
to produce saturation curves. The saturation curve shows a histogram of the
proportion of unique kmers observed per reads processed, and can be used to
assess how deep a sample has been sequenced. The module outputs one plot per
sample.



Taxonomic profiling
*******************

Kaiju
-----
:Tool: `Kaiju`_
:Output folder: ``kaiju``

Run `Kaiju`_ on the trimmed and filtered reads to produce a taxonomic profile.
Outputs several files per sample (one per taxonomic level specified in the
config), plus two files that combine all samples in the run: an HTML Krona
report with the profiles of all samples and a TSV table per taxonomic level.
The output files are::

    <sample>.kaiju
    <sample>.kaiju.<level>.txt
    <sample>.krona
    all_samples.kaiju.krona.html
    all_samples.kaiju.<level>.txt

Kraken2
-------
:Tool: `Kraken2`_
:Output folder: ``kraken2``

Run `Kraken2`_ on the trimmed and filtered reads to produce a taxonomic profile. 
Optionally Bracken can be run to produce abundance profiles for each sample at a
user-specified taxonomic level. The Kraken2 module produces the following files::

    <sample>.kraken
    <sample>.kreport
    all_samples.kraken2.txt
    all_samples.mpa_style.txt

This modules outputs two tables containing the same information in two formats:
one is the default Kraken2 output format, the other is a MetaPhlAn2-like format
(``mpa_style``). The optional Bracken further adds additional output files for
each sample::

    <sample>.<taxonomic_level>.bracken
    <sample>.<taxonomic_level>.filtered.bracken
    <sample>_bracken.kreport
    <sample>.bracken.mpa_style.txt
    all_samples.<taxonomic_level>.bracken.txt
    all_samples.<taxonomic_level>.filtered.bracken.txt
    all_samples.bracken.mpa_style.txt
    
KrakenUniq
----------
:Tool: `KrakenUniq`_
:Output folder: ``krakenuniq``

Run `KrakenUniq`_ on the trimmed and filtered reads to produce a taxonomic profile.
The KrakenUniq module produces the following output files::

   <sample>.kraken.gz
   <sample>.kreport
   all_samples.krakenuniq.txt

MetaPhlAn
----------
:Tool: `MetaPhlAn`_
:Output folder: ``metaphlan``

Please refer to the official `MetaPhlAn`_ installation instructions on how to
install the bowtie2 database in a separate directory outside the conda
environment.

Run `MetaPhlAn`_ on the trimmed and filtered reads to produce a taxonomic profile.
Outputs four files per sample, plus four summaries for all samples::

    <sample>.bowtie2.bz2
    <sample>.sam.bz2
    <sample>.metaphlan.krona
    <sample>.metaphlan.txt
    
    all_samples.metaphlan.krona.html
    all_samples.Species_top50.pdf
    all_samples.metaphlan.txt
    area_plot.metaphlan.pdf

The file called ``all_samples.Species_top50.pdf`` contains a clustered heatmap
plot showing abundances of the top 50 species across all samples. The taxonomic
level and the top ``N`` can be adjusted in the config. The file called
``area_plot.metaphlan.pdf`` is an areaplot summarizing the samples.

StrainPhlAn
------------
:Tool: `StrainPhlAn`_
:Output folder: ``strainphlan``

Run `StrainPhlAn`_ on the ``<sample>.sam.bz2`` output from MetaPhlAn. Generates
these output files of primary interest::

    available_clades.txt
    RAxML_bestTree.{clade_of_interest}.StrainPhlAn3.tre
    {clade_of_interest}.StrainPhlAn3_concatenated.aln

``{clade_of_interest}`` is specified in ``config.yaml``. The outputs are a tree
and an alignment file. Note that StrainPhlAn uses output generated from MetaPhlan
and will thus also need to run MetaPhlAn-associated steps (even if it is not
set to ``True`` in ``config.yaml``).

If the pipeline fails the most common issue will be that
``{clade_of_interest}`` cannot be detected in sufficient quantities in a
sufficient number of samples. Review the file ``available_clades.txt`` to
reveal which clades can be investigated in your set of samples. 


Functional profiling
**************

HUMAnN
----------
:Tool: `HUMAnN`_
:Output folder: ``humann``

Run `HUMAnN`_ on the trimmed and filtered reads to produce a functional profile.
Outputs five files per sample, plus three summaries for all samples::

    <sample>.genefamilies_{unit}.txt
    <sample>.genefamilies.txt
    <sample>.pathabundance_{unit}.txt
    <sample>.pathabundance.txt
    <sample>.pathcoverage.txt
    
    all_samples.humann_genefamilies.txt
    all_samples.humann_pathcoverage.txt
    all_samples.humann_pathabundances.txt

``{unit}`` refers to the normalization method specified in ``config.yaml``, 
the default unit is counts per million (cpm).

.. note::

   HUMAnN uses the taxonomic profiles produced by MetaPhlAn as input,
   so all MetaPhlAn-associated steps are run regardless of whether it is actually
   enabled in ``config.yaml`` or not.

HUMAnN requires large amounts of temporary disk space when processing a sample
and will automatically use a suitable temporary directory from system
environment variable ``$TMPDIR``, using Snakemake's resources feature to
evaluate the variable at runtime (which means it can utilize node-local
temporary disk if executing on a compute cluster).


Mappers
*******
|full_name| allows the use of regular read mapping tools to map the quality
controlled reads to any reference database. All mappers can be used to map
reads against several different databases (see :ref:`Mapping to multiple
databases` below). In addition, all mappers can optionally summarize read
counts per annotated feature via one of two options: 1) supplying a two-column
tab-separated annotation file with one annotation per reference sequence, or 2)
supplying a GTF or SAF format annotation file for features on the reference
sequences. Option number 1 uses a custom Python script
(``scripts/make_count_table.py``) to merge read counts per annotation, which
works well for annotations as large as your memory allows, and option number 2
uses `featureCounts`_ to summarize read counts per annotated feature. Option
number 2 is more flexible and fairly fast for typical annotation scenarios, but
might not work when the number of unique features is much lower than the number
of reference sequences. Read more about these alternatives in :ref:`Summarizing
read counts` below.

.. note::
   The mapper modules are great for mapping reads to databases with e.g.
   antibiotic resistance genes (like `megares`_) or other functionally
   annotated genes of interest.

BBMap
-----
:Tool: `BBMap`_
:Output folder: ``bbmap/<database_name>``

This module maps read using `BBMap`_. The output is in sorted and indexed BAM
format (with an option to keep the intermediary SAM file used to create the
BAM). It is possible to configure the mapping settings almost entirely
according to preference, with the exception of changing the output format. Use
the configuration parameter ``bbmap:extra`` to add any standard BBMap
commandline parameter you want.

Bowtie2
-------
:Tool: `Bowtie2`_
:Output folder: ``bowtie2/<database_name>``

This module maps read using `Bowtie2`_. The output is in BAM format. It
is possible to configure the mapping settings almost entirely according to
preference, with the exception of changing the output format from BAM.
Use the configuration parameter ``bowtie2:extra`` to add any standard Bowtie2
commandline parameter you want.

Mapping to multiple databases
-----------------------------
Note that the configuration settings of all mapper modules are slightly
different from the configuration settings from most other modules. They are
defined as lists in ``config.yaml``, e.g. (note the leading ``-`` that
signifies a list)::

    bbmap:
        - db_name: ""
          db_path: ""
          min_id: 0.76
          keep_sam: False
          keep_bam: True
          extra: ""
          counts_table:
              annotations: ""
          featureCounts:
              annotations: ""
              feature_type: ""
              attribute_type: ""
              extra: ""

This makes it possible to map the reads against several databases, each with
their own mapping options and/or custom annotations. To map against more than
one database, just create another list item underneath, containing all the same
configuration options, but with different settings. For example, to map against
``db1`` and ``db2`` with different annotation files for each::

    bbmap:
        - db_name: "db1"
          db_path: "/path/to/db1"
          min_id: 0.76
          keep_sam: False
          keep_bam: True
          extra: ""
          counts_table:
              annotations: ""
              columns: ""
          featureCounts:
              annotations: ""
              feature_type: ""
              attribute_type: ""
              extra: ""
        - db_name: "db2"
          db_path: "/path/to/db2"
          min_id: 0.76
          keep_sam: False
          keep_bam: True
          extra: ""
          counts_table:
              annotations: "/path/to/db2/annotations.txt"
              columns: "Genus,Phylum"
          featureCounts:
              annotations: ""
              feature_type: ""
              attribute_type: ""
              extra: ""


Summarizing read counts
------------------------

make_count_table.py
...................
:Tool: ``make_count_table.py``
:Output folder: ``<mapper>/<database_name>``

A custom Python script that produces tab-separated count tables with one row
per annotation, and one column per sample. The input is an annotation file that
consists of at least two tab-separated columns. The first line is a header line
with column names (must not contain spaces and avoid strange characters). Here
is an example of column names:: 

    Reference
    Annotation1
    Annotation2
    ...
    AnnotationN

The column names doesn't matter, but the names defined in the annotation file
can be used to select a subset of columns to summarize read counts for (see
more below). The first column should contain the FASTA header for each
reference sequence in the reference database used in the mapping. The count
table script truncates the header string at the first space (because Bowtie2
does this automatically it's easier to just always do it). In practice, since
the script performs truncation of headers, it doesn't matter which mapper was
used or if the annotation file contains entire headers or only the truncated
headers, as long as the bit up until the first space in each reference header
is unique. The script sums counts for each annotation for each sample. 

One parameter for the count summarization is which columns in the annotation
file to summarize on. The column names need to be assigned as a string of
comma-separated column names. They must match exactly to the column names
defined in the annotation file. This is configured in ``config.yaml``. The
script outputs one file per column, with output filename matching
``counts.<column_name>.txt``. The count table feature is activated by entering
an annotation filename in the relevant section of the configuration file,
e.g.::

    bbmap:
        counts_table:
            annotations: "path/to/annotations.tab"
            columns: "Species,Genus,taxid"


featureCounts
.............
:Tool: `featureCounts`_
:Output folder: ``<mapper>/<database_name>``

This uses `featureCounts`_ to summarize read counts per annotation and sample.
The input is a file in `GTF format`_ (or `SAF format`_, read more below).
`featureCounts`_ can summarize read counts on any feature (or meta-feature)
that is defined in your GTF file. Use the featureCounts ``attribute_type`` to
summarize read counts for any attribute defined in your GTF file. To use
`featureCounts`_ to summarize read counts, enter an annotation filename in the
configuration file, e.g.::

    bowtie2:
        featureCounts:
            annotations: "path/to/annotations.gtf"

The featureCounts module outputs several files::

    all_samples.featureCounts
    all_samples.featureCounts.summary
    all_samples.featureCounts.table.txt

The first two files are the default output files from `featureCounts`_, and the
third file is a simplified tab-separated table with count per annotation, in a
format similar to the one described for ``make_count_table.py`` above.

It is also possible to use the simplified annotation format instead of GTF. To
tell `featureCounts`_ you are using a SAF file, add ``-F SAF`` to the
featureCounts ``extra`` configuration setting, e.g.::
    
    bowtie2:
        featureCounts:
            extra: "-F SAF"

Assembly
********

StaG does not offer an assembly workflow at this time.
