.. _BBCountUnique: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/calcuniqueness-guide/
.. _BBDuk:  https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
.. _BBMap: https://sourceforge.net/projects/bbmap/
.. _Centrifuge: https://ccb.jhu.edu/software/centrifuge/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Kaiju: http://kaiju.binf.ku.dk/
.. _groot: https://groot-documentation.readthedocs.io
.. _MetaPhlAn2: https://bitbucket.org/biobakery/metaphlan2/
.. _featureCounts: http://bioinf.wehi.edu.au/featureCounts/
.. _HUMAnN2: https://bitbucket.org/biobakery/humann2/

Modules
=======
|full_name| is a workflow framework that connects several other tools. The
basic assumption is that all analyses start with a quality control of the
sequencing reads (using `FastQC`_), followed by human sequence removal (using
`BBMap`_). This section of the documentation aims to describe useful details
about the separate tools that are used in |full_name|.

The following subsections describe the function of each module included in
|full_name|.  Each tool produces output in a separate subfolder inside the
output folder specified in the configuration file.  The output folders
mentioned underneath each module heading refers to the subfolder inside the
main output folder specified in the configuration file.

Pre-processing
**************

qc_reads
--------------
:Tools: `FastQC`_, `BBDuk`_
:Output folder: ``fastqc``, ``trimmed_qa``

The quality control module uses `FastQC`_ to produce HTML reports of the
quality of the input reads. The quality control module also trims adapter
sequences and performs quality trimming of the input reads. The quality assured
reads are output into ``trimmed_qa``.

remove_human
--------------
:Tool: `BBMap`_
:Output folder: ``filtered_human``

The ``remove_human`` module uses `BBMap`_ to map reads against a specially
filtered and masked version of the human genome to remove reads matching to the
human genome. The output is a pair of paired-end FASTQ files, plus a single
interleaved FASTQ file with all reads that matched the human reference.


assess_depth
--------------
:Tool: `BBCountUnique`_
:Output folder: ``bbcountunique``

The ``assess_depth`` module uses `BBMap`_'s very fast kmer counting algorithms
to produce saturation curves. The saturation curve essentially shows a
histogram of the proportion of unique kmers observed per reads processed, and
can be used to assess how deep a sample has been sequenced. The module outputs
one plot per sample.


Naive sample comparison
***********************

sketch_compare
--------------
:Tools: ``sketch.sh``, ``comparesketch.sh`` from `BBMap`_
:Output folder: ``sketch_compare``

The ``sketch_compare`` module uses the very fast MinHash implementation in
`BBMap`_ to compute MinHash sketches of all samples to do an all-vs-all
comparison of all samples based on their kmer content. The module outputs
gzip-compressed sketches for each sample, as well as a heatmap plot showing the
overall similarity of all samples.


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
(``make_count_table.py``) to merge read counts per annotation which works well
for annotations as large as your memory allows, and option number 2 uses
`featureCounts`_ to summarize read counts per annotated feature. Option number
2 is more flexible and normally faster for typical annotation scenarios, but
might not work when the number of unique features is much lower than the number
of reference sequences. Read more about these alternatives in :ref:`Summarizing
read counts` below.

BBMap
-----
:Tool: `BBMap`_
:Output folder: ``bbmap/<database_name>``

This module maps read using `BBMap`_. The output is in gzipped SAM format. It
is possible to configure the mapping settings almost entirely according to
preference, with the exception of changing the output format from gzipped SAM.
Use the configuration parameter ``bbmap:extra`` to add any standard BBMap
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
          extra: ""
          counts_table:
              annotations: ""
          featureCounts:
              annotations: ""
              feature_type: ""
              attribute_type: ""
              extra: ""
        - db_name: "db2"
          db_path: "/path/to/db2"
          min_id: 0.76
          extra: ""
          counts_table:
              annotations: "/path/to/db2/annotations.txt"
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

A custom Python script produces a tab-separated count table with one row per
annotation, and one column per sample. The input is an annotation file that
consists of two tab-separated columns (no header)::

    Reference sequence
    Annotation

The script sums counts for each annotation for each sample. The output filename
is ``all_samples.counts_table.tab``. To produce the output table, enter an
annotation filename in the configuration file, e.g.::

    bbmap:
        counts_table:
            annotations: "path/to/annotations.tab"

featureCounts
.............
:Tool: `featureCounts`_
:Output folder: ``<mapper>/<database_name>``

This uses the well-known `featureCounts`_ to summarize read counts per
annotation and sample. The input is a file in `GTF format`_ (or `SAF format`_,
read more below). `featureCounts`_ can summarize read counts on any feature (or
meta-feature) that is defined in your GTF file. Use the featureCounts
``attribute_type`` to summarize read counts for any attribute defined in your
GTF file. To use `featureCounts`_ to summarize read counts, enter an annotation
filename in the configuration file, e.g.::

    bowtie2:
        featureCounts:
            annotations: "path/to/annotations.tab"

The featureCounts module outputs several files::

    all_samples.featureCounts
    all_samples.featureCounts.summary
    all_samples.featureCounts.table.tsv

The first two files are the default output files from `featureCounts`_, and the
third file is a simplified tab-separated table with count per annotation, in a
format similar to the one described for ``make_count_table.py`` above.

It is also possible to use the simplified annotation format instead of GTF. To
tell `featureCounts`_ you are using a SAF file, add ``-F SAF`` to the
featureCounts ``extra`` configuration setting, e.g.::
    
    bowtie2:
        featureCounts:
            extra: "-F SAF"


Taxonomic profiling
*******************

Centrifuge
---------
:Tool: `Centrifuge`_
:Output folder: ``centrifuge``

Run `Centrifuge`_ on the trimmed and filtered reads to produce a taxonomic
profile.  Outputs two files per sample: ``<sample>.centrifuge_report.tsv`` and
``<sample>.centrifuge.tsv``.

Kaiju
-----
:Tool: `Kaiju`_
:Output folder: ``kaiju``

Run `Kaiju`_ on the trimmed and filtered reads to produce a taxonomic profile.
Outputs four files per sample, plus a summary HTML Krona report with the
profiles of all samples (``all_samples.kaiju.krona.html``). The four per-sample
output files are::

    <sample>.kaiju
    <sample>.kaiju.summary.family
    <sample>.kaiju.summary.genus
    <sample>.kaiju.summary.species
    <sample>.krona

MetaPhlAn2
----------
:Tool: `MetaPhlAn2`_
:Output folder: ``metaphlan2``

Run `MetaPhlAn2`_ on the trimmed and filtered reads to produce a taxonomic profile.
Outputs three files per sample, plus three summaries for all samples::

    <sample>.bowtie2.bz2
    <sample>.metaphlan2.krona
    <sample>.metaphlan2.txt
    
    all_samples.metaphlan2.krona.html
    all_samples.metaphlan2.pdf
    all_samples.metaphlan2.txt

The file called ``all_samples.metaphlan2.pdf`` contains a standard MetaPhlAn2
clustered heatmap plot containing all samples.


Functional profiling
**************

HUMAnN2
----------
:Tool: `HUMAnN2`_
:Output folder: ``humann2``

Run `HUMAnN2`_ on the trimmed and filtered reads to produce a functional profile.
Outputs three files per sample, plus three summaries for all samples::

    <sample>.genefamilies.tsv
    <sample>.pathcoverage.tsv
    <sample>.pathabundances.tsv
    
    all_samples.humann2_genefamilies.tsv
    all_samples.humann2_pathcoverage.tsv
    all_samples.humann2_pathabundances.tsv

Note that HUMAnN2 uses the taxonomic profiles produced by MetaPhlAn2, so all
MetaPhlAn2-associated steps are run regardless of whether it is actually
enabled in ``config.yaml`` or not.


Antibiotic resistance
*********************

Groot
-------
:Tool: `groot`_
:Output folder: ``groot``

Run `groot`_ to align reads to an antibiotic resistance gene database to
produce antibiotic resistance gene profiles. Outputs one subfolder per sample,
containing two files and two subfolders::

    <sample>/<sample>.groot_aligned.bam
    <sample>/<sample>.groot_report.tsv
    <sample>/<sample>/groot-graphs
    <sample>/<sample>/groot-plots

The ``<sample>.groot.bam`` file contains mapping results against all resistance
gene graphs, and the ``<sample>.groot_report.tsv`` file contains a list of
observed antibiotic resistance genes in the sample. The two subfolders contain 
all mapped graphs and coverage plots of all detected antibiotic resisatance genes.

