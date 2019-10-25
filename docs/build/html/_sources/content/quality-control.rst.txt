##############################
Quality Control
##############################


=================================
Read Distribution Analysis
=================================
| When performing RNA-seq, your sequencing library population is important to assess to ensure a quality sequencing run. Unexpected populations can be indicative of RNA degradation or other effects. In ribosome profiling, the expected footprint size is ~28-30 nucleotides, so you would expect a peak in this region when running your analysis. The following module will run read distribution analysis for all :data:`.fastq` samples within a given directory. It is recommended this analysis be performed on trimmed reads to clean up adapters and get the true distribution of sequence reads in the library. When this is run within the pipeline, it will analyze just the post-trimming :data:`.fastq` files.

| Additionally, if running one of XPRESSpipe's pipelines, you can refer to the MultiQC :data:`html` file for general summary statistics, which include read length distributions for all samples.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe readDistribution --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory of trimmed fastq (or untrimmed fastq) files
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory

.. list-table::
  :widths: 35 50
  :header-rows: 1

  * - Required Arguments
    - Description
  * - :data:`-t \<SE or PE\>, --type \<SE or PE\>`
    - Sequencing type ("SE" for single-end, "PE" for paired-end)
  * - :data:`-e \<experiment_name\>, --experiment \<experiment_name\>`
    - Experiment name
  * - :data:`-m <processors>, --max_processors <processors>`
    - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Analyze read distributions from ribosome profiling libraries**

.. ident with TABs
.. code-block:: python

  $ xpresspipe readDistribution -i riboprof_out/trimmed_fastq -o riboprof_out -e se_test

.. image:: se_test_fastqc_summary.png
  :width: 450px

=================================
Metagene Analysis
=================================
| Analyze each sequencing sample to ensure equal distribution of reads across all transcripts. Can be useful in identifying 5' or 3' biases in sequence preparation.
| Requires a transcriptome-mapped BAM files, which can be output by `STAR <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>`_ and are automatically output during any XPRESSpipe alignment run.

.. code-block:: shell

  $ xpresspipe metagene --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory of transcriptome-mapped BAM files
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to un-modified reference GTF

.. list-table::
  :widths: 35 50
  :header-rows: 1

  * - Optional Arguments
    - Description
  * - :data:`-e \<experiment_name\>, --experiment \<experiment_name\>`
    - Experiment name
  * - :data:`--feature_type \<feature_type\>`
    - Specify feature type (3rd column in GTF file) to be used in calculating metagene coverage (default: exon; alternative: CDS)
  * - :data:`--bam_suffix \<suffix\>`
    - Change from default suffix of toTranscriptome.out.bam if transcriptome-mapped files were processed outside of XPRESSpipe
  * - :data:`-m \<processors\>, --max_processors \<processors\>`
    - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Analyze metagene profiles of sequence libraries**
| - Use default transcript reference (maps to all transcripts, even if non-coding)

.. ident with TABs
.. code-block:: python

  $ xpresspipe metagene -i riboprof_out/alignments/ -o riboprof_out -g se_reference/transcripts.gtf -e se_test

.. image:: se_test_metagene_summary.png
  :width: 450px

NOTE: As you can appreciate, there are systematic 5' biases in these library preparations. A good RNA-seq library should generally have even coverage across all transcripts.


=================================
Intron-collapsed Gene Coverage Analysis
=================================
| Plot the coverage of a given gene for a sample or set of samples with introns collapsed.

.. code-block:: shell

  $ xpresspipe geneCoverage --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory of transcriptome-aligned BAM files
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to reference GTF
   * - :data:`-n \<gene_name\>, --gene_name \<gene_name\>`
     - Gene name (case sensitive)

.. list-table::
  :widths: 35 50
  :header-rows: 1

  * - Optional Arguments
    - Description
  * - :data:`-e \<experiment_name\>, --experiment \<experiment_name\>`
    - Experiment name to save output summaries as
  * - :data:`--bam_suffix \<suffix\>`
    - Change from default suffix of toTranscriptome.out.bam if using a different BAM file
  * - :data:`--type \<type>`
    - Record type to map across (i.e. "exon", "CDS") (case-sensitive)
  * - :data:`--samples \<sample_list\> [<sample_list> ...]`
    - Provide a space-separated list of sample names to include in analysis (will only include those listed, and will plot in the order listed)
  * - :data:`--sample_names \<suffix\>`
    - Provide a space-separated list of sample names to use for labels
  * - :data:`--plot_color \<color>`
    - Indicate plotting color
  * - :data:`-m \<processors\>, --max_processors \<processors\>`
    - Number of max processors to use for tasks (default: No limit)


-----------
Examples
-----------
| **Example 1 -- Analyze gene coverage profile of sequence libraries**
| - Use default transcript reference (will generate a longest transcript-only reference)
| - Analyze SLC1A1
| - Analyze along chosen record type (default: exon, but could also use CDS if looking at ribosome profiling data)
| - Analyzing BAM files ending in :data:`.sort.bam`
| - Specifying names to use in plotting -- if not using :data:`--samples`, these files will be plotted alphabetically, so the listed order should also be alphabetical. If using :data:`--samples`, need to specify names in the same order you provided for this argument.

.. ident with TABs
.. code-block:: python

  $ xpresspipe geneCoverage -i /path/to/bam_files -o ./ -g /path/to/reference.gtf \
    -n SLC1A1 --type exon --bam_suffix .sort.bam \
    --sample_names SRR1795425 SRR1795433 SRR1795435 SRR1795437

.. image:: geneCoverage_IGV_comparison.png
  :width: 750px

NOTE: The coverage estimations use a 20 nt rolling window mean method to smoothen the coverage plots. In both A and B in the image above, the top plot was generated with IGV (https://software.broadinstitute.org/software/igv/) and the bottom with :data:`xpresspipe geneCoverage`. Green boxes show approximately the same region for comparison.



=================================
Codon Periodicity Analysis
=================================
| Analyze periodicity of most abundant read length. Useful in ribosome profiling samples for identifying that ribosomes are taking the expected 3 nucleotide steps along a transcript. If this is not apparent from the analysis, it may be indicative of poor sequence coverage of the ribosome profiling libraries.

.. code-block:: shell

  $ xpresspipe periodicity --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory of transcriptome-aligned BAM files
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to reference GTF

.. list-table::
  :widths: 35 50
  :header-rows: 1

  * - Optional Arguments
    - Description
  * - :data:`-e \<experiment_name\>, --experiment \<experiment_name\>`
    - Experiment name to save output summaries as
  * - :data:`--bam_suffix \<suffix\>`
    - Change from default suffix of toTranscriptome.out.bam if using a different BAM file
  * - :data:`-m \<processors\>, --max_processors \<processors\>`
    - Number of max processors to use for tasks (default: No limit)


-----------
Examples
-----------
| **Example 1 -- Analyze periodicity from ribosome profiling libraries**

.. ident with TABs
.. code-block:: python

  $ xpresspipe periodicity -i riboprof_out/alignments/ -o riboprof_out -g se_reference/transcripts.gtf -e se_test
