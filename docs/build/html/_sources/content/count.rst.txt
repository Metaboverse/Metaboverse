############################
Read Quantification
############################

===============================
Quantifying and Collating Reads
===============================
| In order to quantify aligned reads, they must be counts to a reference transcriptome. This will tell you in relative terms how much of each transcript is expressed in a system. The following sub-module will perform this quantification, as well as compile all sample quantifications into a single data matrix for downstream use.
| XPRESSpipe uses `Cufflinks <http://cole-trapnell-lab.github.io/cufflinks/>`_ as the default, but `HTSeq <https://htseq.readthedocs.io/en/release_0.11.1/>`_ can also be specified. Cufflinks is one of the most accurate read quantifiers currently available, but HTSeq is still widely used and is part of the `TCGA <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>`_ pipeline.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe count --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory of SAM files
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to GTF used for alignment quantification (if a modified GTF was created, this should be provided here; if using Cufflinks and you want isoform abundance estimates, important that you do not provide a longest transcript only GTF)

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-e <experiment_name>`, :data:`--experiment <experiment_name>`
     - Experiment name
   * - :data:`-c <method>`, :data:`--quantification_method <method>`
     - Specify quantification method (default: htseq; other option: cufflinks. If using Cufflinks, no downstream sample normalization is required)
   * - :data:`--feature_type \<feature\>`
     - Specify feature type (3rd column in GTF file) to be used if quantifying with htseq (default: CDS)
   * - :data:`--stranded \<fr-unstranded/fr-firststrand` :data:`/fr-secondstrand||no/yes\>`
     - Specify whether library preparation was stranded (Options before || correspond with Cufflinks inputs, options after correspond with htseq inputs)
   * - :data:`--deduplicate`
     - Include flag to quantify reads with de-duplication (will search for files with suffix :data:`_dedupRemoved.bam`)
   * - :data:`--bam_suffix <suffix>`
     - Change from default suffix of _Aligned.sort.bam
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Count ribosome profiling alignments:**
| - Input points to directory with SAM alignment files that are sorted by name
| - An experiment name is provided to name the final data matrix
| - Reads are quantified only to coding genes and are not counted if mapping to the first x nucleotides of each transcript exon 1 (x being the value provided for truncation when initially creating the reference files)

.. code-block:: shell

  $ xpresspipe count -i riboseq_out/alignments/ -o riboseq_out/ -r se_reference/ -g se_reference/transcripts_codingOnly_truncated.gtf -e se_test

| **Example 2 -- Count paired-end alignments:**
| - Input points to directory with SAM alignment files that are sorted by name
| - An experiment name is not provided and a default name is given to the data matrix using datatime
| - Reads are quantified to the entire transcriptome (coding and non-coding, no truncation)

.. code-block:: shell

  $ xpresspipe count -i pe_out/alignments/ -o pe_out/ -r pe_reference/
