###################
Read Pre-Processing
###################

===================
Read Trimming
===================

| Trimming is a necessary part of RNAseq data processing due to the technological limitations described below:
| - Inherent in RNA-seq library creation, RNA is fragmented and adapter sequences are ligated to the sequence. These adapters include information such as sample batch and act as a primer for the sequencer to recognize the fragment as something to analyze. However, these adapters, once sequenced, prevent alignment to a reference as large chunks of the fragment are synthetic sequence not found in the actual organism's genome/transcriptome.
| - A sequencer's job is to read a fragment base by base and determine the nucleotide species each step of the way. While the technology has greatly improved over the years, a probability of error remains. Mis-called bases can prevent proper alignment of the sequenced fragment to the reference. Therefore, it is important for low confidence base calls to be trimmed from each read.
|
| Trimming is performed by `fastp <https://github.com/OpenGene/fastp>`_

-------------
Arguments
-------------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe trim --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory -- if paired-end, file names should be exactly the same except for :data:`r1/r2.fastq` or similar suffix
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-a \<adapter1 ...\> [\<adapter1 ...\> ...]`, :data:`--adapter \<adapter1 ...\> [\<adapter1 ...\> ...]`
     - Specify adapter(s) in list of strings -- for single-end, only provide one adapter -- if :data:`None` are provided, software will attempt to auto-detect adapters -- if "POLYX" is provided as a single string in the list, polyX adapters will be trimmed. If you want to auto-detect adapters in for paired-end reads, provide :data:`None` twice
   * - :data:`-q \<PHRED_value\>, --quality \<PHRED_value\>`
     - PHRED read quality threshold (default: :data:`28`)
   * - :data:`--min_length \<length_value\>`
     - Minimum read length threshold to keep for reads (default: :data:`17`)
   * - :data:`--umi_location \<location\>`
     - Provide parameter to process UMIs -- provide location (see fastp documentation for more details, generally for single-end sequencing, you would provide 'read1' here; does not work with -a polyX option)
   * - :data:`--umi_length \<length\>`
     - Provide parameter to process UMIs -- provide UMI length (must provide the --umi_location argument); does not work with -a polyX option)
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: Max)

--------------
Examples
--------------
| **Example 1 -- Trim ribosome profiling (or single-end) sequence data using default preferences:**
| - Raw reads are :data:`.fastq`-like and found in the :data:`-i riboprof_test/` directory. Can be uncompressed or compressed via :data:`.gz` or :data:`.zip`
| - A general output directory has been created, :data:`-o riboprof_out/`
| - All other arguments use the default value

.. code-block:: shell

  $ xpresspipe trim -i riboprof_test/ -o riboprof_out/

| **Example 2 -- Predict adapter and trim ribosome profiling (or single-end) sequence data:**
| - A minimum read length of 22 nucleotides after trimming is required in order to keep the read
| - A maximum or 6 processors can be used for the task
| - The :data:`--adapters` argument was not passed, so an attempt to discover adapter sequences will be made (this is not always the most efficient or thorough method of trimming and providing the adapter sequences is recommended)

.. code-block:: shell

  $ xpresspipe trim -i riboprof_test/ -o riboprof_out/ --min_length 22 -m 6

| **Example 3 -- Pass explicit adapter trim ribosome profiling (or single-end) sequence data:**
| - The default minimum read length threshold will be used
| - The maximum number of processors will be used by default
| - The :data:`--adapters` argument was passed, so adapter sequences will trimmed explicitly

.. code-block:: shell

  $ xpresspipe trim -i riboprof_test/ -o riboprof_out/ -a CTGTAGGCACCATCAAT

| **Example 4 -- Predict adapter and trim paired-end sequence data:**
| - The :data:`--adapters` argument was passed as :data:`None None`, so an attempt to discover adapter sequences will be made for paired-end reads. The :data:`-a None None` syntax is essential for :data:`trim` to recognize the reads as paired-end

.. code-block:: shell

  $ xpresspipe trim -i pe_test/ -o pe_out/ -a None None

| **Example 5 -- Pass explicit adapter and trim paired-end sequence data:**
| - The :data:`--adapters` argument was passed, so adapter sequences will trimmed explicitly

.. code-block:: shell

  $ xpresspipe trim -i pe_test/ -o pe_out/ -a ACACTCTTTCCCTACACGACGCTCTTCCGATC GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG

| **Example 6 -- Trim single-end sequence data of polyX adapters:**
| - The :data:`--adapters POLYX` argument was passed, so adapter sequences will trimmed of polyX sequences

.. code-block:: shell

  $ xpresspipe trim -i se_test/ -o se_out/ -a POLYX
