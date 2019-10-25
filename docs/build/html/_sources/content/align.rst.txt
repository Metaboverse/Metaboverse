############################
Alignment
############################
| In order to quantify transcription on a transcript to transcript basis, individual reads called during sequencing must be mapped to the genome. While there are multiple alignment software packages available, XPRESSpipe uses a current version of `STAR <https://github.com/alexdobin/STAR>`_ to perform this step in transcription quantification for several reasons:
| - Performance: While computationally greedy (a human genome alignment requires upwards of 30 Gb RAM), the `performance and accuracy is superior to the majority of other splice aware aligners currently available <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/>`_
| - Splice Junction Aware: STAR is capable of mapping reads spanning a splice junction, where more traditional packages, such as Bowtie, are incapable of doing so and are better suited for tasks such as genome alignment.
| - Standard: The foundation of the pipeline used in XPRESSpipe is based in the `TCGA <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>`_ standards for RNAseq alignment. This method utilizes a guided or 2-pass alignment program. In the guided alignment, a GTF with annotated splice junctions is used to guide the alignments over splice juntions. In the 2-pass alignment, reads are mapped across the genome to identify novel splice junctions. These new annotations are then incorporated into the reference index and reads are re-aligned with this new reference. While more time-intensive, this step can aid in aligning across these junctions, especially in organisms where the transcriptome is not as well annotated.
| - Variant Aware: The user can provide a VCF, such as those provided by `ClinVar <ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/>`_ and `dbSNP <ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/>`_. These files are useful in integrating information about common or disease nucleotide variants into the RNA-Seq alignment stage. The files you use should match the build of the genome you are using (i.e., if using Homo Sapiens GRCh38, these builds should match between curated reference files and VCF file).

============================
Single-End RNAseq Alignment
============================
| The following runs single-end reads alignment using the specified XPRESSpipe-formatted reference directory.
| Notes:
| - For the :data:`--sjdbOverhang` argument, the same value should be entered that was used when creating the STAR reference files.
| - Ribosome profiling data can be run as a single-end RNA-seq

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe align --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-r \<path\>, --reference \<path\>`
     - Path to parent organism reference directory (must have a file called transcripts.gtf within)
   * - :data:`-t \<SE or PE\>, --type \<SE or PE\>`
     - Sequencing type ("SE" for single-end, "PE" for paired-end)

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--two-pass`
     - Use a two-pass STAR alignment for novel splice junction discovery
   * - :data:`--no_multimappers>`
     - Include flag to remove multimapping reads to be output and used in downstream analyses
   * - :data:`--deduplicate`
     - Include flag to quantify reads with de-duplication (will search for files with suffix :data:`_dedupRemoved.bam`)
   * - :data:`--vcf \</path/to/file.vcf\>`
     - Provide full path and file name to VCF file if you would like detect personal variants overlapping alignments
   * - :data:`--output_bed`
     - Include flag to output BED files for each aligned file
   * - :data:`--sjdbOverhang \<sjdbOverhang_amount\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases
   * - :data:`--mismatchRatio \<mismatchRatio\>`
     - Alignment ratio of mismatches to mapped length is less than this value. See STAR documentation for more information on setting this parameter
   * - :data:`--seedSearchStartLmax \<seedSearchStartLmax\>`
     - Adjusting this parameter by providing a lower number will improve mapping sensitivity (recommended value = 15 for reads ~ 25 nts). See STAR documentation for more information on setting this parameter
   * - :data:`genome_size`
     - Only needs to be changed if this argument was provided curing reference building AND using a two-pass alignment. Enter the size of your organism's genome in nucleotides
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Single-end RNAseq alignment:**
| - Raw reads are :data:`.fastq`-like and found in the :data:`-i /path/to/input/files/` directory. Can be uncompressed or compressed via :data:`.gz` or :data:`.zip`
| - A general output directory has been created, :data:`-o riboseq_out/`
| - :data:`--type` is specified as 'SE' and path to parent reference directory is provided
| - The value for :data:`--sjdbOverhang` used in reference creation is provided. Failure to do so will trigger an error
| - BED and BIGWIG files will be output in their own directories in :data:`output`
| - All other arguments use the default value

.. code-block:: shell

  $ xpresspipe align -i /path/to/input/files/ -o riboseq_out/ -t SE -r /path/to/reference/ --sjdbOverhang 49 --output_bed --output_bigwig

============================
Paired-End RNAseq Alignment
============================
| The following runs paired-end reads alignment using the specified XPRESSpipe-formatted reference directory.
| Notes:
| - For the :data:`--sjdbOverhang` argument, the same value should be entered that was used when creating the STAR reference files.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe align --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-r \<path\>, --reference \<path\>`
     - Path to parent organism reference directory
   * - :data:`-t \<SE or PE\>, --type \<SE or PE\>`
     - Sequencing type ("SE" for single-end, "PE" for paired-end)

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--output_bed`
     - Include flag to output BED files for each aligned file
   * - :data:`--output_bigwig`
     - Include flag to output bigwig files for each aligned file
   * - :data:`--sjdbOverhang \<sjdbOverhang_amount\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Paired-end RNAseq alignment:**
| - Raw reads are :data:`.fastq`-like and found in the :data:`-i pe_test/` directory. Can be uncompressed or compressed via :data:`.gz` or :data:`.zip`
| - A general output directory has been created, :data:`-o pe_out/`
| - :data:`--type` is specified as 'PE' and path to parent reference directory is provided
| - The value for :data:`--sjdbOverhang` used in reference creation is provided. Failure to do so will trigger an error. In this case, since the reference was created using default values, the optional flag is not used
| - BED and BIGWIG files are not output
| - All other arguments use the default value

.. code-block:: shell

  $ xpresspipe align -i /path/to/input/files/ -o riboseq_out -t PE -r /path/to/reference/
