###################
Curating References
###################
| In order to quantify transcription levels from RNA-Seq data, reads must be mapped to a reference genome or transcriptome. While there are multiple alignment software packages available, XPRESSpipe performs this step using a current version of `STAR <https://github.com/alexdobin/STAR>`_ for several reasons:
| - **Splice Junction Aware**: STAR is capable of mapping reads spanning a splice junction, where more traditional packages, such as Bowtie, are incapable of doing so and are better suited for tasks such as genome alignment.
| - **Performance**: While computationally greedy (a human genome alignment requires upwards of 30 Gb RAM), the `performance and accuracy is excellent compared to the majority of other splice-aware aligners currently available <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/>`_
| - **Standard**: The foundation of the pipeline used in XPRESSpipe is based in the `TCGA <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>`_ standards for RNA-Seq alignment. This method utilizes a guided or 2-pass alignment program. In the guided alignment, a GTF with annotated splice junctions is used to guide the alignments over splice juntions. In the 2-pass alignment, reads are mapped across the genome to identify novel splice junctions. These new annotations are then incorporated into the reference index and reads are re-aligned with this new reference. While more time-intensive, this step can aid in aligning across these junctions, especially in organisms where the transcriptome is not as well annotated. If mapping to a well-documented organism, this step can be forgone and STAR will use the GTF annotations to determine intronic regions of transcripts for read mapping.

=================================
XPRESSpipe Reference Requirements
=================================
| An XPRESSpipe compatible reference directory must meet some requirements:
| - All chromosomal genome fasta files are in their own directory within the parent reference directory. If a FASTA file with all chromosomes combined is available for your organism, this can be provided, but must be in its own directory.
| - A sub-directory, named :data:`genome`, contains the STAR reference files. If :data:`createReference` is used to curate the reference, and the parent reference directory was provided as output location, this directory creation and file formatting will be handled automatically by XPRESSpipe.
| - A transcript reference (GTF), is located in the reference parent directory and is named :data:`transcripts.gtf`. If a coding-only or truncated reference GTFs are desired for read quantification, these should also be in this directory (:data:`truncate` will handle file naming and formatting so long as the output location is specified as this parent directory). This file will then need to be specified within an XPRESSpipe pipeline.
| **A completed reference directory can be created that follows these requirements by creating a directory, placing the transcripts.gtf and genomic chromosome fasta files in the parent directory and running :data:`curateReference` as described below**

============================
Get Sequence Files
============================
| The following is an example of how to get the reference files needed for generating a human reference:

.. code-block:: shell

  $ mkdir human_reference
  $ mkdir human_reference/genome_fasta
  $ cd human_reference/
  $ curl ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz -o transcripts.gtf.gz
  $ for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do curl -O ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz; done
  $ gzip -d *.gz
  $ mv *fasta genome_fasta

| - The chromosome IDs may vary depending on your organism.


============================================
Perform Full Reference Curation
============================================
| The following will create a XPRESSpipe-formatted reference directory containing all STAR reference files and transcript references needs for quantification and meta-analysis.
| A parent reference directory containing the transcripts.gtf file and all chromosomal genome fasta files must be present.
| *More details as to what each specific parameter is doing can be found in the sections below.*

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe curateReference --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-f \<path\>, --fasta \<path\>`
     - Path to genome fasta files (file names should end in .fa, .fasta, or .txt and no other files should exist in the directory with similar extensions)
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to transcript reference file names 'transcripts.gtf'

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-l, --longest_transcript`
     -  Provide argument to keep only longest transcript per gene record (RECOMMENDED)
   * - :data:`-p, --protein_coding`
     -  Provide argument to keep only gene records annotated as protein coding genes
   * - :data:`-t, --truncate`
     -  Provide argument to truncate gene records
   * - :data:`--truncate_5prime <amount>`
     -  Amount to truncate from 5' end of each transcript, requires --truncate argument (default: 45)
   * - :data:`--truncate_3prime <amount>`
     -  Amount to truncate from 3' end of each transcript, requires --truncate argument (default: 15)
   * - :data:`--sjdbOverhang \<value\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases
   * - :data:`--genome_size \<int\>`
     - If mapping to an organism with a small genome, provide the length in nucleotides. If you are not sure your organism has a small genome, provide the number of bases and XPRESSpipe will decide if this parameter needs to be changed during runtime
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Create XPRESSpipe-formatted reference for single-end alignment:**
| - Creates a star reference for single-end read mapping (1x50bp reads)
| - Keeps the longest transcript for each gene record
| - Keeps only protein_coding annotated transcripts
| - Truncates the first 45 nucleotides from the first exon of every transcript (default)
| - Truncates the last 15 nucleotides from the last exon of every transcript (default)

.. code-block:: shell

  $ xpresspipe curateReference -o /path/to/se/ref/ -f /path/to/se/ref/ -g /path/to/se/ref/transcripts.gtf --longest_transcript --protein_coding --truncate --sjdbOverhang 49

| **Example 2 -- Create refFlat files:**
| - Creates a star reference for paired-end read mapping (2x100bp reads)
| - No modifications are made to the GTF file
| - Processes are limited to 10 cores

.. code-block:: shell

  $ xpresspipe curateReference -o /path/to/pe/ref/ -f /path/to/pe/ref/ -g /path/to/pe/ref/transcripts.gtf -m 10



==========================
STAR Reference Curation
==========================
| The following creates a STAR reference compatible with XPRESSpipe. These files are output in a directory created during curation called :data:`genome` in the specified :data:`--output` directory.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe makeReference --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-f \<path\>, --fasta \<path\>`
     - Path to genome fasta files (file names should end in .fa, .fasta, or .txt and no other files should exist in the directory with similar extensions)
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to transcript reference file names 'transcripts.gtf (DO NOT USE MODIFIED GTF HERE)'

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--sjdbOverhang \<int\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases
   * - :data:`--genome_size \<int\>`
     - If mapping to an organism with a small genome, provide the length in nucleotides. If you are not sure your organism has a small genome, provide the number of bases and XPRESSpipe will decide if this parameter needs to be changed during runtime
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Create a single-end sequencing reference:**
| - Paths to output and location of genome fasta files for each chromosome are provided, as well as path and file name to transcripts.gtf file
| - Default number of threads are used for preparing reference

.. code-block:: shell

  $ xpresspipe makeReference -o /path/to/reference/ -f /path/to/reference/ -g /path/to/reference/transcripts.gtf --sjdbOverhang 49

| **Example 2 -- Create a paired-end sequencing reference:**
| - 12 threads are specified for reference creation
| - The as 2x100bp paired-end sequencing was used, the default value for :data:`--sjdbOverhang` of :data:`100` is appropriate in this case

.. code-block:: shell

  $ xpresspipe makeReference -o /path/to/reference/ -f /path/to/reference/ -g /path/to/reference/transcripts.gtf -t 12

| **Example 3 -- Create a single-end sequencing reference for Saccharomyces cerevisiae:**
| - Paths to output and location of genome fasta files for each chromosome are provided, as well as path and file name to transcripts.gtf file
| - Default number of threads are used for preparing reference\
| - Genome size is specified

.. code-block:: shell

  $ xpresspipe makeReference -o /path/to/reference/ -f /path/to/reference/ -g /path/to/reference/transcripts.gtf --sjdbOverhang 49 --genome_size 3000000

============================================
Transcript Reference Modification
============================================
| At times, quantification of transcripts to a modified transcript reference is desirable. Below are some examples:
| 1. As ribosomal RNA (rRNA) contamination is common in RNA-seq, even when a depletion step was performed prior to library preparation, it is sometimes desirable to not count these and other non-coding RNAs in the quantification and analysis.
| 2. During ribosome profiling library preparation, where a 5' and 3' pile-up of ribosome footprints due to slow initation and termination kinetics of footprints is common, it is suggested to `exclude the first 45-50 nucleotides from the 5' end and 15 nucleotides from the 3' end of each transcript during quantification <https://www.cell.com/cms/10.1016/j.celrep.2016.01.043/attachment/257faf34-ff8f-4071-a642-bfdb531c75b8/mmc1>`_. This command will automatically curate an Ensembl GTF to meet these demands for read quantification.
| 3. Several genes encode multiple isoforms or transcripts. During quantification, many software packages for counting reads to genes consider a read mapping to multiple transcripts of the same gene as a multi-mapper. Unless interested in alternate isoform usage, it is recommended that transcriptome reference files only contain the longest transcript for each gene.
| The :data:`modifyGTF` sub-module provides the ability to make the above-mentioned modifications to a GTF transcriptome reference file. The modified GTF file is output at the end and the filename is labeled with the modifications made. Truncations to each transcript reference are stranded-aware.

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe modifyGTF --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to reference GTF

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`-l, --longest_transcript`
     -  Provide argument to keep only longest transcript per gene record (RECOMMENDED)
   * - :data:`-p, --protein_coding`
     -  Provide argument to keep only gene records annotated as protein coding genes
   * - :data:`-t, --truncate`
     -  Provide argument to truncate gene records
   * - :data:`--truncate_5prime <amount>`
     -  Amount to truncate from 5' end of each transcript, requires --truncate argument (default: 45)
   * - :data:`--truncate_3prime <amount>`
     -  Amount to truncate from 3' end of each transcript, requires --truncate argument (default: 15)
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

-----------
Examples
-----------
| **Example 1 -- Create longest transcript-only, protein coding-only, truncated reference:**
| - Keeps the longest transcript for each gene record
| - Keeps only protein_coding annotated transcripts
| - Truncates the first 45 nucleotides from the first exon of every transcript (default)
| - Truncates the last 15 nucleotides from the last exon of every transcript (default)
| - Each modification desired must be implicitly passed to the sub-module

.. code-block:: shell

  $ xpresspipe modifyGTF -g /path/to/reference/transcripts.gtf --longest_transcript --protein_coding --truncate


.. _curate_link:
