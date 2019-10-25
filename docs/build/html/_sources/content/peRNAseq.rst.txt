############################
Paired-End RNA-seq Pipeline
############################
| The following pipeline will pre-process, align, and quality check paired-end RNA-seq samples using the sub-modules discussed in earlier chapters. For more detailed information concerning these steps, please refer to the appropriate chapter.

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path\>, --input \<path\>`
     - Path to input directory -- if paired-end, file names should be exactly the same except for :data:`r1/r2.fastq` or similar suffix
   * - :data:`-o \<path\>, --output \<path\>`
     - Path to output directory
   * - :data:`-r \<path\>, --reference \<path\>`
     - Path to parent organism reference directory
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to GTF used for alignment quantification (only used for HTSeq quantification)
   * - :data:`-e <experiment_name>`, :data:`--experiment <experiment_name>`
     - Experiment name

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Optional Arguments
     - Description
   * - :data:`--two-pass`
     - Use a two-pass STAR alignment for novel splice junction discovery
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
   * - :data:`--no_multimappers>`
     - Include flag to remove multimapping reads to be output and used in downstream analyses
   * - :data:`--deduplicate`
     - Include flag to quantify reads with de-duplication (will search for files with suffix :data:`_dedupRemoved.bam`)
   * - :data:`--output_bed`
     - Include flag to output BED files for each aligned file
   * - :data:`-c <method>`, :data:`--quantification_method <method>`
     - Specify quantification method (default: htseq; other option: cufflinks. If using Cufflinks, no downstream sample normalization is required)
   * - :data:`--feature_type \<feature\>`
     - Specify feature type (3rd column in GTF file) to be used if quantifying with htseq (default: CDS)
   * - :data:`--stranded \<fr-unstranded/fr-firststrand` :data:`/fr-secondstrand||no/yes\>`
     - Specify whether library preparation was stranded (Options before || correspond with Cufflinks inputs, options after correspond with htseq inputs)
   * - :data:`--method \<RPM, RPKM, FPKM, TPM\>`
     - Normalization method to perform (options: "RPM", "TPM", "RPKM", "FPKM") -- if using either TPM, RPKM, or FPKM, a GTF reference file must be included
   * - :data:`--vcf \</path/to/file.vcf\>`
     - Provide full path and file name to VCF file if you would like detect personal variants overlapping alignments
   * - :data:`--batch \</path/filename.tsv\>`
     - Include path and filename of dataframe with batch normalization parameters
   * - :data:`--sjdbOverhang \<sjdbOverhang_amount\>`
     - Specify length of genomic sequences for constructing splice-aware reference. Ideal length is :data:`read length - 1`, so for 2x100bp paired-end reads, you would use 100 - 1 = 99. However, the default value of :data:`100` should work in most cases
   * - :data:`--mismatchRatio \<mismatchRatio\>`
     - Alignment ratio of mismatches to mapped length is less than this value. See STAR documentation for more information on setting this parameter
   * - :data:`--seedSearchStartLmax \<seedSearchStartLmax\>`
     - Adjusting this parameter by providing a lower number will improve mapping sensitivity (recommended value = 15 for reads ~ 25 nts). See STAR documentation for more information on setting this parameter
   * - :data:`genome_size`
     - Only needs to be changed if this argument was provided curing reference building AND using a two-pass alignment. This should be the length of the organism's genome in nucleotides
   * - :data:`-m <processors>, --max_processors <processors>`
     - Number of max processors to use for tasks (default: No limit)

| Run the following for more details:

.. ident with TABs
.. code-block:: python

  $ xpresspipe peRNAseq --help

-----------
Examples
-----------
| **Example 1 -- Run pipeline on paired-end RNA-seq sample files**

.. ident with TABs
.. code-block:: python

  $ xpresspipe peRNAseq \
                -i pe_test \
                -o pe_out \
                -r pe_reference \
                --gtf transcripts.gtf \
                -e pe_test \
                -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                --method TPM \
                --sjdbOverhang 100
