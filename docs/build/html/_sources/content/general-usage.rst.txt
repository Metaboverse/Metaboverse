.. _general_link:

#############
General Usage
#############

| XPRESSpipe can be run essentially from beginning to end as a pipeline, or as individual sub-modules. We will describe each option in more detail in each section of the documentation. The purpose of XPRESSpipe is to automate the alignment, quality control, and initial analysis of single-end (SE), paired-end (PE), and ribosome profiling data. It is intended that input data is in its own directory and that each file is a properly formatted :data:`.fastq` file. However, the suffix for these files can be :data:`.fq` or :data:`.txt` as well. They can be zipped (:data:`.zip` or :data:`.gz`) or unzipped. When using intermediate sub-modules, such as :data:`align` or :data:`readDistribution`, input will vary and is explicated
in the :data:`--help` menu for each sub-module.

| Further analysis on the resulting datasets can be performed using `XPRESSplot <https://github.com/XPRESSyourself/XPRESSplot>`_.

======================================
File Naming
======================================
| In order for many of the XPRESSpipe functions to perform properly and for the output to be reliable after alignment (except for generation of a raw counts table), file naming conventions must be followed.
| 1. Download your raw sequence data and place in a folder -- this folder should contain all the sequence data and nothing else.

| 2. If you are working with single-end data, the files must be a FASTQ-formatted file and end with the suffix :data:`.fastq`, :data:`.fastq.gz`, :data:`.fq`, :data:`.fq.gz`, :data:`.txt`, :data:`.txt.gz`. We recommend the :data:`.fastq` or :data:`.fastq.gz` suffix.


| 3. If you are working with paired-end data, the rules from :data:`Step 2` apply, but must the suffix must be prefaced by the paired read group number as below:

.. code-block: shell

    ExperimentName_Rep1_a_WT.r1.fastq.gz
    ExperimentName_Rep1_a_WT.r2.fastq.gz
    ExperimentName_Rep2_a_WT.r1.fastq.gz
    ExperimentName_Rep2_a_WT.r2.fastq.gz

or

.. code-block: shell

    ExperimentName_Rep1_a_WT.read1.fastq.gz
    ExperimentName_Rep1_a_WT.read2.fastq.gz
    ExperimentName_Rep2_a_WT.read1.fastq.gz
    ExperimentName_Rep2_a_WT.read2.fastq.gz

===========
Data Output
===========
Running :data:`seRNAseq`, :data:`peRNAseq`, or :data:`riboseq` will output all intermediate and final data files as shown in this schematic:

.. image:: xpresspipe_overview.png
   :width: 600
   :align: center
