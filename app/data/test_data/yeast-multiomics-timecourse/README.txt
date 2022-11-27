USING THE EXAMPLE DATASET WITH METABOVERSE
==========================================

Included in this test dataset are the following files:
- sce_mct1_03hr_counts_diffx.txt
      DESeq2-processed RNA-seq data for the 3 hr time-point of the MCT1 dataset as used in the Metaboverse manuscript vignette. This single time-point provides log2 Fold Change and FDR values for each measured entity, and should be provided on the VARIABLES AND DATA page of Metaboverse.
- metabolomics_timecourse_mct1.txt
      GC-MS metabolomics for the 0min, 15min, 30min, 60min, and 180min time-points of the MCT1 dataset as used in the Metaboverse manuscript vignette. This time-series provides log2 Fold Change and Benjamini-Hochberg adjusted p-values for each measured entity at each time-point, and should be provided on the VARIABLES AND DATA page of Metaboverse.

1) On the CURATION page of Metaboverse, select Saccharomyces cerevisiae as the model organism.
2) On the VARIABLES AND DATA page of Metaboverse, you would also specify the experiment type as "Time-Course", and would provide labels such as "0min, 15min, 30min, 60min, 180min".

To learn more about each option, hover over the "i" icon at each step to display more information, or read the documentation for more complete descriptions of each step: https://metaboverse.readthedocs.io/
