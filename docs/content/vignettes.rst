###############
Data Vignettes
###############

| Test datasets are included with each release of Metaboverse. Once you have opened the executable, unzip the :data:`test_data.zip` file and explore the provided datasets. Each sub-directory is explained below.


===============================
:data:`test_data/human-metabolomics`
===============================
| **Source**: `Metabolomic Markers of Altered Nucleotide Metabolism in Early Stage Adenocarcinoma (Wikoff, et al., Cancer Prev Res (Phila), 2015) <https://aacrjournals.org/cancerpreventionresearch/article/8/5/410/50431/Metabolomic-Markers-of-Altered-Nucleotide>`_.
|
| **Contents**:
| - :data:`lung_tumor_vs_normal_measurements.txt`: GC-TOF metabolomics for paired tumor and normal lung tissue of the LUAD dataset (Wikoff, 2015) as used in the Metaboverse manuscript vignette. This dataset provides log2 Fold Change and Benjamini-Hochberg adjusted p-values for each measured entity at each time-point, and should be provided on the VARIABLES AND DATA page of Metaboverse.
| 
| 1) On the Curation page of Metaboverse, select *Homo sapiens* as the model organism.
| 2) On the VARIABLES AND DATA page of Metaboverse, you would specify the experiment type as "Default (2-condition)".


===============================
:data:`test_data/mouse-multiomics`
===============================
| **Source**: `STATegra, a comprehensive multi-omics dataset of B-cell differentiation in mouse (Gomez-Cabrero, et al., Sci Data, 2019) <https://www.nature.com/articles/s41597-019-0202-7>`_.
| 
| **Contents**:
| - :data:`proteomics_formatted.txt`: Proteomics for the timecourse of B-cell differentiation in mouse. This dataset provides log2 Fold Change and Benjamini-Hochberg adjusted p-values for each measured entity at each time-point, and should be provided on the VARIABLES AND DATA page of Metaboverse.
| - :data:`metabolomics_combined.txt`: A combination of GC- and LC-MS metabolomics data for the timecourse of B-cell differentiation in mouse. This dataset provides log2 Fold Change and Benjamini-Hochberg adjusted p-values for each measured entity at each time-point, and should be provided on the VARIABLES AND DATA page of Metaboverse.
| 
| 1) On the CURATION page of Metaboverse, select Mus musculus as the model organism.   
| 2) On the VARIABLES AND DATA page of Metaboverse, you would specify the experiment type as "Time-Course", and would provide labels such as "0hr, 2hr, 6hr, 12hr, 18hr, 24hr". 


===============================
:data:`test_data/fish-proteomics`
===============================
| **Source**: `Comprehensive and quantitative proteomic analyses of zebrafish plasma reveals conserved protein profiles between genders and between zebrafish and human. (Li, et al., Sci Rep, 2016) <https://www.nature.com/articles/srep24329>`_.
|
| **Contents**:
| - :data:`proteomics_male_female_log.txt`: Proteomics for male vs female zebrafish (Li, 2016). This dataset provides the log10 transformed quantifications for each biological replicate between the two measured conditions. These data should be used as practice for the "Format dataset" tool on the VARIABLES AND DATA page of Metaboverse.
| 
| 1) On the Curation page of Metaboverse, select *Danio rerio* as the model organism.
| 2) On the VARIABLES AND DATA page of Metaboverse, you would launch the "Format dataset" tool and upload the file to generate the log2 Fold Change and p-values for your comparison conditions.
| 3) On the VARIABLES AND DATA page of Metaboverse, you would specify the experiment type as "Default (2-condition)".


===============================
:data:`test_data/yeast-multiomics-timecourse`
===============================
| **Source**: `Metaboverse: Automated discovery and visualization of diverse metabolic regulatory patterns (Berg, et al., bioRxiv, 2020) <https://www.biorxiv.org/content/10.1101/2020.06.25.171850>`_.
| 
| **Contents**:
| - :data:`sce_mct1_03hr_counts_diffx.txt`: DESeq2-processed RNA-seq data for the 3 hr time-point of the MCT1 dataset as used in the Metaboverse manuscript vignette. This single time-point provides log2 Fold Change and FDR values for each measured entity, and should be provided on the VARIABLES AND DATA page of Metaboverse.
| - :data:`metabolomics_timecourse_mct1.txt`: GC-MS metabolomics for the 0min, 15min, 30min, 60min, and 180min time-points of the MCT1 dataset as used in the Metaboverse manuscript vignette. This time-series provides log2 Fold Change and Benjamini-Hochberg adjusted p-values for each measured entity at each time-point, and should be provided on the VARIABLES AND DATA page of Metaboverse.
|
| 1) On the CURATION page of Metaboverse, select *Saccharomyces cerevisiae* as the model organism.
| 2) On the VARIABLES AND DATA page of Metaboverse, you would also specify the experiment type as "Time-Course", and would provide labels such as "0min, 15min, 30min, 60min, 180min".


===============================
:data:`test_data/yeast-multiomics-singletimepoint`
===============================
| **Source**: `Metaboverse: Automated discovery and visualization of diverse metabolic regulatory patterns (Berg, et al., bioRxiv, 2020) <https://www.biorxiv.org/content/10.1101/2020.06.25.171850>`_.
| 
| **Contents**:
| - :data:`sce_mct1_12hr_counts_diffx.txt`: DESeq2-processed RNA-seq data for the 12 hr time-point of the MCT1 dataset as used in the Metaboverse manuscript vignette. This single time-point provides log2 Fold Change and FDR values for each measured entity, and should be provided on the VARIABLES AND DATA page of Metaboverse.
| - :data:`proteomics_mct1_12hr.txt`: Quantitative proteomics for the 12 hr time-point of the MCT1 dataset as used in the Metaboverse manuscript vignette. This single time-point provides log2 Fold Change and Benjamini-Hochberg adjusted p-values for each measured entity, and should be provided on the VARIABLES AND DATA page of Metaboverse.
| - :data:`mct1_12hr_metabolomics.txt`: LC-MS metabolomics for the 12 hr time-point of the MCT1 dataset as used in the Metaboverse manuscript vignette. This single time-point provides log2 Fold Change and Benjamini-Hochberg adjusted p-values for each measured entity at each time-point, and should be provided on the VARIABLES AND DATA page of Metaboverse.
|
| 1) On the CURATION page of Metaboverse, select *Saccharomyces cerevisiae* as the model organism.
| 2) On the VARIABLES AND DATA page of Metaboverse, you would also specify the experiment type as "Default (2-condition)".