# MetaboNet-Analyzer
Biological and Metabolic Networks Analyzer

`WARNING: This repository is currently under development. When the software is complete and stable, the beta tag will be removed.`

## What does MetaboNet-Analyzer do?
### tl;dr
- Curates metabolic model
- Maps networks or sub-networks interactively
- Overlays selected network with transcriptomics, proteomics, and/or metabolomics data
- Performs trend analysis across networks

### For those more details focused
```
Metabolism is a complex network of chemical and enzymatic reactions; however, the past century has seen a
largely reductionist approach to tackling this system. While for decades this was necessary due to
technological limitations, we are now living in the computer age, where entire cells and populations of
cells can be modeled and explored. While scientific fields, such as RNA biology and metabolism, have
experienced massive strides in recent decades with advent of RNA-seq and metabolomics, our ability to
contextualize these massive amounts of data continues to lag. This is problematic as these experiments are
often expensive and time-consuming to produce, yet we only use a fraction of the total data made available
by the experiment. In order to address these limitations, we introduce MetaboNet-Analyzer, a computational
analysis framework for contextualizing -omics datasets within customizable metabolic network
representations. This framework will allow for static and dynamic time-course exploration of these
datasets, and importantly will help contextualize the role of low-expressed analytes within a network that
are normally ignored due to their classification as not “differentially expressed” by most algorithms.
This tool will revolutionize our ability to more holistically understand temporal metabolic shifts and
gene-metabolite inter-cooperativity, as well as ensure we are obtaining the maximum information from these
datasets as possible.
```

## Getting started

### Installation
```
pip install MetaboNet-Analyzer
```

### Getting Started
Requires:
- tab- or comma-separated library-normalized data table
- samples metadata table
