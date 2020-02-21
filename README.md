# ![Metaboverse](https://raw.githubusercontent.com/Metaboverse/Metaboverse/master/docs/images/metaboverse_banner.png)

A Biological and Metabolic Networks Analyzer

`WARNING: This repository is currently under development. When the software is complete and stable enough for beta testing, the alpha tag will be removed.`   
To receive updates, join our mailing list by emailing `metaboverse@gmail.com`

[![Build Status](https://travis-ci.org/Metaboverse/Metaboverse.svg?branch=master)](https://travis-ci.org/Metaboverse/Metaboverse)
[![Documentation Status](https://readthedocs.org/projects/metaboverse/badge/?version=latest)](https://metaboverse.readthedocs.io/en/latest/?badge=latest)

## Demo
A demo version of Metaboverse can be found for MacOS [here](https://github.com/Metaboverse/Metaboverse/releases/download/metaboverse-v0.0.1b/Metaboverse-darwin-x64-demo.zip). Right-click on the Metaboverse icon within the `.zip` folder to launch Metaboverse. A mock human metabolic network is available as a `.json` file with this `.zip` folder. This can be dragged and dropped into the Metaboverse home page to explore the metabolic networks and get an idea for some of the functionality within Metaboverse. However, much more is planned.

## What does Metaboverse do?
### tl;dr
- Curates metabolic model
- Maps networks or sub-networks interactively
- Overlays selected network with transcriptomics, proteomics, and/or metabolomics data
- Performs trend analysis across networks
- Performs nearest-neighbor analysis

### !tl;dr
```
Metabolism is a complex network of chemical and enzymatic reactions; however, the past century has seen
a largely reductionist approach to tackling this system. While for decades this was necessary due to
technological limitations, we are now living in the computer age, where entire cells and populations of
cells can be modeled and explored. While scientific fields, such as RNA biology and metabolism, have
experienced massive strides in recent decades with advent of RNA-seq and metabolomics, our ability to
contextualize these massive amounts of data continues to lag. This is problematic as these experiments
are often expensive and time-consuming to produce, yet we only use a fraction of the total data made
available by the experiment. In order to address these limitations, we introduce Metaboverse, a
computational analysis framework for contextualizing -omics datasets within customizable metabolic
network representations. This framework will allow for static and dynamic time-course exploration of
these datasets, and importantly will help contextualize the role of low-expressed analytes within a
network that are normally ignored due to their classification as not “differentially expressed” by most
algorithms. This tool will revolutionize our ability to more holistically understand temporal metabolic
shifts and gene-metabolite inter-cooperativity, as well as ensure we are obtaining the maximum
information from these datasets as possible.
```

## Getting started

### Requirements
- An internet connection for network curation
- An internet browser (Chrome browser preferred)
- Python
- PyPi

### Installation
0. For alpha-testing, download from this address:
```
https://github.com/Metaboverse/Metaboverse/archive/master.zip
```
1. Navigate to the following web address:
```
https://github.com/Metaboverse/Metaboverse/releases
```
2. Download the latest `.zip` release
3. Unzip the downloaded folder
4. Navigate into the folder and double-click the file named `RUN`

### Other information
- If a page seems to be malfunctioning during the alpha release stage, please press `Ctrl + r` to refresh the page, then try the operation again
