# ![Metaboverse](https://raw.githubusercontent.com/Metaboverse/Metaboverse/master/docs/content/images/metaboverse_banner.png)

[![Build Status](https://travis-ci.org/Metaboverse/Metaboverse.svg?branch=master)](https://travis-ci.org/Metaboverse/Metaboverse)
[![Documentation Status](https://readthedocs.org/projects/metaboverse/badge/?version=latest)](https://metaboverse.readthedocs.io/en/latest/?badge=latest)
[![bioRxiv preprint](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.06.25.171850-BF2636)](https://www.biorxiv.org/content/10.1101/2020.06.25.171850v1)
[![Github All Releases](https://img.shields.io/github/downloads/Metaboverse/Metaboverse/total.svg)]()

## What does Metaboverse do?
Integrating multi- or single-omic metabolic data upon the metabolic network can be challenging for a variety of reasons. Metaboverse seeks to simplify this task for users by providing a simple, user-friendly interface for layering their data on a dynamic representation of the metabolic network. Additionally, Metaboverse provides several new tools to enable the contextualization of metabolic data:
- 90+ model organism networks
- Integrate two-condition, multi-condition, and timecourse data
- Search for regulatory events using our regulatory pattern search engine
- Collapse reactions with intermediate reactions with missing data for easier visualization and analysis of sparse datasets
- Explore super-pathway-specific reaction perturbation networks
- Explore proximal reactions to a specific entity across the global reaction network using our Nearest Neighbors search features
- And more to come!

*During the development (public beta) phase of Metaboverse, we intend to release an updated version of the software weekly as we incorporate feedback from users.*  

Detailed walkthroughs and additional usage information can be found in the [documentation](https://metaboverse.readthedocs.io/en/latest).

[![METABOVERSE](https://yt-embed.herokuapp.com/embed?v=ytTIlBKzq-c)](http://www.youtube.com/watch?v=ytTIlBKzq-c "Metaboverse Video Walkthrough")

## Getting started

### Requirements
- An internet connection for network curation
- The Metaboverse app for your operating system

### Installation
- Download the appropriate Metaboverse app `.zip` file for your operating system from [this location](https://github.com/Metaboverse/Metaboverse/releases/latest).
- Unzip the downloaded folder
- Open the `Metaboverse` app
- Please refer to the [documentation](https://metaboverse.readthedocs.io/en/latest/content/general-usage.html) for more information.

### Getting Help
If you have questions, requests, or bugs to report, please use the [Metaboverse issues forum](https://github.com/Metaboverse/Metaboverse/issues). Please clearly describe the problem, what you have tried, as well as screenshots of any error information. If possible, click on the `View` menu tab, click `Toggle Developer Tools`, click the `Console` tab of the window that opens, and take a screenshot of the output in this panel.

### Feedback
Have any feedback? Let us know [here](https://forms.gle/4z51DMnagWRvKhc38).

### Trying out Metaboverse
You can access some example network curations with biological data [here](https://github.com/Metaboverse/manuscript/tree/master/data/databases) -- download the `.json.zip` files.
With each release archive or Metaboverse, a `test_data.zip` file is included. Unzip this file and read the `README.txt` file for more information on this example dataset.

### How does the underlying code for Metaboverse work?
Metaboverse is currently segmented into two parts:
1. [`metaboverse-cli`](https://github.com/Metaboverse/metaboverse-cli): The network curation and modeling utilities that form the back-end of the Metaboverse app. For each release of Metaboverse, OS-specific binaries are generated of the backend and incorporating into the GUI app.
2. [`Metaboverse` interactive app](https://github.com/Metaboverse/Metaboverse): The platform-independent app for visualizing and exploring data on the metabolic network.

Metaboverse works by doing the following:
1. Generates an organism-specific reaction network using the [Reactome knowledgebase](https://reactome.org/)
2. References the [ChEBI](https://www.ebi.ac.uk/chebi/) and [HMDB](https://hmdb.ca/) databases to cross-reference metabolite synonyms.
3. Generates a network-based reaction model
4. Layers user data onto the global reaction network
5. Optionally broadcasts gene expression data to protein expression nodes where protein values are unavailable
6. Prevent the visualization of certain nodes to create simplified visualizations of pathways
7. Runs just-in-time searches of the global network for regulatory patterns of interest centered around a reaction
8. Generates just-in-time visualizations of global or super-pathway-specific perturbation networks
9. Generates just-in-time general visualization of canonical pathways
