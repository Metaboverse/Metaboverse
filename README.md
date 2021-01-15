# ![Metaboverse](https://raw.githubusercontent.com/Metaboverse/Metaboverse/master/docs/content/images/metaboverse_banner.png)

![Build Status](https://github.com/Metaboverse/Metaboverse/workflows/build/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/metaboverse/badge/?version=latest)](https://metaboverse.readthedocs.io/en/latest/?badge=latest)
[![bioRxiv preprint](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.06.25.171850-BF2636)](https://www.biorxiv.org/content/10.1101/2020.06.25.171850v1)
[![Github All Releases](https://img.shields.io/github/downloads/Metaboverse/Metaboverse/total.svg)](https://github.com/Metaboverse/Metaboverse/releases/)
[![DOI](https://zenodo.org/badge/203264184.svg)](https://zenodo.org/badge/latestdoi/203264184)

## What does Metaboverse do?
Integrating multi- or single-omic metabolic data upon the metabolic network can be challenging for a variety of reasons. Metaboverse seeks to simplify this task for users by providing a simple, user-friendly interface for layering their data on a dynamic representation of the metabolic network and automatically searching the network for interesting regulatory or other patterns. Additionally, Metaboverse provides several tools to enable the contextualization of metabolic data:
- 90+ model organism networks from the [Reactome knowledgebase](https://reactome.org/)
- Integrate two-condition, multi-condition, and timecourse data
- Search for regulatory events using our regulatory pattern search engine
- Collapse reactions with intermediate reactions with missing data for easier visualization and pattern analysis within sparse datasets
- Explore super-pathway-specific reaction perturbation networks
- Explore proximal reactions to a specific entity across the global reaction network using our Nearest Neighbors search features
- And more to come! We have received and implemented many valuable feature requests from the community and encourage these requests!

Detailed walkthroughs and additional usage information can be found in the [documentation](https://metaboverse.readthedocs.io/en/latest).

[![METABOVERSE](https://yt-embed.herokuapp.com/embed?v=U7m78Tbs5KE)](https://youtu.be/U7m78Tbs5KE "Metaboverse Video Walkthrough")

## Getting started

### Requirements
- An internet connection for network curation
- The most current version of the [Metaboverse app](https://github.com/Metaboverse/Metaboverse/releases) for your operating system
- A Linux/macOS/Windows 64-bit operating system
- At least 4 GB RAM and 5 GB of free storage space

### Installation
- Download the appropriate Metaboverse app `.zip` file for your operating system from [this location](https://github.com/Metaboverse/Metaboverse/releases/latest).
- Unzip the downloaded folder
- Open the `Metaboverse` app
- Please refer to the [documentation](https://metaboverse.readthedocs.io/en/latest/content/general-usage.html) for more information in using the app.
- If you would like to use an example dataset, this is labeled `test_data.zip` and can be found within the `Metaboverse` app folder.

### Getting Help
- If you have questions, requests, or bugs to report, please use the [Metaboverse issues forum](https://github.com/Metaboverse/Metaboverse/issues). - Please clearly describe the problem, what you have tried, as well as screenshots of any error information.     
- Generally, for any errors occurring during network building, a file named `metaboverse_session.log` will be output to your specified Output folder. If you receive this file, please upload it to your GitHub Issue. This will output a lot of information, but you can try self-diagnosing by seeing if there is anything in the last ~10-15 lines of this file that might hint at the issue. Otherwise, we are happy to help diagnose the problem!    
- It is also often helpful for us to click on the `View` menu tab, click `Toggle Developer Tools`, click the `Console` tab of the window that opens, and take a screenshot of the output in this panel.

### Feedback
- Have any feedback? Let us know [here](https://forms.gle/4z51DMnagWRvKhc38).
- We also have a discussion forum [here](https://github.com/Metaboverse/Metaboverse/discussions).

### Trying out Metaboverse
With each release archive or Metaboverse, a `test_data.zip` file is included. Unzip this file and read the `README.txt` file for more information on this example dataset.    
