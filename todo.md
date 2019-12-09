# Instructions for to-do list
- If in progress, bold the text using `<b>text</b>`
- Once reviewed and merged, strikethrough using `<del>text</del>`
- FORK THE REPO

# To Do
0. Packaging
- Make executable on workstation
- Make available on cloud (host via AWS)
- Write to output citation info, user variables, etc
- Clear tmp files

1. Front-to-back-end communicaton
	- Variables python to html/js and vice versa
		- Store as JSON formatted variables dictionary?
		```
		var_dict = {
			"transcriptomics": None,
			"proteomics": None,
			"metabolomics": None,
			"output": None,
			"metadata": None,
			"timecourse": false,
			"flux": false,
			"blacklist": [],
			"collapse_missing_reactions": false,
			"split_duplicate_nodes": false,
			"data_min": -5,
			"data_max": 5
		}		
		```
		- https://www.fullstackpython.com/flask.html
	- How to store information when running online (no output doc, AWS bucket?)

2. Front-end
- General
	- Get drop in bars to save and import data
	- Get pages to scale correctly
	- Get metadata to be saved page to page per session
	- Output log of user actions
	- <i><b>(i)</b></i> icons linking to docs

- Home page
	- If curation database given, skip curation page
	- Allow for input of just the network, or network with dataset already integrated

- Curation page
	- Get to run curation based on user input

- Variables page
	- Run motif search or exploration based on user input

- Motif search page
	- Integrate Youjia's module
	- Figure out licensing
	- Allow integration of motif in-pathway panel for exploration module

- Explore page
	- Right side panel (3)
		- more information: On hover, display node name and populate
		- pathway: On side of page, show KEGG or Reactome pathway image, double-click opens new window for that page
		- motif: Integrate motif window for real-time analysis
		- Allow collapse/hide of these side boxes

	- Optional node display
		- option to hide genes, modifiers
		- option to change shapes of different node types

	- Display color bar
		- Need to output in graph info range, color scale, etc
		- https://bl.ocks.org/duspviz-mit/9b6dce37101c30ab80d0bf378fe5e583
		- https://bl.ocks.org/starcalibre/6cccfa843ed254aa0a0d

- Export
	- Give option to export visualizations, data tables, etc.

3. Back-end
- Analyte exploration
	- Drop down of all analytes available in dataset (that have values) and open a k-NN window

- Pathway exploration
	- Fix issue making it so you can't re-explore a network in same session  

- Reactome curation
	- Fix incomplete pathways, missing names
	- Allow for compartment customization
	- Allow custom creation of black list (nodes not to display)
	- Option to split out highly connected entities (so in pathway view you see an individual node for each use)

- BRENDA curation
	- Extract
	- Figure out how to parse relevant information
	- Allow user to customize
	- Search reactions in other species for motifs
		- Display if motif found in putative reaction  

- Data processing
	- Data normalization and merging

- Entity mapping
	- Figure out how to allow for flexible name mapping to prep for Metaboverse
	- Give user page of entities that we couldn't map, have them enter valid alternative ID

- Time-course data
	- Slider bar
	- Integration into network curation
	- Motif search integration
		- ???
	- k-NN integration
		- ???
	- Ideas:
		- https://bl.ocks.org/steveharoz/8c3e2524079a8c440df60c1ab72b5d03
		- http://bl.ocks.org/pranitar/01305d9ad0eba73dbf80
		- https://bl.ocks.org/jrladd/c76799aa63efd7176bd9006f403e854d
		- https://bl.ocks.org/mbostock/6452972
		- https://bl.ocks.org/officeofjane/47d2b0bfeecfcb41d2212d06d095c763

- Flux data
	- Slider bar
	- Integration into network curation
	- Motif search integration
		- ???
	- k-NN integration
		- ???

- k-NN
	- Allow for recursive neighborhood searches
	- Allow custom expansion
	- Fix issue (in curation?) where

- Significance
	- See jActiveModules methodology
	- See M-K trend analysis

- Reaction collapser
	- For reactions with missing data, collapse into bridged reaction for all others
		- During hover, in pop-out or sidebar, show missing reactions in full
	- Remove nodes with missing data

4. Other
- Enable parallel processing
- Implement pipeline, DESeq2 wrappers, etc.
- Travis CI tests
- Docs

5. Validation
- MPC flux
- MetaboNet datasets
- Time-course data 
