# Instructions for to-do list
- If in progress, bold the text using `<b>text</b>`
- Once reviewed and merged, strikethrough using `<del>text</del>`
- FORK THE REPO

# To Do
0. Packaging
- <del>Make executable on workstation</del>
- Make available on cloud (host via AWS)
- Write to output citation info, user variables, etc
- Clear tmp files

1. Front-to-back-end communicaton
	- <del>Variables python to html/js and vice versa</del>
		- <del>Store as JSON formatted variables dictionary?</del>
	- How to store information when running online (no output doc, AWS bucket?)

2. Front-end
- General
	- Get drop in bars to save and import data
	- <del>Get pages to scale correctly</del>
	- <del>Get metadata to be saved page to page per session</del>
	- Output log of user actions
	- <del><i><b>(i)</b></i> icons linking to docs</del>

- Home page
	- <del>If curation database given, skip curation page</del>
	- <del>Allow for input of just the network, or network with dataset already integrated</del>

- Curation page

- Variables page
	- <del>Run curation with loading bar</del>
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
	- <del>Fix issue making it so you can't re-explore a network in same session</del>  

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
	- <b>Figure out how to allow for flexible name mapping to prep for Metaboverse</b>
	- <b>Give user page of entities that we couldn't map, have them enter valid alternative ID</b>

- Time-course data
	- Slider bar
	- Integration into network curation
	- Motif search integration
		- ???
	- k-NN integration
		- ???
	- Ideas:
		- https://bl.ocks.org/steveharoz/8c3e2524079a8c440df60c1ab72b5d03
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
	- <del>Allow for recursive neighborhood searches</del>
	- Allow custom expansion
	- <del>Fix issue (in curation?) where network data variable is mutated</del>
	- <del>Fix hanging object node</del>

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
- Keep feature labels on if user specifies
- Click on one node and show its neighbors
- Info for what pathway a node belongs to
- kNN search by metabolite name
- Select parent pathway (i.e. metabolism) in order to let motif search run


5. Validation
- MPC flux
- MetaboNet datasets
- Time-course data


6. Other  
- Cartesian distortion?




### Notes:
- currently prioritizing based on number of pathways a reaction motif is found
- write js function to extract metabolic sub network formatted in usual way
- or
