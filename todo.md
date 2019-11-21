# Planned Features (in order of perceived priority)

1. Drop-down menu during visualization to shift between different pathways
2. Plot kNNs of metabolite of interest
3. Generate graph from only those metabolites with known values
4. BRENDA database incorporation
5. Timecourse/flux capabilities
6. Node links to other pathways
7. Show pathway illustration in window to give people classical view
8. Show p-value, etc information along with metabolite of interest



# Instructions for to-do list
0. Save as markdown in repo
- if someone starts working on task, just put name in parenthesis next to task
- If ready for review, turn dash into X
- Once reviewed and merged, strikethrough
- FORK THE REPO

1. Exec file
- X launch server (jordan: see metaboverse/visualize/make.sh)
- X open home page
	(JORDAN)
	- X -> see run.sh in home dir for example, launched home page (application/source/index.html)
- outputs: none

2. Metaboverse home page (index.html)
- (JORDAN)
 - -> Build outline, needs work with actually capturing user dragged folder and file
 - -> Would be nice if organization of page would stay consistent if user drags and resizes window

- Get user output location
	- Create output location if not available
	- Create output /tmp folder

- Get info on new curation or use old
	- Button option to import old
	- If not old,
		- run organisms_list.py
			- Contact Reactome.db
			- Parse organisms
			- Output to tmp folder
		- Read in list from tmp folder
		- Display for user the options	  
	- Output to output/variables.json
		- Input output location
		- If using old, database = file name, if not database = False
		- If new, database = False, organism = name
	- Run exec_utils.py
		- Get system info
		- Append to output/tmp/variables.json
		- Or write to output/tmp/system.json

- Launch curate.py
	- Show loading bar

- After curate.py, run Metaboverse options page (vars.html)
- Have links for docs, etc, other external info for user


(JORDAN)
3. If new curation, run curate.py
- Read in output location and organism name
	- Interpret organism ID
- Get files from Reactome
- Parse and curate
- Output organism database json output/
	- Add date curated to json file or append to end of file name

4. Load user variables (vars.html)
- Get user files
- Populate lists from database chosen
- etc
- Append variables to output/tmp/variables.json
- Run 5a or 5b based on user input (motif analysis, or 1 path viz)

5a. Motif analysis
(JORDAN)
Python:
- For each pathway in organism in database
	- Curate pathway JSON
	- Write to output/tmp or output/tmp/tmp1-n
	- Run motif analysis (py or js?)
	- Get output info, append to table in output/tmp
	- Append table location to output/tmp/variables.json
	- Erase JSON file in output/tmp
	- Enable parallel processing

(YOUJIA)
HTML:
- Launch new browser window to display table (table.html)
	- Populate table
	- Display stats
	- Hyperlink pathway names
		- if user clicks, run step 5b
		- Option to come back to this table
	- Button for user to save table
		- Have user provide file name
		- Add output path + filename + .txt (tab delimited)
	- Button to return to vars.html

5b. Viz 1 Path
(JORDAN)
Python:
- Curate pathway JSON
- Write to output/tmp/path.json

HTML:
- Visualize network
	- Read in JSON file from output/tmp/path.json
	- Draw graph (JORDAN)
	- Drawing features (JORDAN)
		- Buttons to toggle names
		- Input boxes for modifying graph settings
		- Display colorbar
	- Options to save current user view
	- Other user button options
- If user exits page, erase JSON file for pathway
- If user was running motif, back button to go to motif page
- Button to return to vars.html

6. Close browser
- Write to output citation info, user variables, etc
- Clear tmp folder
- Exit

7. Other
- X Organize repo for this restructuring to help team effort (JORDAN)
- Timecourse (JORDAN)
	- Viz
		- https://bl.ocks.org/steveharoz/8c3e2524079a8c440df60c1ab72b5d03
		- http://bl.ocks.org/pranitar/01305d9ad0eba73dbf80
		- https://bl.ocks.org/jrladd/c76799aa63efd7176bd9006f403e854d
		- https://bl.ocks.org/mbostock/6452972
		- https://bl.ocks.org/officeofjane/47d2b0bfeecfcb41d2212d06d095c763
		- Colorbar
			- Need to output in graph info range, color scale, etc
			- https://bl.ocks.org/duspviz-mit/9b6dce37101c30ab80d0bf378fe5e583
			- https://bl.ocks.org/starcalibre/6cccfa843ed254aa0a0d
- Compile omics data in preprocessing, then normalize/prep (JORDAN)
	- Implement DESeq2 wrapper
	- Automated processing of proteomics/metabolomics data? (Talk to Ahmad and Alex)
- X Integrate Travis CI, add badges
- X Implement docs pages
- On home page, have link for its my first time and link to walkthrough video or choose-your-own-adventure thingy
- Make metabolite name mapper
	- Way to standardize?
- Analyte search
	- Plug in analyte and return list of all pathways with that analyte
	- How to handle cross pathway information?
	- Be able to combine motif search with analytes focus search
