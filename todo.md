# To Do

## Plan:

### v1
  - two condition
  - basic paired time course
  - explain considerations
  - have proof of principle study in well-described system that shows this is reliable
  - inputs need appropriate statistics for this version
  - pathway motif (highest vs lowest or dead-end node comparisons)
  - collapse time-course motif to assess trends with penalization for large spread
  - co-IP shortest paths pseudo-network?

### v2
  - flux
  - advanced time course
  - advanced stats

## To do:

### In progress:

- Some nodes currently showing up twice with recent curation changes - See glycolysis/GCK1:GKRP complex - two separate species IDs for complex component and reactant - Also, this reaction in kNN seems to be missing modifier or input? - Need to merge
- Validate viz for new curation paradigm - Test some components to make sure they are curating correctly - Are the mappings of data that matter happening?
- Some modifiers seem to be missing their links and when kNN-click only bring up themselves - These cases seem to be consistent and the only time this happens with the kNN function

- Allow custom on-the-fly creation of black list (nodes not to display)
	- display all nodes in order of degree and have a check box -- remove to change, save all as true or false in table for grabbing
- Custom drag and drops that load current session data
- Time-course
- Node collapse
	- https://bl.ocks.org/mbostock/1062288
	- http://bl.ocks.org/GerHobbelt/3071239
- Re-integrate python calls for curation
- Motif design
- include python dependency installs in download
- speed up large networks? Or is it all in the force layout?
- Go back to last network?
- Sometimes protein complexes are registered as metabolite components
	- likely because they also have a CHEBI ID
	- This is probably fine, just include a disclaimer
- Compartments
	- http://bl.ocks.org/GerHobbelt/3071239


### In queue:

- Run motif search or exploration based on user input
- Give user page of entities that we couldn't map, have them enter valid alternative ID - already in table form, just need to show and give option to fill in and rerun?
- Integrate Youjia's module; Figure out licensing
- Allow integration of motif in-pathway panel for exploration module
- Motif: Integrate motif window for real-time analysis
- Option to split out highly connected entities (so in pathway view you see an individual node for each use)
- Export -- give option to export visualizations, data tables, etc.
- Data normalization and merging? Or restrict acceptable inputs
- Speed up curation (archived versions)
- Pre-processing pipelines
- Timecourse motif search integration (simple) - See below
- Use an interactive dial to drag the motif you want (two next to each other) - Allow to specify one is gene, one is metabolite, etc - build that into the regular motif searcher too
	- https://stackoverflow.com/questions/2368784/draw-on-html5-canvas-using-a-mouse
- For now, limit motif search to one node and its behavior - Map them on the network - Look at two adjacent motifs and one goes up. Gradually, one goes down gradually
- Significance - See jActiveModules methodology - See M-K trend analysis
- For reactions with missing data, collapse into bridged reaction for all others
- During hover, in pop-out or sidebar, show missing reactions in full
- Remove nodes with missing data
- Travis CI tests
- Docs
- Info for what pathway a node belongs to
- Select parent pathway (i.e. metabolism) in order to let motif search run
- For MIDAS, co-IP, show paths to connect interactors, prioritize by path length
- MPC flux validation
- MetaboNet datasets
- Time-course data

## Future implementations (non-priority):

- Try again?
- Provide a go back button to get back to the graph they made previously
- Make node of interest bigger than others (optional)
- Cartesian distortion? - Seems to only work as is with d3 <= v3
- Option to hide proteins, other modifiers
- Option to change shapes of different node types
- Flux data - Slider bar - Integration into network curation - Motif search integration - collapse to rate - work in incorporation (max labeling) - k-NN integration?
- BRENDA curation - Extract - Figure out how to parse relevant information - Allow user to customize - Search reactions in other species for motifs - Display if motif found in putative reaction
- Packaging - Make available on cloud (host via AWS) - How to store information when running online (no output doc, AWS bucket?)
- Colorbar - Toggle between expression and significance
- Wrap text for pop-out boxes - change title to div display




- Z-score and FC for motif search

- Can put catalysts in outputs and inhibitors as inputs for motif search -- their expression would be indicative
	- high catalyst means you are presumably getting more output
	- high inhibitor means you are accumulating input

- To modify pathway membership:
	- If a pseudo-pathway is created, add to dictionary
		- Original_reaction_1 => pseudo_reaction_A
		- Original_reaction_2 => pseudo_reaction_A
		- Original_reaction_3 => pseudo_reaction_B
	- Else
		- Original_reaction_4 => Original_reaction_4

- Keep duplicate entities separate
	- Add compartment for reaction and cast to component nodes
	- Make sure mapping is done correctly to apply to both 
