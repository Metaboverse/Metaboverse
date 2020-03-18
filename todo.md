# To do:

### In progress:
- Some nodes act as both catalysts and inhibitors, but only get labeled as one in the metadata
    - label as:
        - catalyst_species_170724
        - inhibitor_species_170724
        - catalyst_170724
        - inhibitor_170724
    - How to ensure values are mapped properly?

- Compartments
	- http://bl.ocks.org/GerHobbelt/3071239
	- Will have to go back and add compartment info
	- On hover, display compartment name

- Z-score and FC for motif search

- Some modifiers seem to be missing their links and when kNN-click only bring up themselves - These cases seem to be consistent and the only time this happens with the kNN function

- Custom drag and drops that load current session data

- Time-course

- Re-integrate python calls for curation

- Motif design

- include python dependency installs in download

- speed up large networks? Or is it all in the force layout?

- Go back to last network?

- Sometimes protein complexes are registered as metabolite components
	- likely because they also have a CHEBI ID
	- This is probably fine, just include a disclaimer

- Broadcasting
	- take absolute minimum? Max? Avg? Let user decide?

### In queue:
- Give user page of entities that we couldn't map, have them enter valid alternative ID - already in table form, just need to show and give option to fill in and rerun?
- Integrate Youjia's module; Figure out licensing
- Allow integration of motif in-pathway panel for exploration module
- Motif: Integrate motif window for real-time analysis
- Option to split out highly connected entities (so in pathway view you see an individual node for each use)
- Export -- give option to export visualizations, data tables, etc.
- Data normalization and merging? Or restrict acceptable inputs
- Speed up curation (archived versions)
- Pre-processing pipelines
- Time-course motif search integration (simple) - See below
- Use an interactive dial to drag the motif you want (two next to each other) - Allow to specify one is gene, one is metabolite, etc - build that into the regular motif searcher too
	- https://stackoverflow.com/questions/2368784/draw-on-html5-canvas-using-a-mouse
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
- Make sure python back-end is Windows compatible 


## Future implementations (non-priority):
- Provide a go back button to get back to the graph they made previously
- Make node of interest bigger than others (optional)
- Cartesian distortion? - Seems to only work as is with d3 <= v3
- Option to hide proteins, other modifiers
- Option to change shapes of different node types
- Flux data - Slider bar - Integration into network curation - Motif search integration - collapse to rate - work in incorporation (max labeling) - k-NN integration?
- BRENDA curation - Extract - Figure out how to parse relevant information - Allow user to customize - Search reactions in other species for motifs - Display if motif found in putative reaction
- Packaging - Make available on cloud (host via AWS) - How to store information when running online (no output doc, AWS bucket?)
- Wrap text for pop-out boxes - change title to div display
