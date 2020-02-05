# To Do

### Completed:
- Wrap text for reactions, notes
- Current selection wraps and pushes below
- Toggle values doesn't update correctly
- Icons linking to docs

### Critical:
- Some nodes currently showing up twice with recent curation changes
- Can't do 2 kNNs
- Some kNNs will just show parent node
- Random arrows with no arm for some links in some pathways
- Some mismatched arrow and arm colors (because it has two roles?)
- Remove pathways with just one reaction?
- Consider super pathways to be those that have more than 200 rxns?
	- Then gather all pathways that fall within those pathways
	- Check that all reactions of pathway must be reactions that are also in the metabolism pathway?
- It looks like genes may be broadcasting expression to proteins, but I think I gave it protein
- No gene components appear to be showing up


0. Packaging
- Make available on cloud (host via AWS)
- Write to output citation info, user variables, etc
- Clear tmp files

1. UI
- How to store information when running online (no output doc, AWS bucket?)
- Get drop in bars to save and import data
- Output log of user actions

2. Curation
- Label genes, proteins, metabolites for selective sorting
- Run motif search or exploration based on user input
- Give user page of entities that we couldn't map, have them enter valid alternative ID
- Update viz for new curation paradigm
- Cast gene expression to proteins if no proteomics data

3. Motif searching
- Integrate Youjia's module; Figure out licensing
- Allow integration of motif in-pathway panel for exploration module

4. Viz
- Show Reactome pathway on hover; double-click opens new window for that page
- Motif: Integrate motif window for real-time analysis
- Allow collapse/hide of these side boxes
- Option to hide genes, modifiers
- Option to change shapes of different node types
- Display color bar
	- When toggling values, label with colorbar
	- Need to output in graph info range, color scale, etc
	- https://bl.ocks.org/duspviz-mit/9b6dce37101c30ab80d0bf378fe5e583
	- https://bl.ocks.org/starcalibre/6cccfa843ed254aa0a0d
- Allow for node-based exploration (list of analytes)
- Allow custom creation of black list (nodes not to display)
- Option to split out highly connected entities (so in pathway view you see an individual node for each use)

5. Processing
- Export -- give option to export visualizations, data tables, etc.
- Data normalization and merging? Or restrict acceptable inputs
- Speed up curation (archived versions)
- Pre-processing pipelines

6. BRENDA curation
- Extract
- Figure out how to parse relevant information
- Allow user to customize
- Search reactions in other species for motifs
	- Display if motif found in putative reaction  

7. Time-course data
- Slider bar
- Integration into network curation

- k-NN integration?
- Ideas:
	- https://bl.ocks.org/steveharoz/8c3e2524079a8c440df60c1ab72b5d03
	- https://bl.ocks.org/jrladd/c76799aa63efd7176bd9006f403e854d
	- https://bl.ocks.org/mbostock/6452972
	- https://bl.ocks.org/officeofjane/47d2b0bfeecfcb41d2212d06d095c763
- Motif search integration
	- See below
- Use an interactive dial to drag the motif you want (two next to each other)
	- Allow to specify one is gene, one is metabolite, etc
	- build that into the regular motif searcher too
	- https://stackoverflow.com/questions/2368784/draw-on-html5-canvas-using-a-mouse
- For now, limit motif search to one node and its behavior
	- Map them on the network
	- Look at two adjacent motifs and one goes up. Gradually, one goes down gradually

8. Flux data
- Slider bar
- Integration into network curation
- Motif search integration
	- collapse to rate
	- work in incorporation (max labeling)
- k-NN integration?

9. Statistics
- Significance
	- See jActiveModules methodology
	- See M-K trend analysis

10. Node collapse
- For reactions with missing data, collapse into bridged reaction for all others
- During hover, in pop-out or sidebar, show missing reactions in full
- Remove nodes with missing data

11. Other
- Travis CI tests
- Docs
- Info for what pathway a node belongs to
- Select parent pathway (i.e. metabolism) in order to let motif search run
- Cartesian distortion?

12. Validation
- MPC flux
- MetaboNet datasets
- Time-course data
