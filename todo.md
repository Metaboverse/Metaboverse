# To Do

## v1
	- two condition
	- basic paired time course
	- explain considerations
	- have proof of principle study in well-described system that shows this is reliable
	- inputs need appropriate statistics for this version
	- pathway motif (highest vs lowest or dead-end node comparisons)
	- collapse time-course motif to assess trends with penalization for large spread
	- co-IP shortest paths pseudo-network?
## v2
	- flux
	- advanced time course
	- advanced stats

### COMPLETED:
- Wrap text for reactions, notes
- Current selection wraps and pushes below
- Toggle values doesn't update correctly
- Icons linking to docs
- Can't do 2 kNNs
- Some kNNs will just show parent node (looks like it fixed itself with some previous edit)
- Random arrows with no arm for some links in some pathways
	- Always protein component
	- Looks like some, but not all have expression layering
	- Solution: arose from edges where source and target were identical
- Remove pathways with just one reaction?
	- Solution: no
- Consider super pathways to be those that have more than 200 rxns?
	- Then gather all pathways that fall within those pathways
	- Check that all reactions of pathway must be reactions that are also in the metabolism pathway?
	- Solution: done
- No gene components appear to be showing up
- Get metabolite components to show
- Genes not getting pulled for reactions -- need to add their ID to the reaction database
- Label genes, proteins, metabolites for selective sorting
- Adding gene components makes the scramble and spacing of nodes worse
	- Solution: Gave option to toggle genes from graph
- Write to output citation info, user variables, etc
- Clear tmp files
- Some mismatched arrow and arm colors (because it has two roles?)
	- Not seeing this anymore, so listing as solved
- Protein component edges now showing up
	- Included sub-type in viz curation side


### CRITICAL:
1. Some nodes currently showing up twice with recent curation changes
	- See glycolysis/GCK1:GKRP complex
		- two separate species IDs for complex component and reactant
	- Also, this reaction in kNN seems to be missing modifier or input?
2. It looks like genes may be broadcasting expression to proteins, but I think I gave it protein
	- Only broadcast to protein nodes that don't have values, but if values, leave it. For broadcasted proteins, outline node as dotted border
	- Will likely need to do so in python for colorbar
	- Export colorbar info to JS
	- Need to do stats coloring too
3. issue where some reaction nodes are called gene_name gene or gene_name(n-n+1)
	- This appears to be a non-issue in mapping as mapping to ID fine
4. Some kNN only bring up one node
5. Some genes are not finding their edges (see Attenuation)
	- In case of HSBP1, two HSBP1 nodes for gene components with two ensembl IDs
		- "ENSG00000230989"
		- "ENSG00000106211" -> this should map to HSPB1
	- HSPH1 does not have any links
6. Labeling "Active mTORC1 complex" as metabolite component in Acetylcholine regulates insulin secretion
	- likely occuring because Active mTORC1 has a CHEBI ID
7. Cast gene expression to proteins if no proteomics data
8. Validate viz for new curation paradigm
	- Test some components to make sure they are curating correctly
9. Entity search
	- kNN
	- Component pathways
		- Drop down with each pathway and target larger node, or all pathways together




### To Do:
1. UI
- Get drop in bars to save and import data
- Output log of user actions
- Wrap text for pop-out boxes

2. Curation
- Run motif search or exploration based on user input
- Give user page of entities that we couldn't map, have them enter valid alternative ID

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

6. Time-course data
- Slider bar
- Integration into network curation
	- https://observablehq.com/@mbostock/the-wealth-health-of-nations
- Ideas:
	- https://bl.ocks.org/steveharoz/8c3e2524079a8c440df60c1ab72b5d03
	- https://bl.ocks.org/jrladd/c76799aa63efd7176bd9006f403e854d
	- https://bl.ocks.org/mbostock/6452972
	- https://bl.ocks.org/officeofjane/47d2b0bfeecfcb41d2212d06d095c763
- Motif search integration (simple)
	- See below

7. Advanced
- Use an interactive dial to drag the motif you want (two next to each other)
	- Allow to specify one is gene, one is metabolite, etc
	- build that into the regular motif searcher too
	- https://stackoverflow.com/questions/2368784/draw-on-html5-canvas-using-a-mouse
- For now, limit motif search to one node and its behavior
	- Map them on the network
	- Look at two adjacent motifs and one goes up. Gradually, one goes down gradually

8. Statistics
- Significance
	- See jActiveModules methodology
	- See M-K trend analysis

9. Node collapse
- For reactions with missing data, collapse into bridged reaction for all others
- During hover, in pop-out or sidebar, show missing reactions in full
- Remove nodes with missing data

10. Other
- Travis CI tests
- Docs
- Info for what pathway a node belongs to
- Select parent pathway (i.e. metabolism) in order to let motif search run
- For MIDAS, co-IP, show paths to connect interactors, prioritize by path length

11. Validation
- MPC flux validation
- MetaboNet datasets
- Time-course data

=================================================================

12. Try again?
- Provide a go back button to get back to the graph they made previously
- Make node of interest bigger than others (optional)
- Cartesian distortion?
	- Seems to only work as is with d3 <= v3

13. Future
- Flux data
	- Slider bar
	- Integration into network curation
	- Motif search integration
		- collapse to rate
		- work in incorporation (max labeling)
	- k-NN integration?
- BRENDA curation
	- Extract
	- Figure out how to parse relevant information
	- Allow user to customize
	- Search reactions in other species for motifs
		- Display if motif found in putative reaction  
- Packaging
	- Make available on cloud (host via AWS)
	- How to store information when running online (no output doc, AWS bucket?)
