""" Sandbox
"""


""" Curation
Include additional info:
- Broad category (metabolism, cell signaling, etc)
- p-value
- Other pathway membership?
"""





""" Name mapping
- genes: default is ensembl gene ID
- proteins: default is gene ID (backend all caps? But save original too for display)
- metabolites: default is chebi ID?

KEEP FETCH SCRIPTS



FROM COMPLEXES DATABASE
=======================
|- should be able to keep this the same, but do first



FROM REACTIONS DATABASE
=======================

pathway_dict (listOfReactions):
|- pathway_1: [reaction_a, reaction_b, ...]
|- pathway_2: [reaction_c, reaction_a, ...]
|- ...


reaction_dict (listOfReactions):
|- reaction_a
|---- id (reaction_392152)
|---- name (Soluble guanylate cyclase converts GTP to cGMP) (remove compartment)
|---- compartment (compartment_70101)
|---- reversible (bool, or from name)
|---- notes (text)
|---- source (reactome, brenda, etc)
|---- reactants:
|------- id: species_390591
|------- reactome_id: from "<bqbiol:is>"
|---------- if reactome_id in complex_dict:
|------------- add componenets? Just genes?
|------- name: ADP:Calcium Bound Sarcomere Protein Complex [cytosol]
|------- components from "<bqbiol:hasPart>":
|------- geneComponents: [Ensembl IDs, ...]
|------- proteinComponents: [uniprot IDs, ...]
|------- metaboliteComponents: [CHEBI IDs, ...]
|------- otherComponents: [anything after last / that doesn't belong with other]
|---- products:
|------- ...
|---- modifiers:
|------- ...
|- reaction_b
|---- ...
|- ...


id_dict (listOfSpecies):
|- ensembleID
|---- ID
|---- name
|---- expression-values
|------- 0: val,
|------- 1: ...
|---- p-value
|---- degree
|---- etc.
|- uniprotID
|---- ...
|- chebiID
|---- ...
|- otherID
|---- ...
|- ...


name_dict (listOfSpecies):
|- nameA-1: id_1
|- nameA-2: id_1
|- nameA-3: id_1
|- nameB-1: id_2
|- nameB-2: id_2
|- ...
Note that proteomics will work best if using UniProt ID
Note that sequencing will work best if using Ensembl ID
Note that metabolomics will work best if using CHEBI ID


metabolite_mapping (from ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz)
compound_id: synonym_name1,
compound_id: synonym_name2,
etc

"""






""" Collapse or high degree node removal
- Do during run time
- If removing nodes, drop down highest connected nodes
- option to split our high degree nodes
"""







""" Time course
- slider bar
- How to incorporate and visualize?
- How to search motifs over time?

- For gene expression, do something like DESeq2 timecourse analysis? Tractable to other omics?
- Or not -- want a network-centric DE method
"""









""" BRENDA
How to know when its too much?
How to show putative relationships?
"""





""" Other
- Proteins are squares with rounded corners
- RNA is diamond

"""
