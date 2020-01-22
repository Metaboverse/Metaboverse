""" Sandbox
"""


""" Curation
Include additional info:
- Broad category (metabolism, cell signaling, etc)
- p-value
- Other pathway membership?

Graph_database
nodes:
    species_id: {
        display_name: ___,
        expression_value: [],
        p_value: [],
        rgba_value: [ [], ],
        rjba_js: [ [], ],
        type: ___,
    }


links:
    {target: species_id,
    source: species_id,
    type: ___
    }

pathways:

reactions:

names:

super_pathways: --> do in interactive mode
    - display all pathway for selection and in parentheses display number of nodes (sort in order of node #)
    - user selects pathway and super pathway for downstream analysis is assembled



1) convert user data ids to values where available
    - convert gene expression data to ensembl ID
    - convert protein expresison data to uniprot ID
    - cross reference metabolite data to CHEBI dictionary
2) return a table for values with missing ids
3) generate graph
4) finalize formatting for viz
5) export

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
Note that proteomics will work best if using UniProt ID
Note that sequencing will work best if using Ensembl ID
Note that metabolomics will work best if using CHEBI ID

metabolite_mapping:
compound_id: synonym_name1,
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
- Display color bar in left side of page
- node search
"""
