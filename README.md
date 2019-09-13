# Metaboverse
Biological and Metabolic Networks Analyzer

`WARNING: This repository is currently under development. When the software is complete and stable, the beta tag will be removed.`

## What does Metaboverse do?
### tl;dr
- Curates metabolic model (currently only for human)
- Maps networks or sub-networks interactively
- Overlays selected network with transcriptomics, proteomics, and/or metabolomics data
- Performs trend analysis across networks

### For those more details focused
```
Metabolism is a complex network of chemical and enzymatic reactions; however, the past century has seen
a largely reductionist approach to tackling this system. While for decades this was necessary due to
technological limitations, we are now living in the computer age, where entire cells and populations of
cells can be modeled and explored. While scientific fields, such as RNA biology and metabolism, have
experienced massive strides in recent decades with advent of RNA-seq and metabolomics, our ability to
contextualize these massive amounts of data continues to lag. This is problematic as these experiments
are often expensive and time-consuming to produce, yet we only use a fraction of the total data made
available by the experiment. In order to address these limitations, we introduce Metaboverse, a
computational analysis framework for contextualizing -omics datasets within customizable metabolic
network representations. This framework will allow for static and dynamic time-course exploration of
these datasets, and importantly will help contextualize the role of low-expressed analytes within a
network that are normally ignored due to their classification as not “differentially expressed” by most
algorithms. This tool will revolutionize our ability to more holistically understand temporal metabolic
shifts and gene-metabolite inter-cooperativity, as well as ensure we are obtaining the maximum
information from these datasets as possible.
```

## Getting started

### Installation
```
pip install metaboverse
```

### Getting Started

#### Network Model Curation
- Requires the reactome database files (will download most current automatically)
- Run as:
```
$ metaboverse curate --output /path/to/output
```

#### Data Curation
- Requires properly formatted and normalized -omics data
  - tab-separated library-normalized data table
  - samples metadata table
- For example, with transcriptomics data, run as:
```
$ metaboverse preprocess -d dataset.tsv -t transcriptomics
```

#### Network Analysis


## To Do:
- Network curation:
1. Add log file prints
2. Decide if both recon and hmdb are required inputs each time
3. Directionality?
4. Allostery via Brenda?
5. Build in network customization at correct spot (curation or analysis?)
6. To network model for analysis, add:
  - Allostery dictionary with name and interaction type
  - If reaction involves protein or complex, mapping for genes
    - i.e. {'complexA': ['gene1', 'gene2', 'gene3']}
    - then take the average for complex expression, with ability to draw out each individual gene
    - maybe add this data to actual data as a re-naming for mapping?
7. Re-write enhancement in C, Cython, or Julia? (too slow)
8. Remove intermediate files from curation to just pass model_dict var between modules






#### node_type 1: Reactions
reaction_name: {
  reactants:,
  products:,
  enzymes:,
  directionality:,
  compartments:,
  processes:,
}


#### node_type 2: Metabolite



#### edge_type: metabolite-reaction relations
- remove high-degree nodes (over 50, i.e. proton to reduce tangle)
-


#### analysis
- for de novo analysis, cycle through all processes and generate graphs, ping if interesting pattern
- allow toggling of compartments and hubs, but by default
  - compartments = false
  - hubs = remove any node with more than 50 relations
- needs to be in networkx pickle object
nodes
  - metabolite
  - reaction
edges
  - reaction-Metabolite
  - metabolite-metabolite
  - reaction-reaction



metaboverse v1.0

- manually curate the central carbon network as proof of principle
  - metabolites
  - complexes (name)
    - components (gene_id, protein_id)
  - relationships
  - allostery
  - highlight breakpoints
    - stats and other analytical tools
    - spatial representation
  - analyze a bunch of metabolomic datasets and use tool to gain better insight
  - noise
  - predictive modeling project
- bioRxiv
- cell systems or better
  - cell metabolism if OA

metaboverse v2.0
- Rehash curator to allow for full model curation with proper format and info
- cycle through each process and
- try to deliver by time of bioRxiv?
