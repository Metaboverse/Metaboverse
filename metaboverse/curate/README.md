# Reactome plan

- curate database based on organism
```
organisms = {
    'Arenicola marina',
    'Bos taurus',
    'Caenorhabditis elegans',
    'Canis familiaris',
    'Cricetulus griseus',
    'Crithidia fasciculata',
    'Danio rerio',
    'Dictyostelium discoideum',
    'Drosophila melanogaster',
    'Escherichia coli',
    'Felis catus',
    'Gallus gallus',
    'Homarus americanus',
    'Homo sapiens',
    'Meleagris gallopavo',
    'Mus musculus',
    'Mycobacterium tuberculosis',
    'Oryctolagus cuniculus',
    'Penicillium chrysogenum',
    'Plasmodium falciparum',
    'Rattus norvegicus',
    'Saccharomyces cerevisiae',
    'Schizosaccharomyces pombe',
    'Sus scrofa',
    'Triticum aestivum',
    'Vigna radiata var. radiata',
    'Xenopus laevis',
    'Xenopus tropicalis'}
```

- read in reactome files for all organisms and parse for organism
- compile the databse
- allostery database and miRNA populate inhibitors and activators
- relations help create edges for network

```
reactome = {

    reaction: {
        id: name
        compartment: []
        processes: []
        directionality: forward, reverse, both
        reactants: []
        products: []
    }

    nodetype1
    'ChEBI2Reactome_PE_All_Levels.txt' --> appears to be most complete for ChEBI
    metabolite: {
        id: name
        compartment: []
        processes: []
        inhibitors: []
        activators: []
        relations: []
    }

    nodetype2
    enzyme: {
        id: name
        compartment: []
        processes: []
        inhibitors: []
        activators: []
        gene_components: []
        protein_components: []
        relations: []
    }
}
```

- pickle
- import into network module for analysis
- have a metabolite name to id mapper, transcriptome, etc




## To Do:
- map metabolites to reactions_db
- map enzymes to reactions_db
  - map genes/proteins to enzymes
