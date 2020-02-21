"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
    Copyright (C) 2019  Jordan A. Berg
    jordan <dot> berg <at> biochem <dot> utah <dot> edu

    This program is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License along with
    this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""Import dependencies
"""
import os
import pandas as pd
import numpy as np
from math import sqrt
from ast import literal_eval
import json
import pickle
import networkx as nx
from networkx.readwrite import json_graph
import matplotlib
import matplotlib.pyplot as plt
cmap = matplotlib.cm.get_cmap('seismic')

"""Graph utils
"""
def name_graph(
        output_file,
        species_id):
    """Name graph
    """

    if output_file[-5:].lower() == '.json':
        graph_name = output_file
    else:
        graph_name = output_file + species_id + '_global_reactions.json'

    return graph_name

"""Graph building
"""
def build_graph(
        network,
        species_reference,
        name_reference,
        protein_reference,
        complexes,
        species_id,
        reverse_genes):
    """Build graph
    - Add nodes and edges
    - Map names to objects in the graph for display
    - Calculate degree for each node
    """

    # Initialize graph object
    G = nx.DiGraph()

    for reactome_id in network.keys():

        G, network = process_reactions(
            graph=G,
            reactome_id=reactome_id,
            network=network,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            complex_reference=complexes,
            species_id=species_id,
            reverse_genes=reverse_genes)

    return G, network

def process_reactions(
        graph,
        reactome_id,
        network,
        species_reference,
        name_reference,
        protein_reference,
        complex_reference,
        species_id,
        reverse_genes):
    """
    """
    new_genes = []

    # Get reaction name
    reaction_id = network[reactome_id]['id']
    reaction_name = network[reactome_id]['name']
    reaction_rev = network[reactome_id]['reversible']
    reaction_notes = network[reactome_id]['notes']

    reactants = network[reactome_id]['reactants']
    products = network[reactome_id]['products']
    modifiers = network[reactome_id]['modifiers'] # ordered list

    # Add reaction node
    graph.add_node(reaction_id)
    graph.nodes()[reaction_id]['id'] = reactome_id
    graph.nodes()[reaction_id]['name'] = reaction_name
    graph.nodes()[reaction_id]['reversible'] = reaction_rev
    graph.nodes()[reaction_id]['notes'] = reaction_notes
    graph.nodes()[reaction_id]['type'] = 'reaction'
    graph.nodes()[reaction_id]['sub_type'] = 'reaction'

    # Add vanilla element nodes and their edges
    for reactant in reactants:

        graph = add_node_edge(
            graph=graph,
            id=reactant,
            name=species_reference[reactant],
            reaction_membership=reaction_id,
            type='reactant',
            sub_type='reactant',
            reversible=reaction_rev,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference)

        graph, gene_components = check_complexes(
            graph=graph,
            id=reactant,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            reverse_genes=reverse_genes)
        for x in gene_components:
            new_genes.append(x)

    for product in products:

        graph = add_node_edge(
            graph=graph,
            id=product,
            name=species_reference[product],
            reaction_membership=reaction_id,
            type='product',
            sub_type='product',
            reversible=reaction_rev,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference)

        graph, gene_components = check_complexes(
            graph=graph,
            id=product,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            reverse_genes=reverse_genes)
        for x in gene_components:
            new_genes.append(x)

    for modifier in modifiers:

        # Extract modifier type
        # Formatted as a list of lists with first index of each sub list being
        # the species ID and the second index being the modifier type
        # ex: [['species_0', 'inhibitor'], ['species_1', 'catalyst']]
        # Labeling the edge should allow for differentiation between the same
        # modifier node acting as a catalyst or inhibitor
        id = modifier[0]
        type = modifier[1]

        graph = add_node_edge(
            graph=graph,
            id=id,
            name=species_reference[id],
            reaction_membership=reaction_id,
            type=type,
            sub_type='modifier',
            reversible='false',
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference)

        graph, gene_components = check_complexes(
            graph=graph,
            id=id,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            reverse_genes=reverse_genes)
        for x in gene_components:
            new_genes.append(x)

    network[reactome_id]['genes'] = new_genes

    return graph, network

def add_node_edge(
        graph,
        id, # node id
        name, # display name
        reaction_membership,
        type,
        sub_type,
        reversible,
        complex_reference,
        species_reference,
        name_reference,
        protein_reference):
    """Add node and edge information to graph
    """

    graph.add_node(id)
    graph.nodes()[id]['id'] = id
    graph.nodes()[id]['name'] = name
    graph.nodes()[id]['type'] = type
    graph.nodes()[id]['sub_type'] = sub_type

    if type == 'reactant':
        graph.add_edges_from([
            (id, reaction_membership)])
        graph.edges()[(id, reaction_membership)]['type'] = type
        graph.edges()[(id, reaction_membership)]['sub_type'] = sub_type

        if reversible == 'true':
            graph.add_edges_from([
                (reaction_membership, id)])
            graph.edges()[(reaction_membership, id)]['type'] = type
            graph.edges()[(reaction_membership, id)]['sub_type'] = sub_type

    elif type == 'product':
        graph.add_edges_from([
            (reaction_membership, id)])
        graph.edges()[(reaction_membership, id)]['type'] = type
        graph.edges()[(reaction_membership, id)]['sub_type'] = sub_type

        if reversible == 'true':
            graph.add_edges_from([
                (id, reaction_membership)])
            graph.edges()[(id, reaction_membership)]['type'] = type
            graph.edges()[(id, reaction_membership)]['sub_type'] = sub_type

    else:
        graph.add_edges_from([
            (id, reaction_membership)])
        graph.edges()[(id, reaction_membership)]['type'] = type
        graph.edges()[(id, reaction_membership)]['sub_type'] = sub_type

    return graph

def check_complexes(
        graph,
        id,
        complex_reference,
        species_reference,
        name_reference,
        protein_reference,
        reverse_genes):
    """Check if species being added is in complex dictionary
    - If record exists, add nodes and edges for the new relationship.
    - If record contains a UniProt ID, cross reference with Ensembl database
    - If complex, label true; else label as false
    """

    add_components = []

    if id in complex_reference.keys():
        graph.nodes()[id]['complex'] = 'true'

        participants = complex_reference[id]['participants']
        for p in participants.keys():

            for x in participants[p]:

                if p.lower() == 'chebi':
                    name = 'CHEBI:' + x
                    sub_type = 'metabolite_component'
                else:
                    name = x

                    if p.lower() == 'uniprot':

                        sub_type = 'protein_component'
                    elif p.lower() == 'mirna':
                        sub_type = 'mirna_component'
                    elif p.lower() == 'ensembl':
                        sub_type = 'gene_component'
                    else:
                        sub_type = 'other'

                try:
                    component_id = name_reference[name]
                    #add_components.append(component_id)

                    graph = add_node_edge(
                        graph=graph,
                        id=component_id,
                        name=species_reference[component_id],
                        reaction_membership=id,
                        type='complex_component',
                        sub_type=sub_type,
                        reversible='false',
                        complex_reference=complex_reference,
                        species_reference=species_reference,
                        name_reference=name_reference,
                        protein_reference=protein_reference)

                    if p.lower() == 'uniprot':

                        # Get protein's corresponding gene ID
                        gene = protein_reference[name]
                        add_components.append(gene)

                        graph = add_node_edge(
                            graph=graph,
                            id=gene,
                            name=reverse_genes[gene],
                            reaction_membership=component_id,
                            type='gene_component',
                            sub_type='gene',
                            reversible='false',
                            complex_reference=complex_reference,
                            species_reference=species_reference,
                            name_reference=name_reference,
                            protein_reference=protein_reference)

                except:
                    pass
                    #print('Could not retrieve components for', name)

    else:
        graph.nodes()[id]['complex'] = 'false'

    return graph, add_components

def uniprot_ensembl_reference(
        uniprot_reference,
        ensembl_reference):
    """Build cross-referencing dictionary to convert uniprot ID to
    corresponding Ensembl ID
    """

    new_dict = {}

    for k, v in uniprot_reference.items():
        try:
            new_dict[v] = ensembl_reference[k]
        except:
            pass

    for k, v in ensembl_reference.items():
        try:
            new_dict[uniprot_reference[k]] = v
        except:
            pass

    return new_dict

def map_attributes(
        graph,
        data,
        stats,
        name_reference,
        degree_dictionary):
    """Data overlay
    - Map repo id to species_id
    - If a node is a complex, take average of neighbors that are not
    To do:
    - Currently, many metabolites that should map are not found in name
    database
    """

    n = len(data.columns.tolist())
    reaction_color = (0.75, 0.75, 0.75, 1)
    missing_color = (1, 1, 1, 1)

    # Re-index data and stats
    data_renamed = data.copy()
    data_renamed.index = data.index.map(network['name_database'])
    data_renamed = data_renamed.loc[data_renamed.index.dropna()]
    data_max = abs(data_renamed).max().max()

    stats_renamed = stats.copy()
    stats_renamed.index = stats.index.map(network['name_database'])
    stats_renamed = stats_renamed.loc[stats_renamed.index.dropna()]
    stats_max = abs(stats_renamed).max().max()

    data_dict = {}
    for index, row in data_renamed.iterrows():
        data_dict[index] = list(row)

    stats_dict = {}
    for index, row in stats_renamed.iterrows():
        stats_dict[index] = list(row)

    for x in list(graph.nodes()):

        current_id = graph.nodes()[x]['id']

        # Add degree
        graph.nodes()[x]['degree'] = degree_dictionary[current_id]

        if current_id in set(data_dict.keys()):

            graph.nodes()[x]['values'] = data_dict[current_id]
            graph.nodes()[x]['values_rgba'] = extract_value(
                value_array=data_dict[current_id],
                max_value=data_max)
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['values_rgba'])

            graph.nodes()[x]['stats'] = stats_dict[current_id]
            graph.nodes()[x]['stats_rgba'] = extract_value(
                value_array=stats_dict[current_id],
                max_value=stats_max)
            graph.nodes()[x]['stats_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['stats_rgba'])

        else:
            if graph.nodes()[x]['type'] == 'reaction':
                colors = [reaction_color for x in range(n)]

            else:
                colors = [missing_color for x in range(n)]

            graph.nodes()[x]['values'] = [None for x in range(n)]
            graph.nodes()[x]['values_rgba'] = colors
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=colors)

            graph.nodes()[x]['stats'] = [None for x in range(n)]
            graph.nodes()[x]['stats_rgba'] = colors
            graph.nodes()[x]['stats_js'] = convert_rgba(
                rgba_tuples=colors)

    return graph

def extract_value(
        value_array,
        max_value):
    """Extract expression value
    """

    rgba = []
    for x in value_array:

        position = (x + max_value) / (2 * max_value)
        rgba_tuple = cmap(position)
        rgba.append(rgba_tuple)

    return rgba

def convert_rgba(
        rgba_tuples):
    """Convert python RGBA tuple to web-friendly tuple for later viz
    """

    js = []
    for x in rgba_tuples:

        rgba_list = list(x)
        rgba_new = []
        for x in rgba_list[:3]:
            rgba_new.append(int(x * 255))

        rgba_new.append(rgba_list[3])

        js.append(tuple(rgba_new))

    return js

def output_graph(
        graph,
        output_name,
        pathway_dictionary,
        super_pathways,
        reaction_dictionary,
        black_list):
    """Output graph and necessary metadata
    """

    data = json_graph.node_link_data(graph)
    data['pathway_dictionary'] = pathway_dictionary
    data['super_pathways'] = super_pathways
    data['reaction_dictionary'] = reaction_dictionary
    data['black_list'] = black_list

    with open(output_name, 'w') as f:
        json.dump(data, f, indent=4) # Parse out as array for javascript

def compile_pathway_degree(
        pathways):
    """Compile database of large pathways
    """

    super_pathways = {}

    for key in list(pathways.keys()):

        if len(pathways[key]["reactions"]) > 200:
            super_pathways[key] = pathways[key]

    return super_pathways

def compile_node_degrees(
        graph):
    """Retrieve degree metrics for each node
    """

    d = {}

    deg_dict = G.degree
    for k,v in deg_dict:

        d[k] = v

    return d

def __main__(
        network,
        data,
        stats,
        species_id,
        output_file,
        black_list=[]):
    """Generate graph object for visualization
    - Place black_list as key in graph object for later parsing
        - Will allow for on-the-fly removal of nodes
    To do:
    - Map average component expression to complex nodes
    - Determine product reactant type -- metabolite, protein, etc
    """

    #############################
    def read_network(
            network_url):
        """Read in network from previous curation module
        - was provided as a URL to the file and saved to args_dict['network']  in "curate" sub-module
        """
        import pickle
        with open(network_url, 'rb') as network_file:
            network = pickle.load(network_file)

        return network

    network = read_network(
        network_url='/Users/jordan/Desktop/HSA_metaboverse_db.pickle')

    data = pd.read_csv(
        '/Users/jordan/Desktop/metaboverse/app/python/analyze/test/cat_data.txt',
        sep='\t',
        index_col=0)

    stats = pd.read_csv(
        '/Users/jordan/Desktop/metaboverse/app/python/analyze/test/cat_stats.txt',
        sep='\t',
        index_col=0)

    output_file = '/Users/jordan/Desktop/HSA_global_reactions.json'

    species_id = 'HSA'
    black_list=[]

    """Start of pancancer-necessary code
    To run, place all graph-making code within the for loop below
    """
    """
    cancers = [
    ['breast', 'brca'],
    ['colon', 'coad'],
    ['colon','read'],
    ['bladder','blca'],
    ['cervix','cesc'],
    ['kidney','kirc'],
    ['kidney','kirp'],
    ['liver','lihc'],
    ['lung','luad'],
    ['lung','lusc'],
    ['prostate','prad'],
    ['stomach','stad'],
    ['thyroid','thca'],
    ['uterus','ucec']]

    for x in [x[1] for x in cancers]:
        output_file = '/Users/jordan/Desktop/HSA_global_reactions_' + x + '.json'
        print(output_file)
        data = stats = pd.read_csv(
            '/Users/jordan/Desktop/cancer_network_topology/tables/' + x + '_z_table.txt',
            sep='\t',
            index_col=0)
    """


    #############################

    print('Preparing metadata...')
    # Generate output file name
    graph_name = name_graph(
        output_file=output_file,
        species_id=species_id)

    # Prepare uniprot to ensembl name mapper
    protein_dictionary = uniprot_ensembl_reference(
        uniprot_reference=network['uniprot_synonyms'],
        ensembl_reference=network['ensembl_synonyms'])
    reverse_genes = {v:k for k,v in network['ensembl_synonyms'].items()}

    # Generate graph
    # Name mapping
    print('Building network...')
    G, network['reaction_database'] = build_graph(
        network=network['reaction_database'],
        species_reference=network['species_database'],
        name_reference=network['name_database'],
        protein_reference=protein_dictionary,
        complexes=network['complex_dictionary'],
        species_id=species_id,
        reverse_genes=reverse_genes)

    # For gene components, add section to reaction database for gene_components and list
    # Pull those in with everything else in JS

    # Overlay data and stats, calculate heatmap values for p-value
    # and expression value
    print('Mapping user data...')
    degree_dictionary = compile_node_degrees(
        graph=G)
    G = map_attributes(
        graph=G,
        data=data,
        stats=stats,
        name_reference=network['name_database'],
        degree_dictionary=degree_dictionary)

    # Generate list of super pathways (those with more than 200 reactions)
    super_pathways = compile_pathway_degree(
        pathways=network['pathway_database'])
    ###
    # For pancancer
    """
    metabolism = network['pathway_database']['R-HSA-1430728']
    network['pathway_database'] = {}
    network['pathway_database']['R-HSA-1430728'] = metabolism

    with open('/Users/jordan/Desktop/cancer_network_topology/network/cancer_reactions.json') as f:
        file = json.load(f)
        for k,v in file.items():

            network['pathway_database'][k] = {}
            network['pathway_database'][k]['id'] = k
            network['pathway_database'][k]['name'] = k
            network['pathway_database'][k]['reactions'] = v
    """
    ###

    # Export graph, pathway membership, pathway degree, black_list, other refs
    print('Exporting graph...')
    output_graph(
        graph=G,
        output_name=graph_name,
        pathway_dictionary=network['pathway_database'],
        super_pathways=super_pathways,
        reaction_dictionary=network['reaction_database'],
        black_list=black_list)
    print('Graphing complete.')
