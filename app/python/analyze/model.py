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

def output_graph(
        graph,
        output_name,
        pathway_dictionary,
        black_list):
    """Output graph and necessary metadata
    """

    data = json_graph.node_link_data(graph)
    data['pathway_dictionary'] = pathway_dictionary
    data['black_list'] = black_list

    with open(output_name, 'w') as f:
        json.dump(data, f, indent=4) # Parse out as array for javascript

"""Graph building
"""
def build_graph(
        network,
        species_reference,
        name_reference,
        complexes,
        species_id):
    """Build graph
    - Add nodes and edges
    - Map names to objects in the graph for display
    - Calculate degree for each node
    """

    # Initialize graph object
    G = nx.DiGraph()

    for reactome_id in network.keys():

        G = process_reactions(
            graph=G,
            reactome_id=reactome_id,
            network=network,
            species_reference=species_reference,
            name_reference=name_reference,
            complex_reference=complexes,
            species_id=species_id)

    # Layer graph metadata for each node


    return G

def process_reactions(
        graph,
        reactome_id,
        network,
        species_reference,
        name_reference,
        complex_reference,
        species_id):
    """
    """

    # Get reaction name
    reaction_id = network[reactome_id]['id']
    reaction_name = network[reactome_id]['name']
    reaction_rev = network[reactome_id]['reversible']
    reaction_notes = network[reactome_id]['notes']

    reactants = network[reactome_id]['reactants']
    products = network[reactome_id]['products']
    modifiers = network[reactome_id]['modifiers'] # ordered list
    modifier_types = network[reactome_id]['modifier_types'] # ordered list

    # Add reaction node
    graph.add_node(reaction_id)
    graph.nodes()[reaction_id]['id'] = reactome_id
    graph.nodes()[reaction_id]['name'] = reaction_name
    graph.nodes()[reaction_id]['reversible'] = reaction_rev
    graph.nodes()[reaction_id]['notes'] = reaction_notes
    graph.nodes()[reaction_id]['type'] = 'reaction'

    # Add vanilla element nodes and their edges
    for reactant in reactants:

        graph = add_node_edge(
            graph=graph,
            id=reactant,
            name=species_reference[reactant],
            reaction_membership=reaction_id,
            type='reactant',
            reversible=reaction_rev)

    for product in products:

        graph = add_node_edge(
            graph=graph,
            id=product,
            name=species_reference[product],
            reaction_membership=reaction_id,
            type='product',
            reversible=reaction_rev)

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
            reversible='false')

    # Expand complexes for components and gene parts
    # Label as component to allow for optional plotting
    # Specify subtypes -- protein_subunit, gene_component, etc.


    return graph

def add_node_edge(
        graph,
        id, # node id
        name, # display name
        reaction_membership,
        type,
        reversible):
    """Add node and edge information to graph
    """

    graph.add_node(id)
    graph.nodes()[id]['id'] = id
    graph.nodes()[id]['name'] = name
    graph.nodes()[id]['type'] = type
    graph.add_edges_from([
        (id, reaction_membership)])
    graph.edges()[(id, reaction_membership)]['type'] = type

    # Add edge direction if reversible
    if reversible == 'true':
        graph.add_edges_from([
            (reaction_membership, id)])
        graph.edges()[(reaction_membership, id)]['type'] = type

    return graph

def check_complexes(
        graph,
        id,
        complex_node,
        complex_reference,
        name_reference,
        species_reference,
        ensembl_reference):
    """Check if species being added is in complex dictionary
    - If record exists, add nodes and edges for the new relationship.
    - If record contains a UniProt ID, cross reference with Ensembl database
    - If complex, label true; else label as false
    """

    if id in complex_reference.keys():
        graph.nodes()[id]['complex'] = 'true'

        participants = complex_reference[id]['participants']
        for p in participants.keys():

            for x in participants[p]:

                if p.lower() == 'chebi':
                    name = 'CHEBI:' + x
                else:
                    name = x

                id = name_reference[name]

                graph = add_node_edge(
                    graph=graph,
                    id=id,
                    name=species_reference[id],
                    reaction_membership=complex_node,
                    type='complex_component',
                    reversible='false')

                if p.lower() == 'uniprot':

                    graph = add_node_edge(
                        graph=graph,
                        id=,
                        name=species_reference[],
                        reaction_membership=id,
                        type='gene_component',
                        reversible='false')

                        ### NEED TO GET NODE NAMES AND DISPLAY NAMES RIGHT FOR EACH




                if p.lower() == 'uniprot':
                    name = x






    else:
        graph.nodes()[id]['complex'] = 'false'

    return graph


network['complex_dictionary']['species_1006173']









"""Data overlay
- Map repo id to species_id
- If a node is a complex, take average of neighbors that are not
"""







def __main__(
        network,
        data,
        stats,
        species_id,
        output_file,
        black_list):
    """Generate graph object for visualization
    - Place black_list as key in graph object for later parsing
        - Will allow for on-the-fly removal of nodes
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
        sep='\t')

    stats = pd.read_csv(
        '/Users/jordan/Desktop/metaboverse/app/python/analyze/test/cat_stats.txt',
        sep='\t')

    species_id = 'HSA'
    output_file = 'HSA_global_reactions.json'
    #############################

    # Generate output file name
    graph_name = name_graph(
        output_file=output_file,
        species_id=species_id)

    # Generate graph
    # Name mapping
    G = build_graph(
        network=network['reaction_database'],
        species_reference=network['species_database'],
        name_reference=network['name_database'],
        complexes=network['complex_dictionary'],
        species_id=species_id)

    # Overlay data and stats, calculate heatmap values for p-value
    # and expression value
    G, names_dictionary = map_attributes(
        graph=G,
        data=data,
        stats=stats)

    # Export graph, pathway membership, pathway degree, black_list, other refs
    output_graph(
        graph=G,
        output_name=graph_name,
        pathway_dictionary=network['pathway_database'],
        black_list=black_list)
