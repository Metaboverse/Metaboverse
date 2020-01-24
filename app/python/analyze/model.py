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

        graph.add_node(reactant)
        graph.nodes()[reactant]['id'] = reactant
        graph.nodes()[reactant]['type'] = 'reactant'
        graph.add_edges_from([
            (reactant, reaction_id)])
        graph.edges()[(reactant, reaction_id)]['type'] = 'reactant'

        # Add edge direction if reversible
        if reaction_rev == 'true':
            graph.add_edges_from([
                (reaction_id, reactant)])
            graph.edges()[(reaction_id, reactant)]['type'] = 'reactant'

    for product in products:

        graph.add_node(product)
        graph.nodes()[product]['id'] = product
        graph.nodes()[product]['type'] = 'product'
        graph.add_edges_from([
            (reaction_id, product)])
        graph.edges()[(reaction_id, product)]['type'] = 'product'

        # Add edge direction if reversible
        if reaction_rev == 'true':
            graph.add_edges_from([
                (product, reaction_id)])
            graph.edges()[(product, reaction_id)]['type'] = 'product'

    for modifier in modifiers:

        graph.add_node(modifier)
        graph.nodes()[modifier]['id'] = modifier
        graph.nodes()[modifier]['type'] = 'modifier'

        # Extract modifier type
        # Labeling the edge should allow for differentiation between the same
        # modifier node acting as a catalyst or inhibitor
        graph.add_edges_from([
            (modifier, reaction_id)])
        graph.edges()[(modifier, reaction_id)]['type'] = 'modifier'



    # Expand complexes for components and gene parts






network['reaction_database'].keys()


network['species_database']['species_54639']




"""Data overlay
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
