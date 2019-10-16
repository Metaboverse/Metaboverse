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

"""Set output directory for graph files
"""
def create_graph_output(
        output):

    graph_directory = output + 'graphs/'

    if not os.path.exists(graph_directory):
        os.makedirs(graph_directory)

    return graph_directory

def build_graph(
        sub_network,
        master_reference,
        complex_reference,
        species_id,
        black_list):

    G = nx.DiGraph()

    # cycle through reactions in process
    for reaction_id in sub_network['reactions'].keys():

        G = process_reactions(
            graph=G,
            reaction_id=reaction_id,
            sub_network=sub_network,
            master_reference=master_reference,
            complex_reference=complex_reference,
            species_id=species_id,
            black_list=black_list)

    return G

def process_reactions(
        graph,
        reaction_id,
        sub_network,
        master_reference,
        complex_reference,
        species_id,
        black_list):

    reaction_name = sub_network['reactions'][reaction_id]['name']

    # Add reaction node
    graph = add_reaction_node(
        graph=graph,
        reaction_id=reaction_id,
        reaction_name=reaction_name)

    # Add reactant nodes and edges
    graph = add_node_edge(
        graph=graph,
        reaction_name=reaction_name,
        sub_network=sub_network['reactions'][reaction_id],
        master_reference=master_reference,
        complex_reference=complex_reference,
        species_id=species_id,
        type='reactants',
        black_list=black_list)

    # Add product nodes and edges
    graph = add_node_edge(
        graph=graph,
        reaction_name=reaction_name,
        sub_network=sub_network['reactions'][reaction_id],
        master_reference=master_reference,
        complex_reference=complex_reference,
        species_id=species_id,
        type='products',
        black_list=black_list)

    # Add modifier nodes and edges
    graph = add_node_edge(
        graph=graph,
        reaction_name=reaction_name,
        sub_network=sub_network['reactions'][reaction_id],
        master_reference=master_reference,
        complex_reference=complex_reference,
        species_id=species_id,
        type='modifiers',
        black_list=black_list)

    return graph

"""Add reaction node config information
"""
def add_reaction_node(
        graph,
        reaction_id,
        reaction_name):

    graph.add_node(reaction_name)
    graph.nodes()[reaction_name]['name'] = reaction_name
    graph.nodes()[reaction_name]['id'] = reaction_id
    graph.nodes()[reaction_name]['type'] = 'reaction'
    graph.nodes()[reaction_name]['stoichiometry'] = None

    return graph

"""Add analyte node config information and relationships
"""
def add_node_edge(
        graph,
        reaction_name, # Use this to link edge back to reaction
        sub_network, # This starts at the reaction level
        master_reference,
        complex_reference,
        species_id,
        type,
        black_list):

    # Will cycle through all the analytes within a given analyte type: i.e., reactant, product, modifier
    for analyte_id in sub_network[type].keys():

        # Get analyte name and check ID
        analyte_name, analyte_id_update = fetch_analyte_info(
            master_reference=master_reference,
            analyte_id=analyte_id,
            species_id=species_id)

        # Add node and edge to reaction
        if analyte_name not in black_list:

            # Add analyte node
            try:
                graph = add_analyte_node(
                    graph=graph,
                    sub_network=sub_network,
                    analyte_name=analyte_name,
                    analyte_id=analyte_id,
                    type=type)
            except:
                graph = add_analyte_node(
                    graph=graph,
                    sub_network=sub_network,
                    analyte_name=analyte_name,
                    analyte_id=analyte_id_update,
                    type=type)

            graph = add_analyte_edge(
                graph=graph,
                sub_network=sub_network,
                analyte_name=analyte_name,
                reaction_name=reaction_name)

            # Check if analyte is a complex
            if analyte_id_update in complex_reference.keys():
                analyte_complex = analyte_id_update

            elif analyte_id in complex_reference.keys():
                analyte_complex = analyte_id

            else:
                analyte_complex = None

            # If a valid complex ID is pinged, search for complex componenets
            if analyte_complex != None:

                # If analyte is a complex, map all participants
                for complex_component in complex_reference[analyte_complex]['participants']:

                    component_name, component_id = fetch_analyte_info(
                        master_reference=master_reference,
                        analyte_id=complex_component,
                        species_id=species_id)

                    if '(' in component_name \
                    and ')' in component_name:
                        component_name = component_name.split('(')[0]

                    # Add node and edge for participant
                    graph = add_complex_entity(
                        graph=graph,
                        component_name=component_name,
                        component_id=component_id,
                        complex_name=analyte_name)

    return graph

"""Add analyte node
"""
def add_analyte_node(
        graph,
        sub_network, # This starts at the reaction level
        analyte_name,
        analyte_id,
        type):

    graph.add_node(analyte_name) # -> figure out why graph object is turning to tuple and not graph

    graph.nodes()[analyte_name]['name'] = analyte_name
    graph.nodes()[analyte_name]['id'] = analyte_id
    graph.nodes()[analyte_name]['type'] = str(type[:-1]) # trim off plurality of type name

    if type == 'modifiers':
        graph.nodes()[analyte_name]['sub_type'] = sub_network[type][analyte_id]['type']

    else:
        graph.node()[analyte_name]['sub_type'] = None

    graph.nodes()[analyte_name]['stoichiometry'] = sub_network[type][analyte_id]['stoichiometry']

    return graph

"""Add analyte edge
"""
def add_analyte_edge(
        graph,
        sub_network, # This starts at the reaction level
        analyte_name,
        reaction_name):

    # Add modifier edge and sub-type ID
    if graph.nodes()[analyte_name]['sub_type'] != None:
        graph.add_edges_from([
            (analyte_name, reaction_name)])
        graph.edges()[(analyte_name, reaction_name)]['type'] = graph.node()[analyte_name]['sub_type']

    # Add reactant edge
    elif graph.nodes()[analyte_name]['type'] == 'reactant':
        graph.add_edges_from([
            (analyte_name, reaction_name)])
        graph.edges()[(analyte_name, reaction_name)]['type'] = graph.node()[analyte_name]['type']

        # Add reversible reaction bi-directional edges
        if '>' in reaction_name and '<' in reaction_name:
            graph.add_edges_from([
                (reaction_name, analyte_name)])
            graph.edges()[(reaction_name, analyte_name)]['type'] = graph.node()[analyte_name]['type']

    # Add product edge
    elif graph.nodes()[analyte_name]['type'] == 'product':
        graph.add_edges_from([
            (reaction_name, analyte_name)])
        graph.edges()[(reaction_name, analyte_name)]['type'] = graph.node()[analyte_name]['type']

        # Add reversible reaction bi-directional edges
        if '>' in reaction_name and '<' in reaction_name:
            graph.add_edges_from([
                (analyte_name, reaction_name)])
            graph.edges()[(analyte_name, reaction_name)]['type'] = graph.node()[analyte_name]['type']

    # All others
    else:
        graph.add_edges_from([
            (analyte_name, reaction_name)])
        graph.edges()[(analyte_name, reaction_name)]['type'] = 'other'

    return graph

"""Add complex entity node and edge
"""
def add_complex_entity(
        graph,
        component_name,
        component_id,
        complex_name):

    graph.add_node(component_name)
    graph.nodes()[component_name]['name'] = component_name
    graph.nodes()[component_name]['id'] = component_id
    graph.nodes()[component_name]['type'] = 'complex_component'
    graph.nodes()[component_name]['sub_type'] = None
    graph.nodes()[component_name]['stoichiometry'] = None

    graph.add_edges_from([
        (component_name, complex_name)])
    graph.edges()[(component_name, complex_name)]['type'] = 'complex_component'

    return graph

"""Get analyte info from master reference
"""
def fetch_analyte_info(
        master_reference,
        analyte_id,
        species_id):

    if species_id + analyte_id.split('-')[-1] in master_reference.keys():
        analyte_name = master_reference[species_id + analyte_id.split('-')[-1]]
        analyte_id = species_id + analyte_id.split('-')[-1]

    elif analyte_id in master_reference.keys():
        analyte_name = master_reference[analyte_id]

    else:
        analyte_name = species_id + analyte_id.split('-')[-1]
        print('Error: Could not find analyte', analyte_name, 'in databases')

    return analyte_name, analyte_id

"""Map node and edge colors and attributes
"""
def map_graph_attributes(
        graph,
        data):

    # Get max value in dataframe, may want to do this omic to omic
    max_value = max(abs(data[0]))

    # Map node color, size
    for node in graph.nodes():

        # Get node color values
        if graph.nodes()[node]['type'] == 'reaction':
            expression_value = None
            color = 'grey'
            rgba = (0.75, 0.75, 0.75, 1)

        else:

            if node in data.index \
            and str(data.loc[node][0]) != 'nan' \
            and data.loc[node][0] != None:
                expression_value = data.loc[node][0]
                color, rgba = extract_value(
                    expression_value=expression_value,
                    max_value=max_value)

            else:
                expression_value = None
                color = 'white'
                rgba = (1, 1, 1, 1)

        # Add attributes
        graph = add_node_attributes(
            graph=graph,
            node=node,
            color=color,
            rgba=rgba,
            expression_value=expression_value)

    # Map edge color, size
    for edge in graph.edges():

        edge = add_edge_attributes(
            graph=graph,
            edge_id=edge)

    return graph

"""Extract expression value
"""
def extract_value(
        expression_value,
        max_value):

    position = (expression_value + max_value) / (2 * max_value)
    rgba = cmap(position)

    return 'custom', rgba

"""Add node attributes
"""
def add_node_attributes(
        graph,
        node,
        color,
        rgba,
        expression_value):

    if graph.nodes()[node]['type'] == 'reaction':
        graph.nodes()[node]['color'] = 'grey'
        graph.nodes()[node]['rgba'] = (0.75, 0.75, 0.75, 1)
        graph.nodes()[node]['rgba_js'] = convert_rgba(graph.nodes()[node]['rgba'])
        graph.nodes()[node]['expression'] = None
        graph.nodes()[node]['size'] = str(1000)

    elif graph.nodes()[node]['type'] == 'reactant' \
    or graph.nodes()[node]['type'] == 'product':
        graph.nodes()[node]['color'] = [color]
        graph.nodes()[node]['rgba'] = rgba
        graph.nodes()[node]['rgba_js'] = convert_rgba(graph.nodes()[node]['rgba'])
        graph.nodes()[node]['expression'] = str(expression_value)
        graph.nodes()[node]['size'] = str(500)

    else:
        graph.nodes()[node]['color'] = color
        graph.nodes()[node]['rgba'] = rgba
        graph.nodes()[node]['rgba_js'] = convert_rgba(graph.nodes()[node]['rgba'])
        graph.nodes()[node]['expression'] = str(expression_value)
        graph.nodes()[node]['size'] = str(300)

    return graph

"""Add edge attributes
"""
def add_edge_attributes(
        graph,
        edge_id):

    # Set color attributes
    if graph.edges()[edge_id]['type'] == 'catalyst'\
    or graph.edges()[edge_id]['type'] == 'positiveregulator':
        graph.edges()[edge_id]['color'] = 'green'
        graph.edges()[edge_id]['rgba'] = (0, .5, 0, 1)
        graph.edges()[edge_id]['rgba_js'] = convert_rgba(graph.edges()[edge_id]['rgba'])

    elif graph.edges()[edge_id]['type'] == 'inhibitor' \
    or graph.edges()[edge_id]['type'] == 'negativeregulator':
        graph.edges()[edge_id]['color'] = 'red'
        graph.edges()[edge_id]['rgba'] = (1, 0, 0, 1)
        graph.edges()[edge_id]['rgba_js'] = convert_rgba(graph.edges()[edge_id]['rgba'])

    elif graph.edges()[edge_id]['type'] == 'complex_component':
        graph.edges()[edge_id]['color'] = 'purple'
        graph.edges()[edge_id]['rgba'] = (0.5, 0, 0.5, 0.9)
        graph.edges()[edge_id]['rgba_js'] = convert_rgba(graph.edges()[edge_id]['rgba'])

    else:
        graph.edges()[edge_id]['color'] = 'grey'
        graph.edges()[edge_id]['rgba'] = (0.75, 0.75, 0.75, 1)
        graph.edges()[edge_id]['rgba_js'] = convert_rgba(graph.edges()[edge_id]['rgba'])

    return graph

"""Output graph information to JSON file
"""
def output_graph(
        graph,
        output_name):

    data = json_graph.node_link_data(G)

    with open(output_name, 'w') as f:
        json.dump([data], f, indent=4) # Parse out as array for javascript

"""Convert python RGBA tuple to web-friendly tuple for later viz
"""
def convert_rgba(
        rgba_tuple):

    rgba_list = list(rgba_tuple)

    rgba_new = []
    for x in rgba_list[:3]:
        rgba_new.append(int(x * 255))

    rgba_new.append(rgba_list[3])
    return tuple(rgba_new)

"""Format graph attributes for networkX plotting
"""
def layout_graph(
        graph):

    node_list = []
    node_color = []
    node_size = []

    for node in graph.nodes():

        node_list.append(node)
        node_color.append(
            literal_eval(
                strip_string(graph.nodes()[node]['rgba'])))
        node_size.append(
            literal_eval(
                strip_string(graph.nodes()[node]['size'])))

    edge_list = []
    edge_color = []

    for edge in graph.edges():

        edge_list.append(edge)
        edge_color.append(
            literal_eval(
                strip_string(graph.edges()[edge]['rgba'])))

    k_value = float(3 / sqrt(len(node_list)))

    node_positions = nx.spring_layout(
        graph,
        k=k_value,
        iterations=20,
        scale=10)

    position_reference = {
        'nodes': {
            'node_list': node_list,
            'node_color': node_color,
            'node_size': node_size
        },
        'edges': {
            'edge_list': edge_list,
            'edge_color': edge_color
        }
    }

    return node_positions, position_reference

"""Unlist a string from list format required for JSON export of graph
"""
def strip_string(
        __str__):

    return str(__str__).strip('[]').replace('\'', '')

"""Plot graph using networkX
"""
def plot_graph(
        graph,
        output):

    positions, position_ref = layout_graph(
        graph=G)

    fig, axes = plt.subplots(
        nrows = 1,
        ncols = 1,
        figsize = (25, 15)) # Create shared axis for cleanliness

    nx.draw(
        graph,
        positions,
        ax=axes,
        nodelist=position_ref['nodes']['node_list'],
        node_color=position_ref['nodes']['node_color'],
        node_size=position_ref['nodes']['node_size'],
        edgelist=position_ref['edges']['edge_list'],
        edge_color=position_ref['edges']['edge_color'],
        with_labels=True,
        font_weight='bold',
        arrowsize=20,
        font_size=8)

    axes.collections[0].set_edgecolor('black')

    plt.savefig(
        output,
        dpi=600,
        bbox_inches='tight')

"""Run graph generation
TODO:
- How to handle metabolite and gene with same name?
    - MAL is metabolite and gene abbreviation
"""
def __main__(
        data, # Assumed to already be name mapped
        network,
        pathways,
        species_id,
        output,
        black_list=[],
        plot=False,
        graph_name='name'):

    # Get output directory for graphs
    graph_directory = create_graph_output(
        output)

    # Make master reference
    master_reference = network['master_reference']
    complex_reference = network['complexes_reference']['complex_mapper']

    # Generate graph
    # Input list of pathway names
    for p in pathways:

        print('Analyzing', p)

        pathway_key = network['pathway_types'][p]

        for k in pathway_key:

            subnetwork = network['pathways'][k]

            name = 'graph-' + subnetwork['reactome_id']

            if graph_name == 'reactome':
                graph_name = graph_directory + name + 'fake.json'
                plot_name = graph_directory + name + '.pdf'

            else:
                graph_name = graph_directory + p.replace(' ', '_') + 'fake.json'
                plot_name = graph_directory + p.replace(' ', '_') + '.pdf'

            species_accession = 'R-' + species_id + '-'

            # Build Graph
                # Add reactions
                # Add analytes
                # Add edges
            G = build_graph(
                sub_network=subnetwork,
                master_reference=master_reference,
                complex_reference=complex_reference,
                species_id=species_accession,
                black_list=black_list)

            # Add node and edge colors and size -> inside runs a function for node and for edge coloring
            G = map_graph_attributes(
                graph=G,
                data=data)

            # Output graph
            output_graph(
                graph=G,
                output_name=graph_name)

            # Plot graph
            if plot != False:
                plot_graph(
                    graph=G,
                    output=plot_name)

"""Test data
"""
def test():

    ################
    output = '/Users/jordan/Desktop/Metaboverse/tests/analysis_tests/'
    with open(output + 'HSA_metaboverse_db.pickle', 'rb') as network_file:
        reference = pickle.load(network_file)
    network = reference

    for x in network['pathway_types'].keys():
        if 'glutamate' in x.lower():
            print(x)

    pathways = [
        'Urea cycle',
        'Citric acid cycle (TCA cycle)',
        'Alanine metabolism',
        'Glutamate and glutamine metabolism']
    species_id = 'HSA'

    trans = pd.read_csv(
        '/Users/jordan/Desktop/danielle_seq/analysis/danielle_count_table_summed_diffx_named.tsv',
        sep='\t',
        index_col=2)
    del trans.index.name
    trans = trans[['log2FoldChange']]
    trans.columns = [0]

    trans.head()

    from random import randint

    addi = []
    for x in range(len(trans[0].tolist())):
        addi.append(randint(-1,1))

    trans[0] = trans[0] + addi
    trans.head()

    metabol = pd.read_csv(
        '/Users/jordan/Desktop/danielle_seq/PGC1_beta_metabolomics.txt',
        sep='\t',
        index_col=0)
    metabol = metabol[['log2foldchange']]
    metabol.columns = [0]
    metabol[0] = [0.053, 1.537, -1.114, -2.709, -1.67, 0.89, -0.3]
    metabol.loc['CAP'] = [1.2]
    metabol.loc['Urea'] = [3]
    metabol

    data = pd.concat([trans, metabol])
    data = data.dropna()

    black_list = ['H2O', 'ATP', 'ADP', 'H+', 'GDP', 'GTP']
    plot = True
    graph_name='name'
    #################
