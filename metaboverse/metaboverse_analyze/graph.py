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
import json
import pickle
import networkx as nx
from networkx.readwrite import json_graph
import matplotlib
import matplotlib.pyplot as plt

cmap = matplotlib.cm.get_cmap('RdYlBu')

def create_graph_output(
        output):

    graph_directory = output + 'graphs/'

    if not os.path.exists(graph_directory):
        os.makedirs(graph_directory)

    return graph_directory

def add_expression_average(
        graph,
        data,
        expression_average,
        analyte):

    if expression_value != 0:
        max_value = max(abs(data[0]))
        position = (expression_average + max_value) / (2 * max_value)
        color_value = cmap(position)
        graph.node()[analyte]['color'] = color_value
        graph.node()[analyte]['expression_value'] = expression_average

    else:
        graph.node()[analyte]['color'] = 'white'
        graph.node()[analyte]['expression_value'] = None

def add_expression_value(
        graph,
        data,
        analyte,
        analyte_id):

    if analyte_id in data.index:
        max_value = max(abs(data[0])) # Get absolute max value in dataset, may want to do this omic to omic
        value = data.loc[analyte_id][0]
        position = (value + max_value) / (2 * max_value)
        color_value = cmap(position)
        graph.node()[analyte]['color'] = color_value
        graph.node()[analyte]['expression_value'] = value

    # If not in user data, color white
    else:
        graph.node()[analyte]['color'] = 'white'
        graph.node()[analyte]['expression_value'] = None

    return graph, value

def add_component_node_edge(
        graph,
        master_reference,
        species_id,
        data,
        analyte,
        complex_component):

     component_name, component_id = fetch_analyte_info(
        master_reference=master_reference,
        species_id=species_id,
        analyte=analyte)

    graph.add_node(complex_component)
    graph.node()[complex_component]['type'] = 'complex_component'
    graph.node()[complex_component]['stoichiometry'] = None

    graph.add_edges_from([
        (component_name, analyte)])
    graph.edges()[(component_name, analyte)]['color'] = 'purple'
    graph.edges()[(component_name, analyte)]['type'] = 'complex_component'

    graph, expression_value = add_expression_value(
        graph=graph,
        data=data,
        analyte=component_name,
        analyte_id=component_id)

    return graph, expression_value

def fetch_analyte_info(
        master_reference,
        species_id,
        analyte):

    if analyte in master_reference.keys():
        analyte_name = master_reference[analyte]
        analyte_id = analyte

    elif species_id + analyte.split('-')[-1] in master_reference.keys():
        analyte_name = master_reference[species_id + analyte.split('-')[-1]]
        analyte_id = species_id + analyte.split('-')[-1]

    else:
        analyte_name = analyte
        analyte_id = analyte

    return analyte_name, analyte_id

def add_node_edge(
        graph,
        reaction_name,
        sub_network,
        master_reference,
        complex_reference,
        data,
        species_id,
        type,
        black_list):

    for key in sub_network['reactions'][reaction_name][type].keys():

        # Get analyte name if possible
        analyte = sub_network['reactions'][reaction_name][type][key]['species_id']

        analyte_name, analyte_id = fetch_analyte_info(
            master_reference=master_reference,
            species_id=species_id,
            analyte=analyte)

        if analyte_name in black_list:
            pass

        else:
            analyte_list.append(analyte_name)

            # Add node and edge to reaction
            graph.add_node(analyte_name)
            graph.node()[analyte_name]['type'] = type[:-1] # trim off plurality of type name
            graph.node()[analyte_name]['stoichiometry'] = sub_network['reactions'][reaction_name][type][key]['stoichiometry']

            if type == 'modifiers':

                mod_type = process['reactions'][reaction_name][type][key]['type']
                graph.node()[analyte_name]['sub_type'] = mod_type.lower()
                graph.node()[analyte_name]['node_size'] = 300

            else:
                graph.node()[analyte_name]['node_size'] = 500

            # Add color for analyte (dark blue, light blue, yellow, orange, red)
            graph, expression_value = add_expression_value(
                graph=graph,
                data=data,
                analyte=analyte_name,
                analyte_id=analyte_id)

            if type == 'modifiers':

                mod_type = process['reactions'][reaction_name][type][key]['type']

                graph.add_edges_from([
                    (analyte_name, reaction_name)])

                if mod_type.lower() == 'catalyst' or mod_type.lower() == 'positiveregulator':
                    graph.edges()[(analyte_name, reaction_name)]['color'] = 'green'
                    graph.edges()[(analyte_name, reaction_name)]['type'] = 'catalyst'

                elif mod_type.lower() == 'inhibitor' or mod_type.lower() == 'negativeregulator':
                    graph.edges()[(analyte_name, reaction_name)]['color'] = 'red'
                    graph.edges()[(analyte_name, reaction_name)]['type'] = 'inhibitor'

                else:
                    graph.edges()[(analyte_name, reaction_name)]['color'] = 'lightblue'
                    graph.edges()[(analyte_name, reaction_name)]['type'] = 'other_modifier'

            elif '>' in reaction_name and '<' in reaction_name:
                graph.add_edges_from([
                    (analyte_name, reaction_name)])
                graph.edges()[(analyte_name, reaction_name)]['color'] = 'grey'
                graph.edges()[(analyte_name, reaction_name)]['type'] = type[:-1] + ' -> reaction'

                graph.add_edges_from([
                    (reactant_name, analyte_name)])
                graph.edges()[(reactant_name, analyte_name)]['color'] = 'grey'
                graph.edges()[(reactant_name, analyte_name)]['type'] = 'reaction -> ' + type[:-1]

            elif '<' in reaction:
                graph.add_edges_from([
                    (reactant_name, analyte_name)])
                graph.edges()[(reaction, analyte_name)]['color'] = 'grey'
                graph.edges()[(reaction, analyte_name)]['type'] = 'reaction -> ' + type[:-1]

            else:
                graph.add_edges_from([
                    (analyte_name, reaction)])
                graph.edges()[(analyte_name, reaction)]['color'] = 'grey'
                graph.edges()[(analyte_name, reaction)]['type'] = type[:-1] + ' -> reaction'

        # Add protein complex nodes and edges
        # Need to figure out mapping here
        if analyte_id in complex_reference.keys():

            value_list = []

            for complex_component in complex_reference[analyte_id]['participants']:

                graph, expression_value = add_component_node_edge(
                    graph=graph,
                    master_reference=master_reference,
                    species_id=species_id,
                    data=data,
                    analyte=analyte_id,
                    complex_component=complex_component)
                value_list.append(expression_value)

            expression_average = sum(value_list) / len(value_list)
            graph, expression_value = add_expression_average(
                graph=graph,
                data=data,
                expression_average=expression_average,
                analyte=component_name)

            #for complex_partner in complex_reference[analyte_id]['partner_complexes']:
            # maybe add this in the future

    return graph

def add_reaction_node(
        graph,
        reaction_name):

    reaction = sub_network['reactions'][reaction_name]['name']
    graph.add_node(reaction)
    graph.node()[reaction]['type'] = 'reaction'
    graph.node()[reaction]['color'] = 'grey'
    graph.node()[reaction]['expression_value'] = None
    graph.node()[reaction]['node_size'] = 1000
    graph.node()[reaction]['stoichiometry'] = None

    return graph, reaction

def process_reactions(
        graph,
        reaction_name,
        sub_network,
        master_reference,
        complex_reference,
        data,
        species_id,
        black_list):

    # add reaction node
    graph, reaction = add_reaction_node(
        graph=graph,
        reaction_name=reaction_name)

    # Add reactant nodes and edges
    graph = add_node_edge(
        graph=graph,
        reaction_name=reaction,
        sub_network=sub_network,
        master_reference=master_reference,
        complex_reference=complex_reference,
        data=data,
        species_id=species_id,
        type='reactants',
        black_list=black_list)

    # Add product nodes and edges
    graph = add_node_edge(
        graph=graph,
        reaction_name=reaction,
        sub_network=sub_network,
        master_reference=master_reference,
        complex_reference=complex_reference,
        data=data,
        species_id=species_id,
        type='products',
        black_list=black_list)

    # Add modifier nodes and edges
    graph = add_node_edge(
        graph=graph,
        reaction_name=reaction,
        sub_network=sub_network,
        master_reference=master_reference,
        complex_reference=complex_reference,
        data=data,
        species_id=species_id,
        type='modifiers',
        black_list=black_list)

    return graph

def layout_graph(
        graph):

    node_list = []
    node_color = []
    node_size = []

    for n in graph.nodes():

        node_list.append(n)
        node_color.append(graph.nodes()[n]['color'])
        node_size.append(graph.nodes()[n]['node_size'])

    edge_list = []
    edge_color = []

    for e in graph.edges():

        edge_list.append(e)
        edge_color.append(graph.edges()[e]['color'])

    node_positions = nx.spring_layout(
        graph,
        k=(3/sqrt(len(node_list))),
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

def build_graph(
        sub_network,
        master_reference,
        complex_reference,
        data,
        species_id,
        black_list):

    G = nx.DiGraph()

    # cycle through reactions in process
    for key in sub_network['reactions'].keys():

        G = process_reactions(
            graph=G,
            reaction_name=key,
            sub_network=sub_network,
            master_reference=master_reference,
            complex_reference=complex_reference,
            data=data,
            species_id=species_id,
            black_list=black_list)

    return G

def output_graph(
        graph,
        output_name):

    data = json_graph.node_link_data(graph)

    with open(output_name, 'w') as f:
        json.dump(data, f, indent=4)

def __main__(
        data,
        network,
        pathways,
        species_id,
        output,
        black_list,
        plot=False):

    # Get output directory for graphs
    graph_directory = create_graph_output(
        output)

    # Make master reference
    master_reference = network['master_reference']
    complex_reference = network['complex_mapper']

    # Generate graph
    for p in pathways:

        print('Analyzing', p)

        pathway_key = network['pathway_types'][p]
        process = network['pathways'][pathway_key]

        name = 'graph_' + species_id + '_' + process
        graph_name = graph_directory + name + '.json'
        plot_name = graph_directory + name + '.pdf'

        species_id = 'R-' + species_id + '-'

        G = build_graph(
            sub_network=process,
            master_reference=master_reference,
            complex_reference=complex_reference,
            data=data,
            species_id=species_id,
            black_list=black_list)

        output_graph(
            graph=G,
            output_name=graph_name)

        if plot != False:

            plot_graph(
                graph=G,
                output=plot_name)

    print('Graphs output to', graph_directory)
