"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
    Copyright (C) 2019  Youjia Zhou, Jordan A. Berg
    zhou325 <at> sci <dot> utah <dot> edu
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
import sys
import re
import json

# This will need to be modified to handle timecourse data
index_number = 0

def process_data(
        model_file):
    """Generate network format for motif searching
    """

    with open(model_file) as json_file:
        network_data = json.load(json_file)

    nodes = network_data['nodes']
    links = network_data['links']
    pathway_dict = network_data['collapsed_pathway_dictionary']
    reactions_dict = network_data['collapsed_reaction_dictionary']

    reaction_nodes = []
    nodes_dict = {}

    for node in nodes:
        if 'notes' in node.keys():
            del node['notes']
        nodes_dict[node['id']] = node
        if node['type'] == 'reaction':
            reaction_nodes.append(node)

    ### only keep reactant & product links ###
    links_new = []
    for link in links:
        if link['type'] in ['reactant', 'product']:
            links_new.append(link)

    ### assign links to reaction nodes ###
    reaction_nodes_dict = {}
    for node in reaction_nodes:
        node['links'] = {'reactant':[], 'product':[]}
        node['pathways'] = []
        reaction_nodes_dict[node['id']] = node
    for link in links_new:
        if link['source'] in reaction_nodes_dict.keys():
            reaction_nodes_dict[link['source']]['links']['product'].append(link)
        elif link['target'] in reaction_nodes_dict.keys():
            reaction_nodes_dict[link['target']]['links']['reactant'].append(link)
    ### assign pathways to reaction nodes ###
    for p_key in pathway_dict:
        for node_id in pathway_dict[p_key]['reactions']:
            if node_id in reaction_nodes_dict.keys():
                reaction_nodes_dict[node_id]['pathways'].append(p_key)
            else:
                pathway_dict[p_key]['reactions'].remove(node_id)

    ### other nodes: only reactant & product nodes ###
    other_nodes_dict = {}
    for node_id in reaction_nodes_dict:
        react_node = reaction_nodes_dict[node_id]

        for link in react_node['links']['reactant']:
            other_node = nodes_dict[link['source']]
            if type(other_node['values']) == list \
            and len(other_node['values']) > 0:
                other_node['expression'] = other_node['values'][index_number]
            else:
                other_node['expression'] = "None"
            other_nodes_dict[link['source']] = other_node

        for link in react_node['links']['product']:
            other_node = nodes_dict[link['target']]
            if type(other_node['values']) == list \
            and len(other_node['values']) > 0:
                other_node['expression'] = other_node['values'][index_number]
            else:
                other_node['expression'] = "None"
            other_nodes_dict[link['target']] = other_node

    network_data['react_nodes'] = reaction_nodes_dict
    network_data['other_nodes'] = other_nodes_dict
    network_data['pathways_motif'] = pathway_dict

    # Write updated network model to same file as before
    with open(model_file,'w') as f:
        json.dump(network_data, f)

def __main__(
        model_file):

    process_data(
        model_file=model_file)
    print("Motif search complete.")
