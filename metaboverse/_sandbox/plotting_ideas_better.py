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
import pandas as pd
import numpy as np
import networkx as nx
import pickle
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt

"""Import internal dependencies
"""

"""Set globals
"""

"""Read in networkx-formatted pickle file from curation step
"""
def read_network(
        file):

    with open(file, 'rb') as network_file:
        network = pickle.load(network_file)

    return network

"""Build custom metabolic model for use in analyses
"""
def __main__(
        args_dict):

    # Read in network file
    network = read_network(
        file=args_dict['model'])

    return network



args_dict = {'model': '/Users/jordan/Desktop/reactome_test/HSA_metaboverse_db.pickle'}
network = __main__(args_dict)

for x in list(network['pathways_types'].keys()):
    if 'tca' in x.lower():
        print(x)

""" Test blacklisting metabolites
"""
black_list = [
    'H+ ',
    'H2O ',
    'CO2 ']


""" Test data
"""
analytes = [
    'Ac-CoA ',
    'OAA ',
    'H2O ',
    'CoA-SH ',
    'CIT ',
    'CS dimer ',
    #'H+ ',
    'NADH ',
    'OAA ',
    'MAL ',
    'NAD+ ',
    'MDH2 dimer ',
    'NADP+ ',
    'MAL ',
    #'H+ ',
    'PYR ',
    'NADPH ',
    'CO2 ',
    'ME3:Mg2+ tetramer ',
    #'H+ ',
    'NAD+ ',
    'NADPH ',
    'NADP+ ',
    #'H+ ',
    'NADH ',
    'NNT dimer ',
    'SUCC-CoA ',
    'Pi ',
    'ADP ',
    'CoA-SH ',
    'ATP ',
    'SUCCA ',
    'SUCLA2:SUCLG1 ',
    'CoA-SH ',
    #'2OG ',
    'NAD+ ',
    'SUCC-CoA ',
    'NADH ',
    'CO2 ',
    'lipo-aKGDH ',
    'MAL ',
    'NAD+ ',
    #'H+ ',
    'NADH ',
    'OAA ',
    'MDH2 dimer ',
    'ISCIT ',
    'CIT ',
    'ACO2 ',
    'MAL ',
    #'FUMA ',
    'H2O ',
    'FH tetramer ',
    #'FUMA ',
    'H2O ',
    'MAL ',
    'FH tetramer ',
    'NAD+ ',
    'ISCIT ',
    #'2OG ',
    #'H+ ',
    'NADH ',
    'CO2 ',
    'IDH3 complex ',
    'ADP ',
    'NADP+ ',
    'ISCIT ',
    #'2OG ',
    #'H+ ',
    'NADPH ',
    'CO2 ',
    'IDH2 dimer ',
    'OAA ',
    'PYR ',
    'CO2 ',
    'FAHD1:Mg2+ dimer ',
    'MAL ',
    'NAD+ ',
    #'H+ ',
    'PYR ',
    'NADH ',
    'CO2 ',
    'ME2:Mg2+ tetramer ',
    #'FUMA ',
    'ATP ',
    'CIT ',
    'ISCIT ',
    'ACO2 ',
    #'SUCCA ',
    #'FUMA ',
    'SDH complex (ox.) ',
    'SUCC-CoA ',
    'Pi ',
    'GDP ',
    'CoA-SH ',
    'GTP ',
    'SUCCA ',
    'SUCLG1:SUCLG2 ']

import random

test_data = []

for x in set(analytes):

    value = random.uniform(-5, 5)

    test_data.append([x, value])

import pandas as pd
test_data = pd.DataFrame(test_data)
test_data.index = test_data[0]
test_data = test_data.drop(0, axis=1)
del test_data.index.name
test_data.columns = [0]

test_data

""" Test combine refs
"""
master_reference = {}

reference_list = [
    'chebi_reference',
    'uniprot_reference',
    'ensembl_reference',
    'ncbi_reference',
    'mirbase_reference']

for x in reference_list:

    for key in network[x].keys():

        master_reference[network[x][key]['analyte_id']] = network[x][key]['analyte']

for key in network['complexes_reference'].keys():
    if key == 'R-HSA-156627':
        print(key)
    master_reference[network['complexes_reference'][key]['complex_id']] = network['complexes_reference'][key]['complex_name']


""" Test run
"""

len(network['pathways'].keys())

for key in network['complexes_reference'].keys():
    complex = network['complexes_reference'][key]['complex_name']
    if 'mico' in complex.lower():
        print(complex)
        print(network['complexes_reference'][key])

network['complexes_reference']['R-HSA-8949617']

network['pathways']['R-HSA-8949613']

master_reference['R-HSA-8949617']
master_reference['R-HSA-8949570']

# init
species_id = 'R-HSA-'
process_name = 'Citric acid cycle (TCA cycle)'
pathway_key = network['pathways_types'][process_name]
database = network['pathways']
analyte_list = []
cmap = matplotlib.cm.get_cmap('RdYlBu')

cmap

if len(pathway_key) == 1:
    process = database[pathway_key[0]]

# Start graphing
G = nx.DiGraph()

# cycle through reactions in process
for key in process['reactions'].keys():

    # add reaction node
    reaction = process['reactions'][key]['name']
    G.add_node(reaction)
    G.node()[reaction]['type'] = 'reaction'
    G.node()[reaction]['color'] = 'grey'
    G.node()[reaction]['size'] = 1000
    G.node()[reaction]['stoichiometry'] = 'N/A'

    # Add reactant nodes and edges
    for sub_key in process['reactions'][key]['reactants'].keys():

        # Get analyte name if possible
        reactant = process['reactions'][key]['reactants'][sub_key]['species_id']

        if reactant in master_reference.keys():
            reactant_name = master_reference[reactant]

        elif species_id + reactant.split('-')[-1] in master_reference.keys():
            reactant_name = master_reference[species_id + reactant.split('-')[-1]]

        else:
            reactant_name = reactant

        if reactant_name in black_list:
            pass
        else:
            analyte_list.append(reactant_name)

            # Add node and edge to reaction
            G.add_node(reactant_name)
            G.node()[reactant_name]['type'] = 'reactant'
            G.node()[reactant_name]['stoichiometry'] = process['reactions'][key]['reactants'][sub_key]['stoichiometry']

            # Add color for analyte (dark blue, light blue, yellow, orange, red)
            if reactant_name in test_data.index:
                max_value = max(abs(test_data[0]))
                value = test_data.loc[reactant_name][0]
                position = (value + max_value) / (2 * max_value)
                color_value = cmap(position)
                G.node()[reactant_name]['color'] = color_value

            # If not in user data, color white
            else:
                G.node()[reactant_name]['color'] = 'white'

            G.node()[reactant_name]['size'] = 500

            if '>' in reaction and '<' in reaction:
                G.add_edges_from([
                    (reactant_name, reaction)])
                G.edges()[(reactant_name, reaction)]['color'] = 'grey'

                G.add_edges_from([
                    (reaction, reactant_name)])
                G.edges()[(reaction, reactant_name)]['color'] = 'grey'

            elif '<' in reaction:
                G.add_edges_from([
                    (reaction, reactant_name)])
                G.edges()[(reaction, reactant_name)]['color'] = 'grey'

            else:
                G.add_edges_from([
                    (reactant_name, reaction)])
                G.edges()[(reactant_name, reaction)]['color'] = 'grey'

    # Add product nodes and edges
    for sub_key in process['reactions'][key]['products'].keys():

        # Get analyte name if possible
        product = process['reactions'][key]['products'][sub_key]['species_id']

        if product in master_reference.keys():
            product_name = master_reference[product]

        elif species_id + product.split('-')[-1] in master_reference.keys():
            product_name = master_reference[species_id + product.split('-')[-1]]

        else:
            product_name = product

        if product_name in black_list:
            pass

        else:
            analyte_list.append(product_name)

            # Add node and edge to reaction
            G.add_node(product_name)
            G.node()[product_name]['type'] = 'product'
            G.node()[product_name]['stoichiometry'] = process['reactions'][key]['products'][sub_key]['stoichiometry']

            if product_name in test_data.index:
                max_value = max(abs(test_data[0]))
                value = test_data.loc[product_name][0]
                position = (value + max_value) / (2 * max_value)
                color_value = cmap(position)
                G.node()[product_name]['color'] = color_value

            # If not in user data, color white
            else:
                G.node()[product_name]['color'] = 'white'

            G.node()[product_name]['size'] = 500

            if '>' in reaction and '<' in reaction:
                G.add_edges_from([
                    (reaction, product_name)])
                G.edges()[(reaction, product_name)]['color'] = 'grey'

                G.add_edges_from([
                    (product_name, reaction)])
                G.edges()[(product_name, reaction)]['color'] = 'grey'

            elif '<' in reaction:
                G.add_edges_from([
                    (product_name, reaction)])
                G.edges()[(product_name, reaction)]['color'] = 'grey'

            else:
                G.add_edges_from([
                    (reaction, product_name)])
                G.edges()[(reaction, product_name)]['color'] = 'grey'

    # Add catalyst and inhibitor nodes and edges
    for sub_key in process['reactions'][key]['modifiers'].keys():

        # Get analyte name if possible
        modifier = process['reactions'][key]['modifiers'][sub_key]['species_id']
        type = process['reactions'][key]['modifiers'][sub_key]['type']

        if modifier in master_reference.keys():
            modifier_name = master_reference[modifier]

        elif species_id + modifier.split('-')[-1] in master_reference.keys():
            modifier_name = master_reference[species_id + modifier.split('-')[-1]]

        else:
            modifier_name = modifier

        if modifier_name in black_list:
            pass

        else:
            analyte_list.append(modifier_name)

            # Add node and edge to reaction
            G.add_node(modifier_name)
            G.node()[modifier_name]['type'] = type.lower() # This will be overwritten if another process uses this as a reactant for example
            G.node()[modifier_name]['stoichiometry'] = process['reactions'][key]['modifiers'][sub_key]['stoichiometry']

            if modifier_name in test_data.index:
                max_value = max(abs(test_data[0]))
                value = test_data.loc[modifier_name][0]
                position = (value + max_value) / (2 * max_value)
                color_value = cmap(position)
                G.node()[modifier_name]['color'] = color_value

            # If not in user data, color white
            else:
                G.node()[modifier_name]['color'] = 'white'

            G.node()[modifier_name]['size'] = 300

            G.add_edges_from([
                (modifier_name, reaction)])

            if type.lower() == 'catalyst' or type.lower() == 'positiveregulator':
                G.edges()[(modifier_name, reaction)]['color'] = 'green'

            elif type.lower() == 'inhibitor' or type.lower() == 'negativeregulator':
                G.edges()[(modifier_name, reaction)]['color'] = 'red'

            else:
                G.edges()[(modifier_name, reaction)]['color'] = 'lightblue'



    # Add protein complex nodes and edges
node_list = []
node_color = []
node_size = []
for n in G.nodes():

    node_list.append(n)
    node_color.append(G.nodes()[n]['color'])
    node_size.append(G.nodes()[n]['size'])

edge_list = []
edge_color = []
for e in G.edges():

    edge_list.append(e)
    edge_color.append(G.edges()[e]['color'])

pos = nx.spring_layout(
    G,
    k=(3/sqrt(len(node_list))),
    iterations=20,
    scale=10)

fig, axes = plt.subplots(
        nrows = 1,
        ncols = 1,
        figsize = (25, 15)) # Create shared axis for cleanliness

nx.draw(
    G,
    pos,
    ax=axes,
    nodelist=node_list,
    node_color=node_color,
    node_size=node_size,
    edgelist=edge_list,
    edge_color=edge_color,
    with_labels=True,
    font_weight='bold',
    arrowsize=20,
    font_size=8)
axes.collections[0].set_edgecolor('black')
plt.savefig('/Users/jordan/Desktop/reactome_test/metaboverse_examples/metaboverse_prototype_TCA.pdf', dpi = 600, bbox_inches='tight')


"""
G.nodes
name = id
type = reaction/reactant/product, etc -> reactions names are not displayed natively, need hover
source = Origin entity
target = target entity
shading = grey or white or expression
edge_type = reaction, genes, catalyst, inhibitor
size = size type
stoichiometry = #

Only run D3 for networks user wants to viz
For unbiased analysis, make network, analyze breakpoints, display interesting break point hubs
"""
import json
from networkx.readwrite import json_graph
data = json_graph.node_link_data(G)
with open('/Users/jordan/Desktop/Metaboverse/metaboverse/_sandbox/graph.json', 'w') as f:
    json.dump(data, f, indent=4)
