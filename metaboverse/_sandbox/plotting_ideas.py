"""License Information
Metabo-verse:
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

if len(pathway_key) == 1:
    process = database[pathway_key[0]]

# Start graphing
G = nx.DiGraph()

# cycle through reactions in process
node_colors = []
edge_colors = []

for key in process['reactions'].keys():

    # add reaction node
    reaction = process['reactions'][key]['name']
    G.add_node(reaction)
    node_colors.append(['reaction', reaction, 'grey', 1000])

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

            # Add color for analyte (dark blue, light blue, yellow, orange, red)
            if reactant_name in test_data.index:
                max_value = max(abs(test_data[0]))
                value = test_data.loc[reactant_name][0]
                position = (value + max_value) / (2 * max_value)
                color_value = cmap(position)
                node_colors.append(['reactant', reactant_name, color_value, 500])

            # If not in user data, color white
            else:
                node_colors.append(['reactant', reactant_name, 'white', 500])

            if '>' in reaction and '<' in reaction:
                G.add_edges_from([
                    (reactant_name, reaction)])
                edge_colors.append([(reactant_name, reaction), 'grey'])

                G.add_edges_from([
                    (reaction, reactant_name)])
                edge_colors.append([(reaction, reactant_name), 'grey'])

            elif '<' in reaction:
                G.add_edges_from([
                    (reaction, reactant_name)])
                edge_colors.append([(reaction, reactant_name), 'grey'])

            else:
                G.add_edges_from([
                    (reactant_name, reaction)])
                edge_colors.append([(reactant_name, reaction), 'grey'])

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

            if product_name in test_data.index:
                max_value = max(abs(test_data[0]))
                value = test_data.loc[product_name][0]
                position = (value + max_value) / (2 * max_value)
                color_value = cmap(position)
                node_colors.append(['product', product_name, color_value, 500])

            # If not in user data, color white
            else:
                node_colors.append(['product', product_name, 'white', 500])

            if '>' in reaction and '<' in reaction:
                G.add_edges_from([
                    (reaction, product_name)])
                edge_colors.append([(reaction, product_name), 'grey'])

                G.add_edges_from([
                    (product_name, reaction)])
                edge_colors.append([(product_name, reaction), 'grey'])

            elif '<' in reaction:
                G.add_edges_from([
                    (product_name, reaction)])
                edge_colors.append([(product_name, reaction), 'grey'])

            else:
                G.add_edges_from([
                    (reaction, product_name)])
                edge_colors.append([(reaction, product_name), 'grey'])

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

            if modifier_name in test_data.index:
                max_value = max(abs(test_data[0]))
                value = test_data.loc[modifier_name][0]
                position = (value + max_value) / (2 * max_value)
                color_value = cmap(position)
                node_colors.append(['modifier', modifier_name, color_value, 300])

            # If not in user data, color white
            else:
                node_colors.append(['modifier', modifier_name, 'white', 300])

            if type.lower() == 'catalyst' or type.lower() == 'positiveregulator':
                edge_colors.append([(modifier_name, reaction), 'green'])

            elif type.lower() == 'inhibitor' or type.lower() == 'negativeregulator':
                edge_colors.append([(modifier_name, reaction), 'red'])

            else:
                edge_colors.append([(modifier_name, reaction), 'lightblue'])

            G.add_edges_from([
                (modifier_name, reaction)])

    # Add protein complex nodes and edges

n = [x[1] for x in node_colors]
n_c = [x[2] for x in node_colors]
n_s = [x[3] for x in node_colors]

e = [x[0] for x in edge_colors]
e_c = [x[1] for x in edge_colors]

plt.figure(121,figsize=(20,15))
pos = nx.spring_layout(
    G,
    k=(3/sqrt(len(n))),
    iterations=20,
    scale=10)
nx.draw(
    G,
    pos,
    nodelist=n,
    node_color=n_c,
    node_size=n_s,
    edgelist=e,
    edge_color=e_c,
    with_labels=True,
    font_weight='bold',
    arrowsize=20,
    font_size=8)
ax = plt.gca() # to get the current axis
ax.collections[0].set_edgecolor('black')
plt.savefig('/Users/jordan/Desktop/reactome_test/metaboverse_examples/metaboverse_prototype_TCA.pdf', dpi = 600, bbox_inches='tight')
