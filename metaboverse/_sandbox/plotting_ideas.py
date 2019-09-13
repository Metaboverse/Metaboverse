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


G = nx.DiGraph()

G.add_node('R-ALL-29356')
G.add_node('R-ALL-73588')

G.add_node('beta-ureidopropionate + H2O => beta-alanine + NH4+ + CO2')

G.add_node('R-ALL-113528')
G.add_node('beta-alanine')
G.add_node('R-ALL-31633')


G.add_edges_from([
    ('R-ALL-29356', 'beta-ureidopropionate + H2O => beta-alanine + NH4+ + CO2'),
    ('R-ALL-73588', 'beta-ureidopropionate + H2O => beta-alanine + NH4+ + CO2'),
    ('beta-ureidopropionate + H2O => beta-alanine + NH4+ + CO2', 'R-ALL-113528'),
    ('beta-ureidopropionate + H2O => beta-alanine + NH4+ + CO2', 'beta-alanine'),
    ('beta-ureidopropionate + H2O => beta-alanine + NH4+ + CO2', 'R-ALL-31633')
    ])

G.add_node('beta-alanine')

G.add_node('Mitochondrial uptake of beta-alanine')

G.add_node('beta-alanine')


G.add_edges_from([
    ('beta-alanine', 'Mitochondrial uptake of beta-alanine'),
    ('Mitochondrial uptake of beta-alanine', 'beta-alanine')
    ])

G.add_node('beta-alanine')
G.add_node('R-ALL-113557')

G.add_node('beta-alanine + pyruvate => 3-oxopropanoate + alanine')

G.add_node('R-ALL-909783')
G.add_node('R-ALL-379697')

G.add_edges_from([
    ('beta-alanine', 'beta-alanine + pyruvate => 3-oxopropanoate + alanine'),
    ('R-ALL-113557', 'beta-alanine + pyruvate => 3-oxopropanoate + alanine'),
    ('beta-alanine + pyruvate => 3-oxopropanoate + alanine', 'R-ALL-909783'),
    ('beta-alanine + pyruvate => 3-oxopropanoate + alanine', 'R-ALL-379697')
    ])

G.add_node('inhibitor')
G.add_edge('inhibitor', 'beta-alanine + pyruvate => 3-oxopropanoate + alanine')

G.add_node('catalyst')
G.add_edge('catalyst', 'beta-alanine + pyruvate => 3-oxopropanoate + alanine')

G.add_node('Gene_A')
G.add_node('Gene_B')
G.add_node('Gene_C')
G.add_edge('Gene_A', 'R-ALL-113557', len="3")
G.add_edge('Gene_B', 'R-ALL-113557', len="3")
G.add_edge('Gene_C', 'R-ALL-113557', len="3")

plt.figure(121,figsize=(10,8))
#pos = nx.circular_layout(G)

edges = G.edges()
edge_colors = []
for u,v in edges:

    if u == 'catalyst':
        edge_colors.append('green')
    elif u == 'inhibitor':
        edge_colors.append('red')
    else:
        edge_colors.append('grey')

nodes = G.nodes()
node_colors = []
for u in nodes:

    if '=>' in u:
        node_colors.append('grey')
    elif 'Gene' in u:
        node_colors.append('purple')
    else:
        node_colors.append('lightblue')

nx.draw(
    G,
    #pos, # add if using pos = nx.circular_layout(G)
    nodes=nodes,
    node_color=node_colors,
    edges=edges,
    edge_color=edge_colors,
    with_labels=True,
    font_weight='bold',
    arrowsize=20,
    font_size=8, )
