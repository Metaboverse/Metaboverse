"""
MetaboNet-Analyzer
A toolkit for navigating and analyzing gene expression datasets
alias: metabalyze
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
import networkx
import pickle

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
    network_file = args_dict['model']
    network = read_network(
        file=args_dict['model'])

    return network



args_dict = {'model': '/Users/jordan/Desktop/_network/network.pickle'}
network = __main__(args_dict)

p = []
for k, v in network['nodes_reactions'].items():
    p.append(network['nodes_reactions'][k]['processes'])

central_carb = [
    'glycolysis',
    'pyruvate oxidation',
    'citric acid cycle',
    'tca cycle',
    'krebs cycle',
    'pentose phosphate pathway'


    'pyruvate metabolism',

]


p = set(p)
p
