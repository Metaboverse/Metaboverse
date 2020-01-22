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

""" 0) Get data
"""
# Import necessary data
output='/Users/jordan/Desktop/'
output_file = output + 'HSA_global_network.json'
species_id = 'HSA'
with open(output + 'HSA_metaboverse_db.pickle', 'rb') as network_file:
    network = pickle.load(network_file)

# Generate output name
if output_file[-5:].lower() == '.json':
    graph_name = output_file
else:
    graph_name = output_file + species_id + '_global_reactions.json'

network.keys()









""" 2) return a table for values with missing ids

nodes:
    species_id: {
        display_name: ___,
        expression_value: [],
        p_value: [],
        rgba_value: [ [], ],
        rjba_js: [ [], ],
        type: ___,
    }


links:
    {target: species_id,
    source: species_id,
    type: ___
    }

pathways:

reactions:

names:

super_pathways:
"""









""" 3) generate graph
"""











""" 4) finalize formatting for viz
"""














""" 5) export
"""
