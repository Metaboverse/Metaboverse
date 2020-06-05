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
import pickle
import pandas as pd

"""Import internal dependencies
"""
from analyze.prepare_data import __main__ as prepare_data
from analyze.model import __main__ as model
from utils import progress_feed

def test():

    args_dict = {
        'network': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE_metaboverse_db.pickle',
        'metabolomics': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/metabolomics_mct1_030min.txt',
        'transcriptomics': 'None',
        'proteomics': 'None',
        'organism_curation': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE_metaboverse_db.pickle',
        'species_id': 'SCE',
        'output_file': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE_global_reactions.json',
    }

    __main__(
        args_dict=args_dict)

def check_db():

    network_url = '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE_metaboverse_db.pickle'
    with open(network_url, 'rb') as network_file:
        network = pickle.load(network_file)

def read_network(
        network_url):
    """Read in network from previous curation module
    - was provided as a URL to the file and saved to args_dict['network'] in
    "curate" sub-module
    """

    with open(network_url, 'rb') as network_file:
        network = pickle.load(network_file)

    return network

def __main__(
        args_dict):
    """Analyze data on network model
    """

    # Get network curation info
    network = read_network(
        network_url=args_dict['network'])
    progress_feed(args_dict, "model", 2)

    if args_dict['organism_curation'] != 'None':
        args_dict['species_id'] = network['species_id']

    # Read in data (if any)
    if str(args_dict['transcriptomics']).lower() != 'none' \
    or str(args_dict['proteomics']).lower() != 'none' \
    or str(args_dict['metabolomics']).lower() != 'none':

        data, stats, unmapped = prepare_data(
            network=network,
            transcriptomics_url=args_dict['transcriptomics'],
            proteomics_url=args_dict['proteomics'],
            metabolomics_url=args_dict['metabolomics'])
        progress_feed(args_dict, "model", 3)
        flag_data = False

    else:
        data = pd.DataFrame()
        data['NoSample'] = [0,0,0]
        data.index = ['dummy_index1', 'dummy_index2', 'dummy_index3']

        stats = pd.DataFrame()
        stats['NoSample'] = [0,0,0]
        stats.index = ['dummy_index1', 'dummy_index2', 'dummy_index3']

        progress_feed(args_dict, "model", 3)
        flag_data = True

    # Generate graph
    graph_name = model(
        args_dict=args_dict,
        network=network,
        data=data,
        stats=stats,
        species_id=args_dict['species_id'],
        output_file=args_dict['output_file'],
        unmapped=unmapped,
        flag_data=flag_data)

    progress_feed(args_dict, "model", 10)
