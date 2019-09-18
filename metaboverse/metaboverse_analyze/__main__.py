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

"""Import internal dependencies
"""
from metaboverse.metaboverse_analyze.curate_data import __main__ as curate_data
from metaboverse.metaboverse_analyze.curate_network import __main__ as curate_network
from metaboverse.metaboverse_analyze.graph import __main__ as graph
from metaboverse.metaboverse_analyze.utils import map_ids
from metaboverse.metaboverse_analyze.utils import retrieve_pathways

"""Analyze data on network model
"""
def __main__(
        args_dict):

    # Read in network
    network = curate_network(
        model=args_dict['model'])

    # Read in data
    data = curate_data(
        metadata=args_dict['metadata'],
        transcriptomics=args_dict['rnaseq'],
        proteomics=args_dict['proteomics'],
        metabolomics=args_dict['metabolomics'],)

    # Map names to work with what metaboverse expects
    data_mapped = map_ids(
        data=data,
        network=network)

    # Get list of pathway to analyze
    pathways = retrieve_pathways(
        args_dict=args_dict,
        network=network)

    # Generate graph(s)
    graph(
        data=data_mapped,
        network=network,
        pathways=pathways,
        species_id=args_dict['species_id'],
        output=args_dict['output'],
        black_list=args_dict['blacklist'])
