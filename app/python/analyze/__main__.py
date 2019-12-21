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
from app.python.analyze.curate_data import __main__ as curate_data
from app.python.analyze.curate_network import __main__ as curate_network
from app.python.analyze.graph import __main__ as graph
from app.python.analyze.utils import map_ids
from app.python.utils import progress_feed

"""Analyze data on network model
"""
def __main__(
        args_dict):

    # Read in network
    progress_feed(args_dict, "graph")
    network = curate_network(
        model=args_dict['model'])
    progress_feed(args_dict, "graph")

    # Read in data
    if args_dict['transcriptomics'].lower() != 'none' \
    or args_dict['proteomics'].lower() != 'none' \
    or args_dict['metabolomics'].lower() != 'none':

        data = curate_data(
            metadata=args_dict['metadata'],
            transcriptomics=args_dict['transcriptomics'],
            proteomics=args_dict['proteomics'],
            metabolomics=args_dict['metabolomics'],
            args_dict=args_dict)
        progress_feed(args_dict, "graph")

        # Map names to work with what metaboverse expects
        data_mapped = map_ids(
            data=data,
            network=network)
        progress_feed(args_dict, "graph")

    else:
        data_mapped = None
        progress_feed(args_dict, "graph")
        progress_feed(args_dict, "graph")

    # Generate graph(s)
    progress_feed(args_dict, "graph")
    graph(
        data=data_mapped,
        network=network,
        species_id=args_dict['species_id'],
        output_file=args_dict['output_file'],
        black_list=args_dict['blacklist'])

    for x in range(10):
        progress_feed(args_dict, "graph")
