"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) 2022 Metaboverse

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
from __future__ import print_function
import pandas as pd
import pickle

"""Import internal dependencies
"""
try:
    from target.build import __main__ as build
    from target.utils import import_midas
    from utils import progress_feed, read_network
except:
    import os
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/target/build.py"))
    build = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build)
    build = build.__main__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/target/utils.py"))
    build_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build_utils)
    import_midas = build_utils.import_midas

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    read_network = utils.read_network


def __main__(
        args_dict):
    """
    INPUTS:
        - Species ID
        - MIDAS database file
        - Output location

    1. Take list of metabolites assayed in MIDAS with IDs, names
    2. Run model.py modules to cross-reference synonyms and find species_ids
    3. Get nodes, find all nearest neighbor reactions (adapt graph.js)
    4. Curate mini-graph with just the metabolites of interest
    5. Add helper dictionary that quickly looks up a metabolite by KEGG ID and
        returns all
    6. For reactions, get list of pathways involved in

    OUTPUT:
        - .eltm file that is the graph and references
    """

    # Get network curation info
    network = read_network(
        network_url=args_dict['network'])
    progress_feed(args_dict, "model", 2)

    if args_dict['organism_curation'] != 'None':
        args_dict['organism_id'] = network['organism_id']

    if str(args_dict['data']).lower() == 'none':
        raise Exception("Missing a valid MIDAS database input.")
    else:
        data, columns = import_midas(
            filename=args_dict['data'])
        progress_feed(args_dict, "model", 3)

        graph_name = build(
            args_dict=args_dict,
            network=network,
            data=data,
            columns=columns,
            species_id=args_dict['organism_id'],
            output_file=args_dict['output_file'])
        progress_feed(args_dict, "model", 10)


def test():

    args_dict = {}
    args_dict['cmd'] = 'electrum'
    args_dict['output'] = 'C:\\Users\\jorda\\Desktop\\'
    output_dir = args_dict['output']
    __main__(args_dict)
