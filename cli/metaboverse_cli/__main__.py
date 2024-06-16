"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) Jordan A. Berg, The University of Utah

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
import argparse
import certifi
import sklearn
import networkx
import scipy
import numpy
import pandas
import requests
import pickle
import sys
import os
# Includes all imports used throughout metaboverse-cli to ensure packaging by pyinstaller

# Set globals
CURATION_DIR = 'mvdb'


"""Import internal dependencies
"""
from __init__ import __version__
from arguments import parse_arguments
from curate.__main__ import main as curate
from analyze.__main__ import main as analyze
from mapper.__main__ import main as mapper
from target.__main__ import main as curate_target
from utils import (
    progress_feed, update_session_vars, get_reference, 
    construct_reference_url, check_url_availability
)


def execute_methods(args_dict, reference_url):
    """Decide the curation path based on command and file availability."""
    if args_dict['cmd'] == 'metaboliteMapper':
        print("Generating metabolite mapper...")
        mapper(args_dict)
    elif args_dict['cmd'] in ['curate', 'electrum']:
        force_new_curation = args_dict.get('force_new_curation', False)

        if force_new_curation:
            print("Forcing new curation of the network model...")
            args_dict = curate(args_dict)
        else:
            if check_url_availability(reference_url):
                try:
                    file = get_reference(args_dict, reference_url)
                    args_dict['organism_curation_file'] = file
                    args_dict['curation'] = file
                    print("Skipping organism network modeling as one was available...")
                    progress_feed(args_dict, "graph", 50)
                except Exception as e:
                    print(e)
                    print("Cannot access pre-curated network model. Falling back to generating new model...")
                    args_dict = curate(args_dict)
            elif 'organism_curation_file' in args_dict and args_dict['organism_curation_file'] not in (None, 'None'):
                args_dict['curation'] = args_dict['organism_curation_file']
                print("Skipping organism network modeling as one was provided by the user...")
                progress_feed(args_dict, "graph", 48)
            else:
                print("Curating network model...")
                args_dict = curate(args_dict)

        if args_dict['cmd'] == 'curate':
            print("Starting network data integration and analysis...")
            args_dict['output_file'] = analyze(args_dict)
        elif args_dict['cmd'] == 'electrum':
            curate_target(args_dict)
    else:
        raise Exception("Invalid sub-module selected")


def main(args=None):
    """Run metaboverse-cli."""
    args, args_dict = parse_arguments(args, __version__)
    progress_feed(args_dict, "graph", 2)

    reference_url = construct_reference_url(args_dict)
    execute_methods(args_dict, reference_url)

    # Finalize session
    args_dict = update_session_vars(args_dict)
    progress_feed(args_dict, "graph", 50)

if __name__ == '__main__':
    sys.exit(main() or 0)