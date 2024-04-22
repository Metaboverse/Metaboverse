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
# Get source url from local `source_db.txt`
SOURCE_URL = open(
    os.path.join(
        os.path.dirname(
            os.path.abspath(__file__)),
        'source_url.txt'),
    'r').read().strip()
# Set curation directory
CURATION_DIR='mvdb'


"""Import internal dependencies
"""
try:
    from __init__ import __version__
    from arguments import parse_arguments
    from curate.__main__ import __main__ as curate
    from analyze.__main__ import __main__ as analyze
    from mapper.__main__ import __main__ as mapper
    from target.__main__ import __main__ as curate_target
    from utils import progress_feed, update_session, \
        safestr, get_metaboverse_cli_version, init_mvrs_file, \
        update_network_vars, update_session_vars
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath(os.path.join(".", "metaboverse_cli", "__init__.py")))
    version = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(version)
    __version__ = version.__version__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath(os.path.join(".", "metaboverse_cli", "arguments.py")))
    arguments = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(arguments)
    parse_arguments = arguments.parse_arguments

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath(os.path.join(".", "metaboverse_cli", "curate/__main__.py")))
    curate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(curate)
    curate = curate.__main__

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath(os.path.join(".", "metaboverse_cli", "analyze/__main__.py")))
    analyze = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyze)
    analyze = analyze.__main__

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath(os.path.join(".", "metaboverse_cli", "mapper/__main__.py")))
    mapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mapper)
    mapper = mapper.__main__

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath(os.path.join(".", "metaboverse_cli", "target/__main__.py")))
    target = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(target)
    curate_target = target.__main__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath(os.path.join(".", "metaboverse_cli", "utils.py")))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    update_session = utils.update_session
    safestr = utils.safestr
    get_metaboverse_cli_version = utils.get_metaboverse_cli_version
    init_mvrs_file = utils.init_mvrs_file
    update_network_vars = utils.update_network_vars
    update_session_vars = utils.update_session_vars


def get_reference(
        args_dict,
        reference_url):
    """Download curation reference
    """

    file_path = os.path.join(
        args_dict['output'],
        args_dict['organism_id'] + '.mvdb')
    
    print(f'Downloading pre-curated .MVDB database from:\n\t{reference_url}')
    response = requests.get(reference_url, allow_redirects=True)
    if response.status_code == 200:
        with open(file_path, 'wb') as file:
            file.write(response.content)
    else:
        raise Exception(f'Failed to download the file. Status code: {response.status_code}')

    return file_path


def construct_reference_url(args_dict):
    """Construct the reference URL based on arguments."""
    this_version = get_metaboverse_cli_version()
    return f"{SOURCE_URL}v{this_version}/{args_dict['organism_id']}/{args_dict['organism_id']}.mvdb"


def check_url_availability(reference_url):
    """Check if the URL is accessible."""
    try:
        print(f'Checking for pre-curated reference file at: {reference_url}')
        response = requests.head(reference_url)
        if response.status_code == 404:
            print("The reference file does not exist at the provided URL.")
            return False
        return True
    except requests.ConnectionError:
        print("Unable to access source files from:", reference_url)
        return False


def decide_curation_path(args_dict, reference_url):
    """Decide the curation path based on command and file availability.
    """
    
    if args_dict['cmd'] == 'metaboliteMapper':
        print('Generating metabolite mapper...')
        mapper(args_dict)
    
    elif args_dict['cmd'] in ['curate', 'electrum']:
        
        force_new_curation = args_dict.get('force_new_curation', False)

        if force_new_curation:
            print('Forcing new curation of the network model...')
            args_dict = curate(args_dict)

        else:
            if check_url_availability(reference_url):
                # Logic for using pre-curated MVDB file
                try:
                    file = get_reference(args_dict, reference_url)
                    args_dict['organism_curation_file'] = file
                    args_dict['curation'] = file
                    print('Skipping organism network modeling as one was available...')
                    progress_feed(args_dict, "graph", 50)
        
                except Exception as e:
                    print(e)
                    print('Cannot access pre-curated network model. Falling back to generating new model...')
                    args_dict = curate(args_dict)
            
            elif 'organism_curation_file' in args_dict \
            and args_dict['organism_curation_file'] not in (None, 'None'):
                # Logic for provided MVDB file
                args_dict['curation'] = args_dict['organism_curation_file']
                print('Skipping organism network modeling as one was provided by the user...')
                progress_feed(args_dict, "graph", 48)

            else:
                # Curate MVDB file from scratch
                print('Curating network model...')
                args_dict = curate(args_dict)

        
        if args_dict['cmd'] == 'curate':
            print('Starting network data integration and analysis...')
            args_dict['output_file'] = analyze(args_dict)
    
        elif args_dict['cmd'] == 'electrum':
            curate_target(args_dict)
    
    else:
        raise Exception('Invalid sub-module selected')


def main(args=None):
    """Run metaboverse-cli."""
    args, args_dict = parse_arguments(args, __version__)
    progress_feed(args_dict, "graph", 2)

    reference_url = construct_reference_url(args_dict)
    decide_curation_path(args_dict, reference_url)

    # Finalize session
    args_dict = update_session_vars(args_dict)
    progress_feed(args_dict, "graph", 50)

if __name__ == '__main__':
    sys.exit(main() or 0)