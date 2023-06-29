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


# Set globals
SOURCE_URL='https://rutter.chpc.utah.edu/Metaboverse/source/'
CURATION_DIR='mvdb'


def get_reference(
        args_dict,
        reference_url):
    """Download curation reference
    """

    file = os.path.join(
        args_dict['output'],
        args_dict['organism_id'] + '.mvdb')
    
    print('Downloading pre-curated .MVDB database...', '\n\t', reference_url)
    os.system('curl -kL ' + reference_url + ' -o \"' + file + '\"')

    return file


def main(
        args=None):
    """Run metaboverse-cli
    """

    # Read in arguments
    args, args_dict = parse_arguments(
        args,
        __version__)
    progress_feed(args_dict, "graph", 2)

    # Get info on archived database versions available for direct download
    this_version = get_metaboverse_cli_version()
    reference_url = (
        SOURCE_URL
        + 'v' + this_version + '/'
        + CURATION_DIR + '/'
        + args_dict['organism_id'] + '.mvdb')
    
    # If unable to access pre-curated network, force new curation
    if args_dict['force_new_curation'] != True:
        try:
            url_response = requests.head(reference_url)
        except:
            print("Unable to access source files from: " + str(reference_url))
            print("Will force a new curation of source files instead...")
            args_dict['force_new_curation'] = True
            url_response = ''
    else:
        url_response = ''
        
    if args_dict['cmd'] == 'metaboliteMapper':
        print('Generating metabolite mapper...')
        mapper(args_dict)

    # Run metaboverse-curate
    elif args_dict['cmd'] == 'curate' or args_dict['cmd'] == 'electrum':

        if args_dict['cmd'] == 'curate':
            print('Generating Metaboverse-compatible database...')
        elif args_dict['cmd'] == 'electrum':
            print('Generating Electrum-compatible database...')

        # MVDB file provided by user
        if 'organism_curation_file' in args_dict \
        and safestr(args_dict['organism_curation_file']) != 'None' \
        and safestr(args_dict['organism_curation_file']) != None \
        and safestr(
                args_dict['organism_curation_file']).split('.')[-1] != 'xml' \
        and safestr(
                args_dict['organism_curation_file']).split('.')[-1] != 'sbml' \
        and safestr(
                args_dict['organism_curation_file']).split('.')[-1] != 'json':
            # Update args_dict with path for network model
            args_dict = update_network_vars(args_dict)
            args_dict = update_session_vars(args_dict)
            print('Skipping organism network modeling as one was provided by the user...')
            progress_feed(
                args_dict=args_dict,
                process="graph",
                amount=48)

        # MVDB file exists in repo
        elif (args_dict['force_new_curation'] == False \
        or args_dict['force_new_curation'] == "False") \
        and url_response.status_code != 404 and url_response.status_code != 10054:
            try:
                file = get_reference(
                    args_dict=args_dict,
                    reference_url=reference_url)
                args_dict['organism_curation_file'] = file
                args_dict = update_network_vars(args_dict)
                args_dict = update_session_vars(args_dict)
                print('Skipping organism network modeling as one was found...')
                progress_feed(
                    args_dict=args_dict,
                    process="graph",
                    amount=50)

            except:
                print('Curating network model...')
                args_dict = curate(args_dict)

        # Curate MVDB file from scratch
        else:
            if safestr(
                    args_dict['organism_curation_file']).split('.')[-1] == 'json' \
            and safestr(
                    args_dict['database_source']) == 'custom':
                args_dict['curation'] = args_dict['organism_curation_file']
            
            print('Curating network model...')
            args_dict = curate(args_dict)

        args_dict = init_mvrs_file(args_dict)

        # Curate data overlaid on organism network
        print('Curating data onto the network model...')
        if args_dict['cmd'] == 'curate':
            args_dict['output_file'] = analyze(args_dict)
        elif args_dict['cmd'] == 'electrum':
            curate_target(args_dict)

    # Print some error messaging
    else:
        raise Exception('Invalid sub-module selected')

    args_dict = update_session_vars(args_dict)
    progress_feed(
        args_dict=args_dict,
        process="graph",
        amount=50)


if __name__ == '__main__':
    """Run main
    """
    sys.exit(main() or 0)
