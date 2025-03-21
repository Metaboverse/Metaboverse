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
import json
import zipfile
from pathlib import Path
# Includes all imports used throughout metaboverse-cli to ensure packaging by pyinstaller

# Set globals
try:
    # First try normal package imports
    from metaboverse_cli import __version__, SOURCE_URL
    from metaboverse_cli.arguments import parse_arguments
    from metaboverse_cli.curate.__main__ import __main__ as curate
    from metaboverse_cli.analyze.__main__ import __main__ as analyze
    from metaboverse_cli.mapper.__main__ import __main__ as metaboliteMapper
    from metaboverse_cli.target.__main__ import __main__ as target
    from metaboverse_cli.utils import (
        update_session_vars,
        progress_feed,
        get_metaboverse_cli_version,
        safestr,
        read_network,
        write_database,
        write_database_json,
        prepare_output
    )
except ImportError:
    try:
        # Then try relative imports
        from __init__ import __version__, SOURCE_URL
        from arguments import parse_arguments
        from curate.__main__ import __main__ as curate
        from analyze.__main__ import __main__ as analyze
        from mapper.__main__ import __main__ as metaboliteMapper
        from target.__main__ import __main__ as target
        from utils import (
            update_session_vars,
            progress_feed,
            get_metaboverse_cli_version,
            safestr,
            read_network,
            write_database,
            write_database_json,
            prepare_output
        )
    except ImportError:
        try:
            # Finally try direct imports with frozen path handling
            import importlib.util
            import os
            import sys

            def get_base_path():
                """Get the base path for the application, handling both frozen and non-frozen cases."""
                if getattr(sys, 'frozen', False):
                    # Running in a PyInstaller bundle
                    return os.path.dirname(sys.executable)
                else:
                    # Running in development
                    return os.path.dirname(os.path.abspath(__file__))

            def load_module(name, path):
                """Load a module from a file path, handling both frozen and non-frozen cases."""
                try:
                    spec = importlib.util.spec_from_file_location(name, path)
                    if spec is None:
                        raise ImportError(f"Could not find module at {path}")
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
                    return module
                except Exception as e:
                    print(f"Error loading module {name} from {path}: {str(e)}")
                    raise

            base_path = get_base_path()
            
            # Load all required modules
            init = load_module("__init__", os.path.join(base_path, "__init__.py"))
            __version__ = init.__version__
            SOURCE_URL = init.SOURCE_URL

            arguments = load_module("arguments", os.path.join(base_path, "arguments.py"))
            parse_arguments = arguments.parse_arguments

            curate_main = load_module("curate_main", os.path.join(base_path, "curate", "__main__.py"))
            curate = curate_main.__main__

            analyze_main = load_module("analyze_main", os.path.join(base_path, "analyze", "__main__.py"))
            analyze = analyze_main.__main__

            mapper_main = load_module("mapper_main", os.path.join(base_path, "mapper", "__main__.py"))
            metaboliteMapper = mapper_main.__main__

            target_main = load_module("target_main", os.path.join(base_path, "target", "__main__.py"))
            target = target_main.__main__

            utils = load_module("utils", os.path.join(base_path, "utils.py"))
            update_session_vars = utils.update_session_vars
            progress_feed = utils.progress_feed
            get_metaboverse_cli_version = utils.get_metaboverse_cli_version
            safestr = utils.safestr
            read_network = utils.read_network
            write_database = utils.write_database
            write_database_json = utils.write_database_json
            prepare_output = utils.prepare_output

        except Exception as e:
            print(f"Error during module loading: {str(e)}")
            raise

# Set curation directory
CURATION_DIR='mvdb'

def get_project_root():
    """Get the path to the project root directory"""
    current_file = Path(__file__).resolve()
    for parent in current_file.parents:
        if parent.name == 'cli':
            return parent
        if parent.name == 'Metaboverse':
            return parent / 'cli'
    return current_file.parent.parent.parent

# Add project root to path
project_root = get_project_root()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
if str(project_root.parent) not in sys.path:
    sys.path.insert(0, str(project_root.parent))

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
    
    if args_dict['cmd'] in ['curate', 'electrum']:
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
            target(args_dict)
    
    else:
        raise Exception('Invalid sub-module selected')


def main(args=None):
    """Run metaboverse-cli."""
    args, args_dict = parse_arguments(args, __version__)
    progress_feed(args_dict, "graph", 2)

    if args_dict['cmd'] == 'metaboliteMapper':
        # For metaboliteMapper, we don't need organism-specific files
        print('Generating metabolite mapper...')
        metaboliteMapper(args_dict)
    else:
        # For other commands that need organism-specific files
        reference_url = construct_reference_url(args_dict)
        decide_curation_path(args_dict, reference_url)

    # Finalize session
    args_dict = update_session_vars(args_dict)
    progress_feed(args_dict, "graph", 50)

if __name__ == '__main__':
    sys.exit(main() or 0)