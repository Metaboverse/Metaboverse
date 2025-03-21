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
from ftplib import FTP
from contextlib import closing
import pickle
import json
import math
import sys
import os
import zipfile
from pathlib import Path

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

try:
    # First try normal package imports
    from metaboverse_cli.analyze.utils import update_progress
except ImportError:
    try:
        # Then try relative imports
        from analyze.utils import update_progress
    except ImportError:
        try:
            # Finally try direct imports
            import importlib.util
            def load_module(name, path):
                spec = importlib.util.spec_from_file_location(name, path)
                if spec is None:
                    raise ImportError(f"Could not find module at {path}")
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                return module

            base_path = os.path.dirname(__file__)
            analyze_utils = load_module("analyze_utils", os.path.join(base_path, "analyze", "utils.py"))
            update_progress = analyze_utils.update_progress
        except ImportError as e:
            print(f"Error importing dependencies: {e}")
            print(f"Current sys.path: {sys.path}")
            print(f"Current directory: {os.getcwd()}")
            print(f"File location: {__file__}")
            raise

try:
    # First try normal package imports
    from metaboverse_cli import __version__
except ImportError:
    try:
        # Then try relative imports
        from __init__ import __version__
    except ImportError:
        try:
            # Finally try direct imports
            import importlib.util
            init = load_module("__init__", os.path.join(base_path, "__init__.py"))
            __version__ = init.__version__
        except ImportError as e:
            print(f"Error importing version: {e}")
            print(f"Current sys.path: {sys.path}")
            print(f"Current directory: {os.getcwd()}")
            print(f"File location: {__file__}")
            raise


def download_file_ftp(url, output_path):
    """Download a file from an FTP URL to a given output path."""
    # Parse the URL to extract the FTP server address and file path
    if url.startswith('ftp://'):
        url = url[6:]
    server_address, file_path = url.split('/', 1)

    # Use ftplib to connect to the FTP server and download the file
    with closing(FTP(server_address)) as ftp:
        ftp.login()  # Log in as anonymous user
        with open(output_path, 'wb') as file:
            ftp.retrbinary(f'RETR {file_path}', file.write)


def init_mvrs_file(args_dict):

    if args_dict['cmd'] == 'electrum':
        if 'output_file' in args_dict \
                and safestr(args_dict['output_file']) == 'None':
            args_dict['output_file'] = args_dict['output'] \
                + args_dict['organism_id'] \
                + '-latest.eldb'
    else:
        if 'output_file' in args_dict \
                and safestr(args_dict['output_file']) == 'None':
            args_dict['output_file'] = args_dict['output'] \
                + args_dict['organism_id'] \
                + '.mvrs'

    return args_dict


def get_metaboverse_cli_version():
    """Get this version of metaboverse-cli
    """

    return __version__


def update_network_vars(args_dict):
    """Update internal network variables when a pre-curated file is provided
    """

    # check if file exists
    if os.path.isfile(args_dict['organism_curation_file']):
        with open(args_dict['organism_curation_file'], 'rb') as network_file:
            try:
                network = pickle.load(network_file)
                args_dict['organism_id'] = network['organism_id']
                if args_dict['output_file'] == None \
                        or args_dict['output_file'] == "None" \
                        or args_dict['output_file'] == "find":
                    args_dict['output_file'] = args_dict['output'] \
                        + args_dict['organism_id'] \
                        + '.mvrs'
                args_dict['curation'] = args_dict['organism_curation_file'].split(os.path.sep)[-1]
            except:
                print(
                    "Warning: Unable to open organism reference file: " \
                    + args_dict['organism_curation_file'])

    return args_dict


def update_session_vars(args_dict):
    """Update session variables when a pre-curated file is provided
    """

    if 'session_data' in args_dict:
        session_file = args_dict['session_data']
    else:
        return args_dict

    update_session(
        session_file=session_file,
        key='organism_id',
        value=args_dict['organism_id'])
    update_session(
        session_file=session_file,
        key='output_file',
        value=args_dict['output_file'])
    update_session(
        session_file=session_file,
        key='curation',
        value=args_dict['curation'])
    update_session(
        session_file=session_file,
        key='database_url',
        value=args_dict['output_file'])

    return args_dict


def read_network(
        file_path,
        network_url):
    """Read in network from previous curation module
    - was provided as a URL to the file and saved to args_dict['network'] in
    "curate" sub-module
    """

    with open(os.path.join(file_path, network_url), 'rb') as network_file:
        network = pickle.load(network_file)

    return network


def prepare_output(
        output):
    """Get output directory prepared
    """

    # Check provided path exists
    if not os.path.isdir(output):
        os.makedirs(output)

    # Clean up path
    if os.path.abspath(output).endswith(os.path.sep):
        dir = os.path.abspath(output)
    else:
        dir = os.path.abspath(output) + os.path.sep

    return dir


def write_database(
        output,
        file,
        database):
    """Write reactions database to pickle file
    """

    dir = prepare_output(
        output=output)

    # Write information to file
    with open(os.path.join(dir, file), 'wb') as file_product:
        pickle.dump(database, file_product)


def write_database_json(
        output,
        file,
        database):
    """Write reactions database to JSON file
    """

    dir = prepare_output(
        output=output)

    # Write information to file
    with open(dir + file, 'w') as file_product:
        json.dump(database, file_product, indent=4)


def safestr(obj):
    """Covert ascii text if needed
    """

    return str(obj).encode('ascii', 'ignore').decode('ascii')


def update_session(
        session_file,
        key,
        value):
    """Update session information
    """

    if os.path.exists(str(session_file)) and str(session_file) != 'None':

        with open(session_file) as json_file:
            session = json.load(json_file)
            session[key] = value

        with open(session_file, 'w') as outfile:
            json.dump(session, outfile)

    else:
        print("Session file not found: " + str(session_file))


def get_session_value(
        session_file,
        key):

    if os.path.exists(str(session_file)) and str(session_file) != 'None':

        with open(session_file) as json_file:
            session = json.load(json_file)
            if key in session:
                return session[key]
            else:
                print("Cannot find specified key: " + key)
                return 'unknown'

    else:
        print("File at " + str(session_file) + " does not exist.")
        return 'unknown'


def progress_feed(
        args_dict=None,
        process="graph",
        amount=1):
    """JS progress feed
    """

    if args_dict != None:
        if 'progress_log' in args_dict \
                and str(args_dict['progress_log']) != 'None':
            feed_file = args_dict['progress_log']

            if os.path.exists(feed_file) and process != None:

                with open(feed_file) as json_file:
                    data = json.load(json_file)
                    data[process] += amount
                    if data[process] >= 100:
                        data[process] = 100

                with open(feed_file, 'w') as outfile:
                    json.dump(data, outfile)
    else:
        print('Could not access local variables during progress_feed() update.')


def track_progress(
        args_dict,
        _counter,
        _number,
        _total):
    """Keep track of progress of long collapse step
    """

    _counter += 1

    if _counter % max(1, math.floor(_number / _total)) == 0:
        progress = math.floor(_total * (_counter / _number))
        progress_feed(args_dict, "graph", min(1, progress))

    return _counter


def check_directories(
        input,
        argument):
    """Check directory formatting
    """

    # Check that a file wasn't passed in
    if os.path.isdir(input) != True:
        raise Exception(safestr(argument) + ': ' +
                        safestr(input) + ' is not a directory')

    # Check input directory name is formatted correctly and fix if necessary
    input = safestr(os.path.abspath(input))

    if not input.endswith(os.path.sep):
        input += os.path.sep

    return input


def check_files(
        input,
        argument):
    """Check file formatting
    """

    # Check that a file wasn't passed in
    if os.path.isfile(input) != True:
        raise Exception(safestr(argument) + ': ' +
                        safestr(input) + ' is not a file')

    # Check input directory name is formatted correctly and fix if necessary
    input = safestr(os.path.abspath(input))

    return input


def check_curate(
        args_dict):
    """Check curation arguments
    """

    should_exit = False

    if args_dict['organism_id'] == None \
    or args_dict['organism_id'].lower() == 'none' \
    or args_dict['organism_id'].lower() == 'null':
        print('\nIncorrect species identifier provided: ' +
              safestr(args_dict['organism_id']))
        print('Please refer to https://reactome.org/ContentService/data/species/all for a valid list of organisms')
        should_exit = True

    if args_dict['output'] == None \
    or not os.path.exists(os.path.dirname(args_dict['output'])):
        print('\nIncorrect output parameter provided: ' +
              safestr(args_dict['output']))
        should_exit = True

    if 'organism_curation_file' in args_dict \
    and safestr(args_dict['organism_curation_file']) != 'None' \
    and safestr(args_dict['organism_curation_file']) != None:
        if safestr(
                args_dict['organism_curation_file']).split('.')[-1] == 'mvdb':
            pass
        elif safestr(
                args_dict['organism_curation_file']).split('.')[-1] == 'xml' \
        or safestr(
                args_dict['organism_curation_file']).split('.')[-1] == 'sbml':
            pass
        elif safestr(
                args_dict['organism_curation_file']).split('.')[-1] == 'json' \
        and safestr(
                args_dict['database_source']) == 'custom':
            pass
        else:
            print('\nIncorrect organism curation file type provided : ' +
                  safestr(args_dict['organism_curation_file']))
            should_exit = True

    if 'neighbor_dictionary_file' in args_dict \
    and safestr(args_dict['neighbor_dictionary_file']) != 'None' \
    and safestr(args_dict['neighbor_dictionary_file']) != None \
    and safestr(args_dict['neighbor_dictionary_file']).split('.')[-1] != 'nbdb':
        print('\nIncorrect neighbor dictionary file type provided : ' +
              safestr(args_dict['neighbor_dictionary_file']))
        should_exit = True

    if 'graph_template_file' in args_dict \
    and safestr(args_dict['graph_template_file']) != 'None' \
    and safestr(args_dict['graph_template_file']) != None \
    and safestr(args_dict['graph_template_file']).split('.')[-1] != 'mvrs':
        print('\nIncorrect graph template file type provided : ' +
              safestr(args_dict['graph_template_file']))
        should_exit = True
        
    elif 'graph_template_file' in args_dict \
    and safestr(args_dict['graph_template_file']) != 'None' \
    and safestr(args_dict['graph_template_file']) != None \
    and safestr(args_dict['graph_template_file']).split('.')[-1] == 'mvrs':
        should_exit = False

    if should_exit == True:
        sys.exit(1)


def argument_checks(
        args_dict):
    """Run general checks on arguments
    Not sub-module-specific
    """

    # Check output file
    if 'output' in args_dict \
            and args_dict['output'] == None:
        args_dict['output'] = os.getcwd()
    args_dict['output'] = safestr(args_dict['output'])

    if not args_dict['output'].endswith(os.path.sep):
        args_dict['output'] = args_dict['output'] + os.path.sep

    # Check user-provided directory formatting
    for key, value in args_dict.items():

        if key == 'cmd' or key == 'organism_id':
            pass

        elif os.path.isfile(str(value)) == True:
            args_dict[key] = check_files(
                args_dict[key],
                key)

        elif os.path.isdir(str(value)) == True:
            args_dict[key] = check_directories(
                args_dict[key],
                key)

        else:
            pass

    return args_dict
