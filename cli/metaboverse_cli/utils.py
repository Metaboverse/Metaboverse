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
import requests
import pickle
import json
import math
import sys
import os

from __init__ import __version__
SOURCE_URL = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'source_url.txt'), 'r').read().strip()


def download_file_ftp(url, output_path):
    """Download a file from an FTP URL to a given output path."""
    if url.startswith('ftp://'):
        url = url[6:]
    server_address, file_path = url.split('/', 1)

    with closing(FTP(server_address)) as ftp:
        ftp.login()  # Log in as anonymous user
        with open(output_path, 'wb') as file:
            ftp.retrbinary(f'RETR {file_path}', file.write)



def init_mvrs_file(args_dict):
    """Initialize MVRS file based on command."""
    suffix = '-latest.eldb' if args_dict['cmd'] == 'electrum' else '.mvrs'
    if 'output_file' in args_dict and safestr(args_dict['output_file']) == 'None':
        args_dict['output_file'] = os.path.join(args_dict['output'], f"{args_dict['organism_id']}{suffix}")
    return args_dict


def get_metaboverse_cli_version():
    """Get this version of metaboverse-cli."""
    return __version__


def update_network_vars(args_dict):
    """Update internal network variables when a pre-curated file is provided."""
    curation_file = args_dict.get('organism_curation_file')
    if curation_file and os.path.isfile(curation_file):
        with open(curation_file, 'rb') as network_file:
            try:
                network = pickle.load(network_file)
                args_dict['organism_id'] = network['organism_id']
                if not args_dict.get('output_file') or args_dict['output_file'] in ["None", "find"]:
                    args_dict['output_file'] = os.path.join(args_dict['output'], f"{args_dict['organism_id']}.mvrs")
                args_dict['curation'] = os.path.basename(curation_file)
            except Exception as e:
                print(f"Warning: Unable to open organism reference file: {curation_file}")
                print(f"Error: {e}")
    return args_dict


def update_session_vars(args_dict):
    """Update session variables when a pre-curated file is provided."""
    session_file = args_dict.get('session_data')
    for key in ['organism_id', 'output_file', 'curation', 'database_url']:
        update_session(session_file, key, args_dict.get(key))
    return args_dict


def read_network(file_path, network_url):
    """Read in network from previous curation module."""
    with open(os.path.join(file_path, network_url), 'rb') as network_file:
        return pickle.load(network_file)


def prepare_output(output):
    """Get output directory prepared."""
    if not os.path.isdir(output):
        os.makedirs(output)
    return os.path.abspath(output) + os.path.sep


def write_database(output, file, database):
    """Write reactions database to pickle file."""
    dir = prepare_output(output)
    with open(os.path.join(dir, file), 'wb') as file_product:
        pickle.dump(database, file_product)


def write_database_json(output, file, database):
    """Write reactions database to JSON file."""
    dir = prepare_output(output)
    with open(os.path.join(dir, file), 'w') as file_product:
        json.dump(database, file_product, indent=4)


def safestr(obj):
    """Convert to ASCII text if needed."""
    return str(obj).encode('ascii', 'ignore').decode('ascii')


def update_session(session_file, key, value):
    """Update session information."""
    if os.path.exists(session_file) and session_file != 'None':
        with open(session_file) as json_file:
            session = json.load(json_file)
            session[key] = value
        with open(session_file, 'w') as outfile:
            json.dump(session, outfile)
    else:
        print(f"Session file not found: {session_file}")


def get_session_value(session_file, key):
    """Get a session value by key."""
    if os.path.exists(session_file) and session_file != 'None':
        with open(session_file) as json_file:
            session = json.load(json_file)
            return session.get(key, 'unknown')
    print(f"File at {session_file} does not exist.")
    return 'unknown'


def progress_feed(args_dict=None, process="graph", amount=1):
    """JS progress feed."""
    if args_dict:
        feed_file = args_dict.get('progress_log')
        if feed_file and os.path.exists(feed_file) and process:
            with open(feed_file) as json_file:
                data = json.load(json_file)
                data[process] = min(100, data[process] + amount)
            with open(feed_file, 'w') as outfile:
                json.dump(data, outfile)
    else:
        print('Could not access local variables during progress_feed() update.')


def track_progress(args_dict, _counter, _number, _total):
    """Keep track of progress of long collapse step."""
    _counter += 1
    if _counter % max(1, math.floor(_number / _total)) == 0:
        progress = math.floor(_total * (_counter / _number))
        progress_feed(args_dict, "graph", min(1, progress))
    return _counter


def check_directories(input, argument):
    """Check directory formatting."""
    if not os.path.isdir(input):
        raise Exception(f"{safestr(argument)}: {safestr(input)} is not a directory")
    input = safestr(os.path.abspath(input))
    return input if input.endswith(os.path.sep) else input + os.path.sep


def check_files(input, argument):
    """Check file formatting."""
    if not os.path.isfile(input):
        raise Exception(f"{safestr(argument)}: {safestr(input)} is not a file")
    return safestr(os.path.abspath(input))


def check_curate(args_dict):
    """Check curation arguments."""
    should_exit = False

    def print_error(message):
        nonlocal should_exit
        print(message)
        should_exit = True

    organism_id = args_dict.get('organism_id')
    if not organism_id or organism_id.lower() in ['none', 'null']:
        print_error(f"\nIncorrect species identifier provided: {safestr(organism_id)}")
        print_error('Please refer to https://reactome.org/ContentService/data/species/all for a valid list of organisms')

    output = args_dict.get('output')
    if not output or not os.path.exists(os.path.dirname(output)):
        print_error(f"\nIncorrect output parameter provided: {safestr(output)}")

    curation_file = args_dict.get('organism_curation_file')
    if curation_file and safestr(curation_file).split('.')[-1] not in ['mvdb', 'xml', 'sbml', 'json']:
        print_error(f"\nIncorrect organism curation file type provided: {safestr(curation_file)}")

    neighbor_file = args_dict.get('neighbor_dictionary_file')
    if neighbor_file and safestr(neighbor_file).split('.')[-1] != 'nbdb':
        print_error(f"\nIncorrect neighbor dictionary file type provided: {safestr(neighbor_file)}")

    template_file = args_dict.get('graph_template_file')
    if template_file and safestr(template_file).split('.')[-1] != 'mvrs':
        print_error(f"\nIncorrect graph template file type provided: {safestr(template_file)}")

    if should_exit:
        sys.exit(1)


def argument_checks(args_dict):
    """Run general checks on arguments."""
    output = args_dict.get('output', os.getcwd())
    args_dict['output'] = safestr(output if output.endswith(os.path.sep) else output + os.path.sep)

    for key, value in args_dict.items():
        if key in ['cmd', 'organism_id']:
            continue
        if os.path.isfile(value):
            args_dict[key] = check_files(value, key)
        elif os.path.isdir(value):
            args_dict[key] = check_directories(value, key)

    return args_dict


def get_reference(args_dict, reference_url):
    """Download curation reference."""
    file_path = os.path.join(args_dict['output'], f"{args_dict['organism_id']}.mvdb")
    
    print(f"Downloading pre-curated .MVDB database from:\n\t{reference_url}")
    response = requests.get(reference_url, allow_redirects=True)
    if response.status_code == 200:
        with open(file_path, 'wb') as file:
            file.write(response.content)
    else:
        raise Exception(f"Failed to download the file. Status code: {response.status_code}")

    return file_path


def construct_reference_url(args_dict):
    """Construct the reference URL based on arguments."""
    version = get_metaboverse_cli_version()
    return f"{SOURCE_URL}v{version}/{args_dict['organism_id']}/{args_dict['organism_id']}.mvdb"


def check_url_availability(reference_url):
    """Check if the URL is accessible."""
    try:
        print(f"Checking for pre-curated reference file at: {reference_url}")
        response = requests.head(reference_url)
        if response.status_code == 404:
            print("The reference file does not exist at the provided URL.")
            return False
        return True
    except requests.ConnectionError:
        print("Unable to access source files from:", reference_url)
        return False