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
import networkx as nx
import pandas as pd
from datetime import date
import requests
import json
import os
import sys
import pickle
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
    from metaboverse_cli import __version__, SOURCE_URL
    from metaboverse_cli.analyze.utils import remove_defective_reactions
    from metaboverse_cli.utils import progress_feed, track_progress, read_network, \
        get_metaboverse_cli_version, write_database, safestr, update_session_vars
except ImportError:
    try:
        # Then try relative imports
        from __init__ import __version__, SOURCE_URL
        from analyze.utils import remove_defective_reactions
        from utils import progress_feed, track_progress, read_network, \
            get_metaboverse_cli_version, write_database, safestr, update_session_vars
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
            
            init = load_module("__init__", os.path.join(base_path, "..", "__init__.py"))
            __version__ = init.__version__
            SOURCE_URL = init.SOURCE_URL

            utils_analyze = load_module("utils_analyze", os.path.join(base_path, "utils.py"))
            remove_defective_reactions = utils_analyze.remove_defective_reactions

            utils = load_module("utils", os.path.join(base_path, "..", "utils.py"))
            progress_feed = utils.progress_feed
            track_progress = utils.track_progress
            read_network = utils.read_network
            get_metaboverse_cli_version = utils.get_metaboverse_cli_version
            write_database = utils.write_database
            safestr = utils.safestr
            update_session_vars = utils.update_session_vars

        except ImportError as e:
            print(f"Error importing dependencies: {e}")
            print(f"Current sys.path: {sys.path}")
            print(f"Current directory: {os.getcwd()}")
            print(f"File location: {__file__}")
            raise


# Set curation directory
NEIGHBOR_DIR='nbdb'
TEMPLATE_DIR='mvrs'


"""Import internal dependencies
"""
try:
    from analyze.prepare_data import __main__ as prepare_data
    from analyze.model import __template__
    from analyze.model import __model__
    from analyze.model import load_references
    from analyze.model import load_metabolite_synonym_dictionary
    from analyze.utils import remove_defective_reactions
    from utils import progress_feed, track_progress, read_network, \
                      get_metaboverse_cli_version, write_database, safestr, \
                      update_session_vars
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
    prepare_data = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prepare_data)
    prepare_data = prepare_data.__main__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/analyze/model.py"))
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    __template__ = model.__template__
    __model__ = model.__model__
    load_references = model.load_references
    load_metabolite_synonym_dictionary = model.load_metabolite_synonym_dictionary

    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "analyze", "utils.py"
                     ))
    spec = importlib.util.spec_from_file_location("", module_path)
    analyze_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyze_utils)
    remove_defective_reactions = analyze_utils.remove_defective_reactions

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    track_progress = utils.track_progress
    read_network = utils.read_network
    get_metaboverse_cli_version = utils.get_metaboverse_cli_version
    write_database = utils.write_database
    safestr = utils.safestr
    update_session_vars = utils.update_session_vars


def process_data(
        network,
        args_dict):
    """
    """
    if str(args_dict['transcriptomics']).lower() != 'none' \
            or str(args_dict['proteomics']).lower() != 'none' \
            or str(args_dict['metabolomics']).lower() != 'none':

        data, stats, unmapped = prepare_data(
            network=network,
            transcriptomics_url=args_dict['transcriptomics'],
            proteomics_url=args_dict['proteomics'],
            metabolomics_url=args_dict['metabolomics'],
            database_source=args_dict['database_source'])
        flag_data = False

    else:
        data = pd.DataFrame()
        data['NoSample'] = [0, 0, 0]
        data.index = ['dummy_index1', 'dummy_index2', 'dummy_index3']

        stats = pd.DataFrame()
        stats['NoSample'] = [0, 0, 0]
        stats.index = ['dummy_index1', 'dummy_index2', 'dummy_index3']

        unmapped = {
            'transcriptomics_unmapped': [],
            'proteomics_unmapped': []
        }
        flag_data = True

    return data, stats, unmapped, flag_data


def download_template(args_dict, url, user_provided=False):

    def get_template_file(url, args_dict):
        # Construct the file path where the template will be saved
        file_path = os.path.join(
            args_dict['output'],
            args_dict['organism_id'] + '_template.mvrs')
        
        print(f'Downloading Metaboverse graph template for organism from:\n\t{url}')
        
        # Make the request to download the file
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # This will raise an HTTPError if an error occurs
            
            # Open the file path for writing in binary mode
            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192): 
                    file.write(chunk)
                    
            print('\tDownload completed successfully.')
            return file_path
        except requests.exceptions.RequestException as e:
            print(f"\tAn error occurred while downloading the file: {e}")
            return None

    if user_provided == False:
        file = get_template_file(url, args_dict)
    else:
        file = url

    with open(file) as graph_template:
        try:
            return json.load(graph_template), file
        except:
            return None, None


def read_template(
        args_dict,
        graph_data,
        network,
        file):
    """
    """

    print('Parsing Metaboverse graph template for organism...')
    graph = nx.readwrite.json_graph.node_link_graph(
        {
            'nodes': graph_data['nodes'],
            'links': graph_data['links']
        },
        directed=graph_data['directed'],
        multigraph=graph_data['multigraph'])
    network['reaction_database'] = graph_data['reaction_dictionary']
    network['pathway_database'] = graph_data['pathway_dictionary']
    super_pathways = graph_data['super_pathways']
    degree_dictionary = graph_data['degree_dictionary']

    reverse_genes, protein_dictionary, chebi_dictionary, \
    name_reference, uniprot_mapper = load_references(
        args_dict=args_dict,
        ensembl=network['ensembl_synonyms'],
        uniprot=network['uniprot_synonyms'],
        chebi=network['chebi_mapper'],
        uniprot_metabolites=network['uniprot_metabolites'])
    metabolite_mapper = load_metabolite_synonym_dictionary()

    args_dict["curation_version"] = network["metaboverse-curate_version"]
    args_dict["curation_date"] = network["curation_date"]
    args_dict["database_version"] = network["database_version"]
    args_dict['template_url'] = file
    args_dict['template_version'] = graph_data['metadata']['template_version']
    args_dict['template_date'] = graph_data['metadata']['template_date']

    progress_feed(args_dict, "graph", 9)

    return graph, args_dict, network, name_reference, \
        degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper


def download_neighbors_dictionary(
        args_dict,
        url,
        user_provided=False):
    """
    """
    print('Downloading Metaboverse neighbors dictionary for organism...')

    def get_neighbor_file(url, args_dict):
        file = os.path.join(
            args_dict['output'],
            args_dict['organism_id'] + '.nbdb')
        print('Downloading nearest neighbors database...', '\n\t', url)
        os.system('curl -kL ' + url + ' -o \"' + file + '\"')
        return file

    if user_provided == False:
        file = get_neighbor_file(url, args_dict)
    else:
        file = url

    neighbors_dictionary = read_network(
        file_path=args_dict['output'],
        network_url=file)

    neighbors_dictionary['nbdb-Metaboverse-url'] = file

    return neighbors_dictionary


def make_neighbors_dictionary(
        args_dict,
        graph,
        reaction_dictionary):
    """
    """
    reaction_ids = set(reaction_dictionary.keys())

    print('Generating Metaboverse neighbors dictionary for organism...')

    counter = 0
    edges = list(graph.edges)
    edge_len = len(edges)
    neighbors_dictionary = {}
    for e in edges:
        one = e[0]
        two = e[1]
        if one in neighbors_dictionary.keys():
            neighbors_dictionary[one].add(two)
        else:
            neighbors_dictionary[one] = set()
            neighbors_dictionary[one].add(two)
        if two in neighbors_dictionary.keys():
            neighbors_dictionary[two].add(one)
        else:
            neighbors_dictionary[two] = set()
            neighbors_dictionary[two].add(one)
        counter = track_progress(args_dict, counter, edge_len, 3)

    print('Tuning neighbors dictionary...')
    reaction_neighbors_dictionary = {}
    counter = 0
    neighbors = list(neighbors_dictionary.keys())
    neighbors_number = len(neighbors)
    for neighbor in neighbors:
        if neighbor in reaction_ids:
            components = neighbors_dictionary[neighbor]
            connected_reactions = set()
            for _c in components:
                for _c_ in neighbors_dictionary[_c]:
                    if _c_ in reaction_ids:
                        connected_reactions.add(_c_)
            connected_reactions = list(connected_reactions)
            reaction_neighbors_dictionary[neighbor] = connected_reactions
        counter = track_progress(args_dict, counter, neighbors_number, 3)

    reaction_neighbors_dictionary['nbdb-Metaboverse-version'] = get_metaboverse_cli_version()
    reaction_neighbors_dictionary['nbdb-Metaboverse-date'] = date.today().strftime('%Y-%m-%d')
    reaction_neighbors_dictionary['nbdb-Metaboverse-url'] = os.path.join(args_dict['output'], args_dict['organism_id'] + '.nbdb')

    print('Writing neighbors dictionary to database file...')
    write_database(
        output=args_dict['output'],
        file=args_dict['organism_id'] + '.nbdb',
        database=reaction_neighbors_dictionary)

    return reaction_neighbors_dictionary


def __main__(
        args_dict):
    """Analyze data on network model
    """

    # Get network curation info
    network = read_network(
        file_path=args_dict['output'],
        network_url=args_dict['organism_curation_file'])
    progress_feed(args_dict, "graph", 1)

    if args_dict['organism_curation_file'] != 'None':
        args_dict['organism_id'] = network['organism_id']

    # Read in data (if any)
    data, stats, unmapped, flag_data = process_data(
        network=network,
        args_dict=args_dict)
    progress_feed(args_dict, "graph", 2)
    print("Processed data entries: " + str(data.shape[0]))
    
    # Generate graph template
    this_version = get_metaboverse_cli_version()
    test_url = f"{SOURCE_URL}v{this_version}/{args_dict['organism_id']}/{args_dict['organism_id']}_template.mvrs"    
    force_new_curation = args_dict.get('force_new_curation', False) in [True, "true", "True", "TRUE"]
    
    # If unable to access pre-curated network, force new curation
    url_response = None  # Default to None
    if not force_new_curation:
        print(f'Checking for pre-curated network template at:\n\t{test_url}')
        try:
            url_response = requests.head(test_url)
            if url_response.status_code == 404:
                url_response = None
            else:
                print("\tFound pre-curated network template.")
        except:
            print("\tUnable to access source files from: " + str(test_url))
            print("\tWill force a new curation of source files instead...")
            force_new_curation = True

    graph_data, file = None, None
    if not force_new_curation and 'graph_template_file' in args_dict:
        if args_dict['graph_template_file'] not in (None, 'None'):
            graph_data, file = download_template(args_dict, args_dict['graph_template_file'], user_provided=True)
        elif url_response is not None:
            graph_data, file = download_template(args_dict, test_url)

    if graph_data is not None:
        print(f'Download of graph template successful. Parsing...')
        graph, args_dict, network, name_reference, \
        degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper = read_template(
            args_dict=args_dict,
            graph_data=graph_data,
            network=network,
            file=file)
    else:
        print(f'Curating new network model...')
        graph, args_dict, network, name_reference, \
        degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper = __template__(
            args_dict=args_dict,
            network=network,
            species_id=args_dict['organism_id'],
            output_file=args_dict['output_file'])


    if len(graph.nodes) == 0 or len(graph.edges) == 0:
        raise Exception("Unable to generate a reaction-based network based on the input organism template.")
    else:
        print("Successfully loaded network with " + str(len(graph.nodes)) + " nodes and " + str(len(graph.edges)) + " edges")
    
    # Generate graph template
    neighbors_url = f"{SOURCE_URL}v{this_version}/{args_dict['organism_id']}/{args_dict['organism_id']}.nbdb"
    
    # Initial check for pre-curated network, forcing new curation if necessary
    neighbor_response = None
    force_neighbors = False
    if not force_new_curation:
        print(f'Checking for pre-curated neighbors dictionary at:\n\t{neighbors_url}')
        try:
            neighbor_response = requests.head(neighbors_url)
        except Exception as e:
            print(f"\tUnable to access source files from: {neighbors_url}\n\tError: {e}")
            print("\tWill force a new curation of source files instead...")
            force_neighbors = True

    if not force_neighbors:
        if 'neighbor_dictionary_file' in args_dict and safestr(args_dict['neighbor_dictionary_file']) not in (None, 'None'):
            try:
                neighbors_dictionary = download_neighbors_dictionary(args_dict=args_dict, url=args_dict['neighbor_dictionary_file'], user_provided=True)
            except Exception as e:
                print(f"Error downloading user-provided neighbors dictionary: {e}")
                force_neighbors = True

        elif neighbor_response is not None and neighbor_response.status_code != 404:
            try:
                neighbors_dictionary = download_neighbors_dictionary(args_dict=args_dict, url=neighbors_url)
            except Exception as e:
                print(f"Error downloading neighbors dictionary from {neighbors_url}: {e}")
                force_neighbors = True
        else:
            force_neighbors = True

    if force_neighbors:
        no_defective_reactions = remove_defective_reactions(network=network)
        neighbors_dictionary = make_neighbors_dictionary(args_dict=args_dict, graph=graph, reaction_dictionary=no_defective_reactions)
    else:
        progress_feed(args_dict, "graph", 6)
                

    # Overlay data on graph and collapse as able
    print("Modeling data onto network...")
    graph_name = __model__(
        graph=graph,
        args_dict=args_dict,
        network=network,
        data=data,
        stats=stats,
        species_id=args_dict['organism_id'],
        output_file=args_dict['output_file'],
        neighbors_dictionary=neighbors_dictionary,
        name_reference=name_reference,
        degree_dictionary=degree_dictionary,
        chebi_dictionary=chebi_dictionary,
        uniprot_mapper=uniprot_mapper,
        metabolite_mapper=metabolite_mapper,
        super_pathways=super_pathways,
        unmapped=unmapped,
        flag_data=flag_data)

    args_dict = update_session_vars(args_dict)

    return graph_name


def test():
    args_dict = {
        'output': "C:\\Users\\jorda\\Desktop",
        'curation': "MMU.mvdb",
        'organism_id': 'MMU',
        'output_file': "C:\\Users\\jorda\\Desktop\\MMU.mvrs",
        'labels': "0",
        'blocklist': "H+",
        'database_date': "",
        'curation_date': ""}
    args_dict = {
        'output': "C:\\Users\\jorda\\Desktop",
        'curation': "HSA.mvdb",
        'graph_template_file': "C:\\Users\\jorda\\Desktop\\HSA_template.mvrs",
        'organism_id': 'HSA'}

    network = read_network(
        file_path="C:\\Users\\jorda\\Desktop",
        network_url="MMU.mvdb")

    neighbors = read_network(
        file_path="C:\\Users\\jorda\\Desktop",
        network_url="MODEL1604210000.nbdb")

    neighbors_dictionary = download_neighbors_dictionary(
        args_dict=args_dict,
        url="C:\\Users\\jorda\\Desktop\\HSA.nbdb",
        user_provided=True)

    graph, args_dict, network, name_reference, \
    degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper = read_template(
        args_dict=args_dict,
        network=network,
        url=args_dict['graph_template_file'],
        user_provided=True)
