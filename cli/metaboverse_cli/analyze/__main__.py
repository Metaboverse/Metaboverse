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


# Set globals
# Get source url from local `source_db.txt`
SOURCE_URL = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'source_url.txt'), 'r').read().strip()

# Set curation directory
NEIGHBOR_DIR='nbdb'
TEMPLATE_DIR='mvrs'


# Import internal dependencies
from analyze.prepare_data import prepare_data
from analyze.model import __template__
from analyze.model import __model__
from analyze.model import load_references
from analyze.model import load_metabolite_synonym_dictionary
from analyze.utils import remove_defective_reactions
from utils import (progress_feed, track_progress, read_network, 
                    get_metaboverse_cli_version, write_database, 
                    update_session_vars)


def process_data(network, args_dict):
    if any(str(args_dict.get(key)).lower() != 'none' for key in ['transcriptomics', 'proteomics', 'metabolomics']):
        data, stats, unmapped = prepare_data(
            network=network,
            transcriptomics_url=args_dict['transcriptomics'],
            proteomics_url=args_dict['proteomics'],
            metabolomics_url=args_dict['metabolomics'],
            database_source=args_dict['database_source']
        )
        flag_data = False
    else:
        data = pd.DataFrame(index=['dummy_index1', 'dummy_index2', 'dummy_index3'], data={'NoSample': [0, 0, 0]})
        stats = pd.DataFrame(index=['dummy_index1', 'dummy_index2', 'dummy_index3'], data={'NoSample': [0, 0, 0]})
        unmapped = {'transcriptomics_unmapped': [], 'proteomics_unmapped': []}
        flag_data = True

    return data, stats, unmapped, flag_data


def download_template(args_dict, url, user_provided=False):
    def get_template_file(url, args_dict):
        file_path = os.path.join(args_dict['output'], f"{args_dict['organism_id']}_template.mvrs")
        print(f'Downloading Metaboverse graph template for organism from:\n\t{url}')
        
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192): 
                    file.write(chunk)
                    
            print('\tDownload completed successfully.')
            return file_path
        except requests.exceptions.RequestException as e:
            print(f"\tAn error occurred while downloading the file: {e}")
            return None

    return get_template_file(url, args_dict) if not user_provided else url


def read_template(args_dict, graph_data, network, file):
    print('Parsing Metaboverse graph template for organism...')
    graph = nx.readwrite.json_graph.node_link_graph(
        {'nodes': graph_data['nodes'], 'links': graph_data['links']},
        directed=graph_data['directed'],
        multigraph=graph_data['multigraph']
    )
    network['reaction_database'] = graph_data['reaction_dictionary']
    network['pathway_database'] = graph_data['pathway_dictionary']
    super_pathways = graph_data['super_pathways']
    degree_dictionary = graph_data['degree_dictionary']

    _, _, chebi_dictionary, name_reference, uniprot_mapper = load_references(
        args_dict=args_dict,
        ensembl=network['ensembl_synonyms'],
        uniprot=network['uniprot_synonyms'],
        chebi=network['chebi_mapper'],
        uniprot_metabolites=network['uniprot_metabolites']
    )
    metabolite_mapper = load_metabolite_synonym_dictionary()

    args_dict.update({
        "curation_version": network["metaboverse-curate_version"],
        "curation_date": network["curation_date"],
        "database_version": network["database_version"],
        'template_url': file,
        'template_version': graph_data['metadata']['template_version'],
        'template_date': graph_data['metadata']['template_date']
    })

    progress_feed(args_dict, "graph", 9)

    return (graph, args_dict, network, name_reference, degree_dictionary, super_pathways, 
           chebi_dictionary, uniprot_mapper, metabolite_mapper)


def download_neighbors_dictionary(args_dict, url, user_provided=False):
    def get_neighbor_file(url, args_dict):
        file = os.path.join(args_dict['output'], f"{args_dict['organism_id']}.nbdb")
        print(f'Downloading nearest neighbors database...\n\t{url}')
        os.system(f'curl -kL {url} -o \"{file}\"')
        return file

    return read_network(args_dict['output'], get_neighbor_file(url, args_dict) if not user_provided else url)


def make_neighbors_dictionary(args_dict, graph, reaction_dictionary):
    print('Generating Metaboverse neighbors dictionary for organism...')
    reaction_ids = set(reaction_dictionary.keys())
    neighbors_dictionary = {}

    for counter, e in enumerate(graph.edges, 1):
        neighbors_dictionary.setdefault(e[0], set()).add(e[1])
        neighbors_dictionary.setdefault(e[1], set()).add(e[0])
        track_progress(args_dict, counter, len(graph.edges), 3)

    reaction_neighbors_dictionary = {}
    for counter, neighbor in enumerate(neighbors_dictionary, 1):
        if neighbor in reaction_ids:
            connected_reactions = {nbr for comp in neighbors_dictionary[neighbor] for nbr in neighbors_dictionary[comp] if nbr in reaction_ids}
            reaction_neighbors_dictionary[neighbor] = list(connected_reactions)
        track_progress(args_dict, counter, len(neighbors_dictionary), 3)

    reaction_neighbors_dictionary.update({
        'nbdb-Metaboverse-version': get_metaboverse_cli_version(),
        'nbdb-Metaboverse-date': date.today().strftime('%Y-%m-%d'),
        'nbdb-Metaboverse-url': os.path.join(args_dict['output'], f"{args_dict['organism_id']}.nbdb")
    })

    write_database(args_dict['output'], f"{args_dict['organism_id']}.nbdb", reaction_neighbors_dictionary)
    return reaction_neighbors_dictionary


def main(args_dict):
    network = read_network(args_dict['output'], args_dict['organism_curation_file'])
    progress_feed(args_dict, "graph", 1)

    if args_dict['organism_curation_file'] != 'None':
        args_dict['organism_id'] = network['organism_id']

    data, stats, unmapped, flag_data = process_data(network, args_dict)
    progress_feed(args_dict, "graph", 2)
    print(f"Processed data entries: {data.shape[0]}")
    
    this_version = get_metaboverse_cli_version()
    test_url = f"{SOURCE_URL}v{this_version}/{args_dict['organism_id']}/{args_dict['organism_id']}_template.mvrs"
    force_new_curation = args_dict.get('force_new_curation', False) in [True, "true", "True", "TRUE"]
    
    url_response = None
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
    if not force_new_curation:
        if args_dict.get('graph_template_file') not in (None, 'None'):
            graph_data, file = download_template(args_dict, args_dict['graph_template_file'], user_provided=True)
        elif url_response:
            graph_data, file = download_template(args_dict, test_url)

    if graph_data:
        print(f'Download of graph template successful. Parsing...')
        graph, args_dict, network, name_reference, degree_dictionary, super_pathways, chebi_dictionary, uniprot_mapper, metabolite_mapper = read_template(args_dict, graph_data, network, file)
    else:
        print(f'Curating new network model...')
        graph, args_dict, network, name_reference, degree_dictionary, super_pathways, chebi_dictionary, uniprot_mapper, metabolite_mapper = __template__(args_dict, network, args_dict['organism_id'], args_dict['output_file'])

    if not graph.nodes or not graph.edges:
        raise Exception("Unable to generate a reaction-based network based on the input organism template.")
    else:
        print(f"Successfully loaded network with {len(graph.nodes)} nodes and {len(graph.edges)} edges")
    
    neighbors_url = f"{SOURCE_URL}v{this_version}/{args_dict['organism_id']}/{args_dict['organism_id']}.nbdb"
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
        if args_dict.get('neighbor_dictionary_file') not in (None, 'None'):
            try:
                neighbors_dictionary = download_neighbors_dictionary(args_dict, args_dict['neighbor_dictionary_file'], user_provided=True)
            except Exception as e:
                print(f"Error downloading user-provided neighbors dictionary: {e}")
                force_neighbors = True
        elif neighbor_response and neighbor_response.status_code != 404:
            try:
                neighbors_dictionary = download_neighbors_dictionary(args_dict, neighbors_url)
            except Exception as e:
                print(f"Error downloading neighbors dictionary from {neighbors_url}: {e}")
                force_neighbors = True
        else:
            force_neighbors = True

    if force_neighbors:
        no_defective_reactions = remove_defective_reactions(network)
        neighbors_dictionary = make_neighbors_dictionary(args_dict, graph, no_defective_reactions)
    else:
        progress_feed(args_dict, "graph", 6)

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
        flag_data=flag_data
    )

    args_dict = update_session_vars(args_dict)
    return graph_name

