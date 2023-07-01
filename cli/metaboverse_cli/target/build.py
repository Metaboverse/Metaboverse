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
from collections import Counter
from networkx.readwrite import json_graph
import networkx as nx
import itertools
import pickle
import json
import numpy as np
import pandas as pd
from datetime import date
import zipfile
import requests
import re
import os

try:
    from analyze.model import build_chebi_reference, build_name_reference, load_metabolite_synonym_dictionary, uniprot_ensembl_reference, gather_synonyms, name_graph, compile_node_degrees
    from utils import progress_feed
except:
    import importlib.util
    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "analyze", "model.py"
                     ))
    spec = importlib.util.spec_from_file_location("", module_path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    build_chebi_reference = model.build_chebi_reference
    build_name_reference = model.build_name_reference
    load_metabolite_synonym_dictionary = model.load_metabolite_synonym_dictionary
    uniprot_ensembl_reference = model.uniprot_ensembl_reference
    gather_synonyms = model.gather_synonyms
    name_graph = model.name_graph
    compile_node_degrees = model.compile_node_degrees

    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "utils.py"
                     ))
    spec = importlib.util.spec_from_file_location("progress_feed", module_path)
    progress_feed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(progress_feed)
    progress_feed = progress_feed.progress_feed


def get_metabolites(
        data,
        columns):
    """Make metabolite dictionary with the available synonyms from the MIDAS
        database
    """

    common_key = columns['Common_metabolite_name']
    common_names = set()
    for index, row in data.iterrows():
        common_names.add(row[common_key])

    return list(common_names)


def fetch_chebi_mappers(
        metabolites=[
            "1,3-Diaminopropane",
            "2-Ketobutyric acid",
            "2-Hydroxybutyric acid",
            "2-Methoxyestrone"],
        metaboanalyst_api="http://api.xialab.ca/mapcompounds"
):
    """Adapted from R from MetaboAnalyst 5.0
    - https://www.metaboanalyst.ca/docs/RTutorial.xhtml#3.2%20Compound%20name%20mapping
    - https://www.metaboanalyst.ca/docs/APIs.xhtml
    """

    input_str = "{\n\t\"queryList\": \""
    for v in metabolites:
        input_str = input_str + v + ";"
    input_str = input_str + "\",\n\t\"inputType\": \"name\"\n}"

    query_results = requests.post(
        metaboanalyst_api,
        data=input_str,
        headers={
            'Content-Type': "application/json",
            'cache-control': "no-cache",
        }
    )

    return query_results.json()


def define_mapper(
        mapper,
        data,
        columns):
    """Make metabolite dictionary with the available synonyms from the MIDAS
        database
    """

    _mapper = {}
    for x in mapper:
        _mapper[x['query']] = x

    name_key = columns['Metabolite']
    common_key = columns['Common_metabolite_name']
    kegg_key = columns['KEGG_ID']
    hmdb_key = columns['HMDB_ID']
    midas_key = columns['MIDAS_ID']
    smiles_key = columns['SMILES_metabolite']

    metabolites = {}
    for index, row in data.iterrows():
        _name = str(row[name_key])
        _common = str(row[common_key])
        _kegg = str(row[kegg_key])
        _hmdb = str(row[hmdb_key])
        _midas = str(row[midas_key])
        _smiles = str(row[smiles_key])

        if _name in metabolites.keys():
            metabolites[_kegg]['hmdb_id'].add(_hmdb)
            metabolites[_kegg]['name'].add(_name)
            metabolites[_kegg]['common_name'].add(_common)
            metabolites[_kegg]['kegg_id'].add(_kegg)
            metabolites[_kegg]['midas_id'].add(_midas)
            metabolites[_kegg]['smiles'].add(_smiles)
        else:
            metabolites[_kegg] = {
                'hmdb_id': set([_hmdb]),
                'name': set([_name]),
                'common_name': set([_common]),
                'kegg_id': set([_kegg]),
                'midas_id': set([_midas]),
                'smiles': set([_smiles]),
                'chebi_id': set()
            }

        if _common in _mapper.keys():
            _hmdb_mod = _mapper[_common]['hmdb_id']
            _kegg_mod = _mapper[_common]['kegg_id']
            _common_mod = _mapper[_common]['hit']
            _chebi = _mapper[_common]['chebi_id']

            metabolites[_kegg]['hmdb_id'].add(_hmdb_mod)
            metabolites[_kegg]['kegg_id'].add(_kegg_mod)
            metabolites[_kegg]['common_name'].add(_common_mod)
            metabolites[_kegg]['chebi_id'].add(_chebi)
            metabolites[_kegg]['chebi_id'].add('CHEBI:' + _chebi)

    return metabolites


def fix_species_ids(
        id_list):
    """Fix specific IDs that are not mapping all synonyms' species IDs
    """

    # pyruvic acid -> pyruvate
    if 'species_5357717' in id_list:
        id_list.add('species_29398')

    return id_list


def targeted_graph(
        metabolites,
        reactions,
        pathways,
        species_reference,
        reversed_species,
        name_database,
        metabolite_mapper,
        uniprot_mapper,
        component_database):
    """Build graph
    - Add nodes and edges
    - Map names to objects in the graph for display
    - Calculate degree for each node
    """

    # Initialize graph object
    reference = {}

    for _m in metabolites.keys():
        # Find species_id
        species_ids, parsed_syns = get_species(
            metabolite=metabolites[_m],
            species_reference=species_reference,
            reversed_species=reversed_species,
            name_database=name_database,
            metabolite_mapper=metabolite_mapper,
            uniprot_mapper=uniprot_mapper)
        species_ids = fix_species_ids(
            id_list=species_ids)

        # Find nearest neighbor reactions
        _reaction_list = set()
        _reactome_list = set()
        for _s in species_ids:
            for _k, _v in reactions.items():
                if _s in _v['reactants'] \
                        or _s in _v['products']:
                    _reaction_list.add(_k)
                    _reactome_list.add(reactions[_k]['reactome'])
                for _mod in _v['modifiers']:
                    if _mod[0] == _s:
                        _reaction_list.add(_k)
                        _reactome_list.add(reactions[_k]['reactome'])

        _pathway_list = set()
        for _p in pathways.keys():
            _reactions = pathways[_p]['reactions']
            _reactome_id = pathways[_p]['reactome']
            for _r in _reaction_list:
                if _r in _reactions:
                    _pathway_list.add(_reactome_id)

        reference[_m] = {
            'id': _m,
            'name': list(metabolites[_m]['name']),
            'common_name': list(metabolites[_m]['common_name']),
            'species_ids': list(species_ids),
            'synonyms': list(parsed_syns),
            'smiles': list(metabolites[_m]['smiles']),
            'pathways': list(_pathway_list),
            'reactions': list(_reaction_list),
            'reactome_reactions': list(_reactome_list)
        }

    return reference


def get_species(
        metabolite,
        species_reference,
        reversed_species,
        name_database,
        metabolite_mapper,
        uniprot_mapper):
    """
    """

    species_ids = set()

    chebi_id = metabolite['chebi_id']
    hmdb_id = metabolite['hmdb_id']
    name = metabolite['name']
    common_name = metabolite['common_name']
    kegg_id = metabolite['kegg_id']
    identifiers = itertools.chain(
        chebi_id,
        hmdb_id,
        name,
        common_name,
        kegg_id)
    try:
        _mapper = next(iter(metabolite['name'])).lower()
        _mapper = metabolite_mapper['mapping_dictionary'][_mapper]
    except:
        try:
            _mapper = next(iter(metabolite['chebi_id']))
            if 'CHEBI:' not in _mapper:
                _mapper = 'CHEBI:' + _mapper
        except:
            _mapper = next(iter(metabolite['name']))

    mapper_id, parsed_syns_list = gather_synonyms(
        map_id=_mapper,
        init_syns=list(identifiers),
        metabolite_mapper=metabolite_mapper,
        uniprot_mapper=uniprot_mapper,
        ignore_enantiomers=True)

    # cross-reference HMDB IDs in relevant databases
    for _id in parsed_syns_list:

        if _id in species_reference:
            species_ids.add(species_reference[_id])
        elif _id in reversed_species:
            species_ids.add(reversed_species[_id])
        elif _id in name_database:
            species_ids.add(name_database[_id])

    return species_ids, parsed_syns_list


def reverse_object(
        data):
    """
    """

    reverse_dictionary = {}
    for k, v in data.items():
        reverse_dictionary[v] = k

    return reverse_dictionary


def test():

    import os
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    read_network = utils.read_network
    network = read_network(
        network_url='C:\\Users\\jorda\\Desktop\\HSA.mvdb')

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/target/utils.py"))
    build_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build_utils)
    import_midas = build_utils.import_midas
    data, columns = import_midas(
        filename='C:\\Users\\jorda\\Desktop\\projects\\Electrum\\_data\\MIDAS-latest.txt')

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/analyze/model.py"))
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    gather_synonyms = model.gather_synonyms

    output_file = 'C:\\Users\\jorda\\Desktop\\HSA-latest.eldb'
    species_id = 'HSA'
    args_dict = {}


def __main__(
        args_dict,
        network,
        data,
        columns,
        species_id,
        output_file):
    """Generate graph object for visualization
    """

    print('Preparing metadata...')
    # Generate output file name
    graph_name = name_graph(
        output_file=output_file,
        species_id=species_id)

    reverse_genes = {v: k for k, v in network['ensembl_synonyms'].items()}
    protein_dictionary = uniprot_ensembl_reference(
        uniprot_reference=network['uniprot_synonyms'],
        ensembl_reference=reverse_genes)

    chebi_dictionary = build_chebi_reference(
        chebi=network['chebi_mapper'],
        uniprot=network['uniprot_metabolites'])

    name_reference = build_name_reference(
        ensembl=network['ensembl_synonyms'],
        uniprot=network['uniprot_synonyms'])

    # add any mapping IDs
    # Add synonyms
    # Change name to user provided if available
    metabolite_mapper = load_metabolite_synonym_dictionary()

    u = {}
    for k, v in network['uniprot_metabolites'].items():
        u[v] = k

    common_metabolites = get_metabolites(
        data=data,
        columns=columns)

    metaboanalyst_chebi_table = fetch_chebi_mappers(
        metabolites=common_metabolites)

    metabolites = define_mapper(
        mapper=metaboanalyst_chebi_table,
        data=data,
        columns=columns)

    reversed_species = reverse_object(
        data=network['species_database'])

    # Generate graph
    # Name mapping
    print('Building network...')
    reference = targeted_graph(
        metabolites=metabolites,
        reactions=network['reaction_database'],
        pathways=network['pathway_database'],
        species_reference=network['species_database'],
        reversed_species=reversed_species,
        name_database=network['name_database'],
        metabolite_mapper=metabolite_mapper,
        uniprot_mapper=u,
        component_database=network['components_database'])

    # May need a reaction id mapper for reactome IDs
    reactome_mapper = {}
    for _k in network['reaction_database'].keys():
        _id = network['reaction_database'][_k]['id']
        _reactome = network['reaction_database'][_k]['reactome']
        reactome_mapper[_id] = _reactome
        reactome_mapper[_reactome] = _id

    no_defective_reactions = {}
    for key in network['reaction_database'].keys():
        rxn_name = network['reaction_database'][key]['name'].lower()
        if 'defective' not in rxn_name \
                and 'mutant' not in rxn_name:
            no_defective_reactions[key] = network['reaction_database'][key]

    no_defective_pathways = {}
    for key in network['pathway_database'].keys():
        rxn_name = network['pathway_database'][key]['name'].lower()
        if 'defective' not in rxn_name \
                and 'mutant' not in rxn_name:
            no_defective_pathways[key] = network['pathway_database'][key]

    output_database = {}
    output_database['neighbors_dictionary'] = reference
    output_database['pathway_dictionary'] = no_defective_pathways
    output_database['reaction_database'] = no_defective_reactions
    output_database['reactome_mapper'] = reactome_mapper
    output_database['curation_date'] = date.today().strftime('%Y-%m-%d')

    with open(graph_name, 'w') as f:
        json.dump(output_database, f, indent=4)

    print('Graphing complete.')

    return graph_name
