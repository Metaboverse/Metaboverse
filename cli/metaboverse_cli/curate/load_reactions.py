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
import xml.etree.ElementTree as et
import os
import glob
import tarfile
import json
import re
from utils import progress_feed, track_progress, update_session, safestr

# Global variables
rdf_namespace = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
bqbiol_namespace = '{http://biomodels.net/biology-qualifiers/}'
chebi_split = 'CHEBI:'
uniprot_split = 'uniprot:'
reactome_split = 'reactome:'
gene_split = 'gene='
mirbase_split = 'acc='
other_split = '/'


def get_namespace(sbml_tree):
    """Extract namespace from SBML tree."""
    return re.search('{(.*)}', sbml_tree.tag).group(0)


def handle_folder_contents(directory):
    """Remove all files and the directory itself."""
    if os.path.exists(directory):
        files = glob.glob(os.path.join(directory, "*"))
        for f in files:
            try:
                os.remove(f)
            except:
                print(f'Unable to remove file: {f} ... skipping...')
        try:
            os.rmdir(directory)
        except:
            print(f'Unable to remove directory: {directory} ... skipping...')


def unpack_pathways(output_dir, url='https://reactome.org/download/current/all_species.3.1.sbml.tgz'):
    """Load tarballed SBML reactome pathway files from Reactome site."""
    file = os.path.join(output_dir, os.path.basename(url))
    pathways_dir = file[:-4] + os.path.sep
    handle_folder_contents(pathways_dir)

    os.system(f'curl -kL {url} -o "{file}"')
    os.makedirs(pathways_dir)

    with tarfile.open(file, "r:gz") as tar:
        tar.extractall(path=pathways_dir)
    
    os.remove(file)
    return pathways_dir


def get_pathways(species_id, pathways_dir):
    """Get list of pathways to parse."""
    if not os.path.isdir(pathways_dir):
        raise Exception(f'{pathways_dir} does not exist')

    dir = os.path.abspath(pathways_dir) + os.path.sep
    file_list = os.listdir(dir)
    pathways_list = [f for f in file_list if species_id in f]
    return [f.split('.')[:-1][0] for f in pathways_list]


def get_database(pathways_dir, pathway_name, extension='.sbml'):
    """Import SBML reaction data."""
    pathway_file = os.path.join(pathways_dir, pathway_name + extension)
    pathway_contents = et.parse(pathway_file)
    return pathway_contents.getroot()


def get_metadata(reaction, sbml_namespace):
    """Get basic metadata for a reaction."""
    compartment = safestr(reaction.attrib['compartment'])
    id = safestr(reaction.attrib['id'])
    name = safestr(reaction.attrib['name'])
    reversible = safestr(reaction.attrib['reversible'])
    reversible = 'true' if reversible == 'false' and '<' in name and '>' in name else reversible

    try:
        notes = safestr(reaction.findall(f'{sbml_namespace}notes')[0][0].text)
    except:
        notes = ''
        print(f'No notes available for {name}')

    try:
        reactome = next(
            _rank.attrib[f'{rdf_namespace}resource'].split(reactome_split)[1]
            for rank in reaction.iter(f'{bqbiol_namespace}is')
            for _rank in rank.iter(f'{rdf_namespace}li')
            if 'reactome' in _rank.attrib[f'{rdf_namespace}resource'].lower()
        )
    except:
        reactome = ''
        print(f'No Reactome ID available for {name}')

    return compartment, id, reactome, name, reversible, notes


def add_reaction(pathway_database, reaction, pathway):
    """Add reactions to pathway."""
    _id = reaction.attrib['id']
    pathway_database[pathway]['reactions'].add(_id)
    return pathway_database, _id


def add_reaction_components(type, reaction, sbml_namespace):
    """Add reaction components to reactions database."""
    component_list = reaction.findall(f'{sbml_namespace}{type}')
    if component_list:
        component_list = component_list[0]

    items = [
        (child.attrib.get('species', ''), 'catalyst' if 'catalyst' in child.attrib['id'] or 'positive' in child.attrib['id'] else 'inhibitor' if 'inhibitor' in child.attrib['id'] or 'negative' in child.attrib['id'] else 'other')
        if 'modifier' in child.attrib['id'] else child.attrib.get('species', '')
        for child in component_list
    ]
    return items


def add_reaction_components_manual(type, reaction, sbml_namespace):
    """Add reaction components to reactions database."""
    component_list = reaction.findall(f'{sbml_namespace}{type}')
    if component_list:
        component_list = component_list[0]

    items = [
        (child.attrib.get('species', ''), 'modifier')
        if type == 'listOfModifiers' else child.attrib.get('species', '')
        for child in component_list
    ]
    return items


def add_names(name_database, child, specie, search_string='is'):
    """Add names to dictionary to map species ID."""
    for rank in child.iter(f'{bqbiol_namespace}{search_string}'):
        for _rank in rank.iter(f'{rdf_namespace}li'):
            item = _rank.attrib[f'{rdf_namespace}resource']
            _id = item.split('/')[-1]
            if 'chebi' in item.lower():
                _id = check_chebi(item=_id).split(' ')[0]
            name_database[_id] = specie
            name_database[specie] = specie
            if '(' in _id and ')' in _id:
                name_database = add_alternative_names(name_database, _id, specie)
    return name_database


def add_bigg_names(name_database, child, specie, sbml_namespace, search_string='notes'):
    """Add names to dictionary to map species ID."""
    for c in child:
        if c.tag == f'{sbml_namespace}notes':
            for z in c[0]:
                _z = z.text
                for db in ['bigg', 'biopath', 'kegg', 'metacyc', 'mxnref']:
                    if db in _z.lower():
                        ids = _z.split(':')[1].replace(' ', '').split(',')
                        for _id in ids:
                            name_database[_id] = specie
    name_database[specie] = specie
    return name_database


def add_alternative_names(name_database, item, specie):
    """Add alternative names to name database for mapping."""
    _remove = item[item.find('('): item.find(')') + 1]
    mod_item = item.replace(_remove, '')
    name_database[mod_item] = specie
    return name_database


def check_chebi(item):
    """Special formatting handling for CHEBI entries."""
    item_parsed = item.lower().split(chebi_split.lower())[1]
    return f'CHEBI:{item_parsed}'


def process_item(item):
    """Extract ID and entity type from "is" tag."""
    if 'chebi' in item.lower():
        _id = f'{chebi_split}{item.split(chebi_split)[1]}'
        _type = 'metabolite_component'
    elif 'uniprot' in item.lower():
        _id = item.split(uniprot_split)[1]
        _type = 'protein_component'
    elif 'gene' in item.lower():
        _id = item.split(gene_split)[1]
        _type = 'gene_component'
    elif 'mirbase' in item.lower():
        _id = item.split(mirbase_split)[1]
        _type = 'mirna_component'
    else:
        _id = item.split(other_split)[-1]
        _type = 'other'
    return _id, _type


def add_species(species_database, name_database, compartment_database, components_database, pathway_record, sbml_namespace):
    """Add species records for pathway to database."""
    species = pathway_record.findall(f'{sbml_namespace}listOfSpecies')[0]
    for child in species:
        specie = child.attrib.get('id', '')
        name = child.attrib.get('name', '')
        compartment = child.attrib.get('compartment', '')

        if '[' in name:
            name = name.split(' [')[0]

        species_database[specie] = name
        compartment_database[specie] = compartment
        name_database[name] = specie
        components_database[specie] = {
            'id': specie,
            'reactome_id': '',
            'name': name,
            'is': '',
            'isEncodedBy': '',
            'hasPart': [],
            'type': '',
            'compartment': compartment
        }

        for rank in child.iter(f'{bqbiol_namespace}is'):
            for _rank in rank.iter(f'{rdf_namespace}li'):
                item = _rank.attrib[f'{rdf_namespace}resource']
                _id, _type = process_item(item) if 'reactome' not in item.lower() else (item.split(reactome_split)[1], components_database[specie]['reactome_id'])
                components_database[specie]['is'] = _id
                components_database[specie]['type'] = _type if 'reactome' not in item.lower() else 'complex_component'

        for rank in child.iter(f'{bqbiol_namespace}hasPart'):
            for _rank in rank.iter(f'{rdf_namespace}li'):
                item = _rank.attrib[f'{rdf_namespace}resource']
                _id = item.split(other_split)[-1]
                components_database[specie]['hasPart'].append(_id)

        name_database = add_names(name_database, child, specie, 'is')
        name_database = add_bigg_names(name_database, child, specie, sbml_namespace, 'notes')

    return species_database, name_database, compartment_database, components_database


def process_components(output_dir, pathways_dir, pathways_list, species_id, args_dict=None):
    """Process species-specific pathways."""
    pathway_database = {}
    reaction_database = {}
    species_database = {}
    name_database = {}
    compartment_database = {}
    compartment_dictionary = {}
    components_database = {}

    print(f'Extracting pathway-level reaction data for: {species_id}')

    for counter, pathway in enumerate(pathways_list, 1):
        counter = track_progress(args_dict, counter, len(pathways_list), 7)

        db = get_database(pathways_dir, pathway)
        sbml_namespace = get_namespace(db)

        pathway_record = db.findall(f'{sbml_namespace}model')[0]
        pathway_info = pathway_record.attrib

        id = pathway_info['id']
        pathway_database[pathway] = {
            'id': id,
            'reactome': pathway,
            'name': pathway_info['name'],
            'reactions': set()
        }

        reactions = pathway_record.findall(f'{sbml_namespace}listOfReactions')[0]

        compartments = pathway_record.findall(f'{sbml_namespace}listOfCompartments')[0]
        for c in compartments:
            compartment_dictionary[c.attrib['id']] = c.attrib['name']

        for reaction in reactions:
            compartment, id, reactome, name, reversible, notes = get_metadata(reaction, sbml_namespace)

            pathway_database, reaction_id = add_reaction(pathway_database, reaction, pathway)
            name_database[name] = reaction_id
            reaction_database[id] = {
                'compartment': compartment,
                'id': id,
                'reactome': reactome,
                'name': name,
                'reversible': reversible,
                'notes': notes,
                'reactants': add_reaction_components('listOfReactants', reaction, sbml_namespace),
                'products': add_reaction_components('listOfProducts', reaction, sbml_namespace),
                'modifiers': add_reaction_components('listOfModifiers', reaction, sbml_namespace)
            }

        pathway_database[pathway]['reactions'] = list(pathway_database[pathway]['reactions'])

        species_database, name_database, compartment_database, components_database = add_species(
            species_database, name_database, compartment_database, components_database, pathway_record, sbml_namespace)

    return args_dict, pathway_database, reaction_database, species_database, name_database, compartment_database, compartment_dictionary, components_database


def load_sbml(sbml_url):
    """Load supported SBML file for organism network curation."""
    parsed_url = sbml_url.split(os.path.sep)
    parsed_path = os.path.sep.join(parsed_url[:-1]) + os.path.sep
    parsed_file = parsed_url[-1].split('.')[0]
    parsed_extension = '.' + parsed_url[-1].split('.')[1]

    return get_database(parsed_path, parsed_file, extension=parsed_extension)


def load_custom_json(sbml_url):
    """Load custom SBML file for organism network curation."""
    parsed_url = sbml_url.split(os.path.sep)
    parsed_path = os.path.sep.join(parsed_url[:-1]) + os.path.sep
    parsed_file = parsed_url[-1].split('.')[0]
    parsed_extension = '.' + parsed_url[-1].split('.')[1]

    pathway_file = os.path.join(parsed_path, parsed_file + parsed_extension)
    with open(pathway_file) as json_file:
        data = json.load(json_file)
    return data


def update_model_metadata(sbml_db, args_dict):
    """Get model metadata and update session info."""
    session_file = args_dict['session_data']
    organism_id = sbml_db[0].attrib['id']
    args_dict['organism_id'] = organism_id

    update_session(session_file, 'organism_id', organism_id)
    update_session(session_file, 'organism', sbml_db[0].attrib.get('name', 'unknown'))

    metaid = sbml_db[0].attrib.get('metaid')
    _ver = f'{metaid} ({args_dict["database_source"]})' if metaid else 'N/A'
    update_session(session_file, 'database_version', _ver)
    args_dict['database_version'] = _ver

    return args_dict


def process_manual(sbml_db, args_dict):
    """Parse network curation elements."""
    sbml_namespace = get_namespace(sbml_db)
    pathway_database = {'All': {'id': 'All', 'reactome': 'All', 'name': 'All', 'reactions': set()}}
    reaction_database = {}
    name_database = {}
    compartment_dictionary = {}
    compartment_database = {}
    species_database = {}
    components_database = {}

    args_dict = update_model_metadata(sbml_db, args_dict)
    elements = list(sbml_db[0])

    for x in elements:
        if x.tag == f'{sbml_namespace}listOfCompartments':
            for child in x:
                if child.tag == f'{sbml_namespace}compartment':
                    compartment_dictionary[child.attrib.get('id', '')] = child.attrib.get('name', '')

        if x.tag == f'{sbml_namespace}listOfSpecies':
            for child in x:
                if child.tag == f'{sbml_namespace}species':
                    specie = child.attrib.get('id', '')
                    name = child.attrib.get('name', specie)
                    compartment = child.attrib.get('compartment', '')
                    species_database[specie] = name
                    compartment_database[specie] = compartment
                    name_database[name] = specie
                    components_database[specie] = {
                        'id': specie,
                        'reactome_id': '',
                        'name': name,
                        'is': '',
                        'isEncodedBy': '',
                        'hasPart': [],
                        'type': '',
                        'compartment': compartment
                    }

                    for rank in child.iter(f'{bqbiol_namespace}is'):
                        for _rank in rank.iter(f'{rdf_namespace}li'):
                            item = _rank.attrib[f'{rdf_namespace}resource']
                            _id, _type = process_item(item) if 'reactome' not in item.lower() else (item.split(reactome_split)[1], 'complex_component')
                            components_database[specie]['is'] = _id
                            components_database[specie]['type'] = _type

                    for rank in child.iter(f'{bqbiol_namespace}hasPart'):
                        for _rank in rank.iter(f'{rdf_namespace}li'):
                            item = _rank.attrib[f'{rdf_namespace}resource']
                            components_database[specie]['hasPart'].append(item.split(other_split)[-1])

                    name_database = add_names(name_database, child, specie, 'is')
                    name_database = add_bigg_names(name_database, child, specie, sbml_namespace, 'notes')

    for x in elements:
        if x.tag == f'{sbml_namespace}listOfReactions':
            for child in x:
                if child.tag == f'{sbml_namespace}reaction':
                    _id = child.attrib.get('id', '')
                    _name = child.attrib.get('name', _id)
                    _reversible = child.attrib.get('reversible', 'false')

                    name_database[_name] = _id
                    pathway_database['All']['reactions'].add(_id)
                    reaction_database[_id] = {
                        'compartment': '',
                        'id': _id,
                        'name': _name,
                        'reversible': _reversible,
                        'notes': '',
                        'reactants': add_reaction_components_manual('listOfReactants', child, sbml_namespace),
                        'products': add_reaction_components_manual('listOfProducts', child, sbml_namespace),
                        'modifiers': add_reaction_components_manual('listOfModifiers', child, sbml_namespace)
                    }

    return args_dict, pathway_database, reaction_database, species_database, name_database, compartment_database, compartment_dictionary, components_database


def update_model_metadata_custom(sbml_url, sbml_db, args_dict):
    """Get custom model metadata and update session info."""
    session_file = args_dict['session_data']
    organism_id = os.path.basename(sbml_url).split(".json")[0]
    args_dict['organism_id'] = organism_id
    args_dict['organism'] = organism_id

    update_session(session_file, 'organism_id', organism_id)
    update_session(session_file, 'organism', organism_id)
    update_session(session_file, 'database_version', 'N/A')
    args_dict['database_version'] = 'N/A'

    return args_dict


def add_names_custom(name_database, metabolite_dictionary, name, specie):
    """Add custom names to the name database for mapping."""
    name_database[specie] = specie

    alt_name = name.replace("_", " ").replace("?", "").replace("\u00b1", "").replace("0.001", "").replace("0.999", "").replace("2.0", "").replace("4.0", "").replace("6.0", "")
    if alt_name and alt_name[0] == "-":
        alt_name = alt_name[1:]
    name_database[alt_name] = specie

    if alt_name[-1] == ")":
        alt_name2 = alt_name[alt_name.find('('): alt_name.find(')') + 1]
        name_database[alt_name2] = specie

    if alt_name in metabolite_dictionary:
        name_database[metabolite_dictionary[alt_name]] = specie

    return name_database


def process_custom(sbml_db, sbml_url, args_dict):
    """Parse custom network curation elements."""
    pathway_database = {'All': {'id': 'All', 'reactome': 'All', 'name': 'All', 'reactions': set()}}
    reaction_database = {}
    name_database = {}
    compartment_dictionary = {}
    compartment_database = {}
    species_database = {}
    components_database = {}

    args_dict = update_model_metadata_custom(sbml_url, sbml_db, args_dict)
    compartment_dictionary["N/A"] = "N/A"

    for x in sbml_db["species"]:
        specie = sbml_db["species"][x]["id"]
        name = sbml_db["species"][x]["name"]
        species_type = sbml_db["species"][x]["type"]
        isEncodedBy = ''
        hasPart = []
        if species_type == "modifier":
            isEncodedBy = specie
            hasPart.append(specie)
            species_type = "catalyst"
        else:
            species_type = "metabolite_component"

        species_database[specie] = name
        compartment_database[specie] = "N/A"
        name_database[name] = specie
        components_database[specie] = {
            'id': specie,
            'reactome_id': specie,
            'name': name,
            'is': specie,
            'isEncodedBy': isEncodedBy,
            'hasPart': hasPart,
            'type': species_type,
            'compartment': "N/A"
        }

        name_database = add_names_custom(name_database, sbml_db["synonyms"], name, specie)

    for x in sbml_db["reactions"]:
        _id = x
        _name = sbml_db["reactions"][x]["name"]
        _reversible = "N/A"
        reactants = sbml_db["reactions"][x]['reactants']
        products = sbml_db["reactions"][x]['products']
        modifiers = sbml_db["reactions"][x]['modifiers']

        name_database[_name] = _id
        pathway_database['All']['reactions'].add(_id)
        reaction_database[_id] = {
            'compartment': '',
            'id': _id,
            'name': _name,
            'reversible': _reversible,
            'notes': '',
            'reactants': reactants,
            'products': products,
            'modifiers': modifiers
        }

    return args_dict, pathway_database, reaction_database, species_database, name_database, compartment_database, compartment_dictionary, components_database


def load_reactions(species_id, output_dir, database_source, sbml_url, args_dict):
    """Fetch all reactions for a given organism."""
    if database_source.lower() == 'reactome':
        pathways_dir = unpack_pathways(output_dir)
        progress_feed(args_dict, "graph", 10)

        pathways_list = get_pathways(species_id, pathways_dir)
        progress_feed(args_dict, "graph", 5)

        args_dict, pathway_database, reaction_database, species_database, name_database, compartment_database, compartment_dictionary, components_database = process_components(output_dir, pathways_dir, pathways_list, species_id, args_dict)

        if 'sbml' in pathways_dir:
            handle_folder_contents(pathways_dir)
        else:
            print('Could not find SBML file directory, skipping removal of this directory...')

    elif database_source.lower() == 'biomodels/bigg' and sbml_url != "None":
        sbml_db = load_sbml(sbml_url)
        progress_feed(args_dict, "graph", 10)

        args_dict, pathway_database, reaction_database, species_database, name_database, compartment_database, compartment_dictionary, components_database = process_manual(sbml_db, args_dict)
        progress_feed(args_dict, "graph", 13)

    elif database_source.lower() == 'custom' and sbml_url != "None":
        sbml_db = load_custom_json(sbml_url)
        progress_feed(args_dict, "graph", 10)

        args_dict, pathway_database, reaction_database, species_database, name_database, compartment_database, compartment_dictionary, components_database = process_custom(sbml_db, sbml_url, args_dict)
        progress_feed(args_dict, "graph", 13)

    else:
        raise Exception('Input database type not supported by Metaboverse. If you would like the database type included, please submit an issue at <https://github.com/Metaboverse/Metaboverse/issues>.')

    return args_dict, pathway_database, reaction_database, species_database, name_database, compartment_dictionary, components_database
