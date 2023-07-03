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
import xml.etree.ElementTree as et
import hashlib
import tarfile
import time
import glob
import stat
import json
import re
import sys
import os


"""Import internal dependencies
"""
try:
    from utils import progress_feed, track_progress, update_session, safestr
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    track_progress = utils.track_progress
    update_session = utils.update_session
    safestr = utils.safestr


"""Global variables
"""
rdf_namespace = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
bqbiol_namespace = '{http://biomodels.net/biology-qualifiers/}'
chebi_split = 'CHEBI:'
uniprot_split = 'uniprot:'
reactome_split = 'reactome:'
gene_split = 'gene='
mirbase_split = 'acc='
other_split = '/'


"""Functions
"""


def test():

    output_dir = '/Users/jordan/Desktop/'
    pathways_dir = '/Users/jordan/Desktop/metaboverse_data/records/HSA/'
    species_id = 'HSA'
    args_dict = None


def get_namespace(
        sbml_tree):
    """
    """

    return re.search('{(.*)}', sbml_tree.tag).group(0)


def handle_folder_contents(dir):

    if os.path.exists(dir):
        try:
            files = glob.glob(dir + "*")
            for f in files:
                os.remove(f)
        except:
            print('Unable to remove files from: ' +
                  str(dir) + ' ... skipping...')

        try:
            os.rmdir(dir)
        except:
            print('Unable to remove directory named: ' +
                  str(dir) + ' ... skipping...')


def unpack_pathways(
        output_dir,
        url='https://reactome.org/download/current/all_species.3.1.sbml.tgz'):
    """Load tarballed sbml reactome pathway files from reactome site
    """

    file = output_dir + url.split('/')[-1]
    pathways_dir = file[:-4] + os.path.sep
    handle_folder_contents(
        dir=pathways_dir)

    os.system('curl -kL ' + url + ' -o \"' + file + '\"')
    os.makedirs(pathways_dir)

    tar = tarfile.open(file, "r:gz")
    tar.extractall(path=pathways_dir)
    tar.close()
    os.remove(file)

    return pathways_dir


def get_pathways(
        species_id,
        pathways_dir):
    """Get list of pathways to parse
    """

    # Check provided path exists
    if not os.path.isdir(pathways_dir):
        raise Exception(pathways_dir, 'does not exist')

    # Clean up path
    if os.path.abspath(pathways_dir).endswith(os.path.sep):
        dir = os.path.abspath(pathways_dir)
    else:
        dir = os.path.abspath(pathways_dir) + os.path.sep

    # Get list of files and their reaction name
    file_list = os.listdir(dir)
    pathways_list = [f for f in file_list if species_id in f]
    pathways_list = [f.split('.')[:-1][0] for f in pathways_list]

    return pathways_list


def get_database(
        pathways_dir,
        pathway_name,
        extension='.sbml'):
    """Import sbml reaction data
    """

    if not pathways_dir.endswith(os.path.sep):
        pathways_dir = pathways_dir + os.path.sep

    pathway_file = pathways_dir + pathway_name + extension
    pathway_contents = et.parse(pathway_file)
    contents = pathway_contents.getroot()

    return contents


def get_metadata(
        reaction,
        sbml_namespace):
    """Get basic metadata for a reaction
    """

    compartment = safestr(reaction.attrib['compartment'])
    id = safestr(reaction.attrib['id'])
    name = safestr(reaction.attrib['name'])

    reversible = safestr(reaction.attrib['reversible'])
    if reversible == 'false':
        if '<' in name and '>' in name:
            reversible = 'true'

    try:
        notes = safestr(reaction.findall(
            str(sbml_namespace + 'notes')
        )[0][0].text)
    except:
        notes = ''
        print('No notes available for', name)

    try:
        for rank in reaction.iter(str(bqbiol_namespace + 'is')):
            for _rank in rank.iter(str(rdf_namespace + 'li')):
                item = _rank.attrib[str(rdf_namespace + 'resource')]
                if 'reactome' in item.lower():
                    reactome = item.split(reactome_split)[1]
    except:
        reactome = ''
        print('No notes available for', name)

    return compartment, id, reactome, name, reversible, notes


def add_reaction(
        pathway_database,
        reaction,
        pathway,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Add reactions to pathway
    """

    _id = reaction.attrib['id']
    pathway_database[pathway]['reactions'].add(_id)

    return pathway_database, _id


def add_reaction_components(
        type,
        reaction,
        sbml_namespace):
    """Add reaction components to reactions database
    For type, options are "listOfReactants", "listOfProducts", or
    "listOfModifiers"
    """

    # Collect modifiers for a given reaction by species ID
    component_list = reaction.findall(
        str(sbml_namespace + type)
    )

    if len(component_list) > 0:
        component_list = component_list[0]

    items = []
    for child in component_list:

        if 'modifier' in child.attrib['id']:

            type = None
            if 'catalyst' in child.attrib['id'] \
                    or 'positive' in child.attrib['id']:
                type = 'catalyst'

            elif 'inhibitor' in child.attrib['id'] \
                    or 'negative' in child.attrib['id']:
                type = 'inhibitor'

            else:
                type = 'other'

            items.append([child.attrib['species'], type])

        else:
            items.append(child.attrib['species'])

    return items


def add_reaction_components_manual(
        type,
        reaction,
        sbml_namespace):
    """Add reaction components to reactions database
    For type, options are "listOfReactants", "listOfProducts", or
    "listOfModifiers"
    """

    # Collect modifiers for a given reaction by species ID
    component_list = reaction.findall(
        str(sbml_namespace + type)
    )

    if len(component_list) > 0:
        component_list = component_list[0]

    items = []
    for child in component_list:

        if type == 'listOfModifiers':
            _type = 'modifier'
            items.append([child.attrib['species'], _type])

        else:
            items.append(child.attrib['species'])

    return items


def add_names(
        name_database,
        child,
        specie,
        search_string='is',
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Add names to dictionary to map species ID
    """

    for rank in child.iter(str(bqbiol_namespace + search_string)):
        for _rank in rank.iter(str(rdf_namespace + 'li')):

            item = _rank.attrib[str(rdf_namespace + 'resource')]
            _id = item.split('/')[-1]
            if 'chebi' in item.lower():
                _id = check_chebi(item=_id)
                _id = _id.split(' ')[0]
            name_database[_id] = specie
            name_database[specie] = specie

            # If element has parentheses, remove what's in between as
            # additional key
            if '(' in _id and ')' in _id:
                name_database = add_alternative_names(
                    name_database=name_database,
                    item=_id,
                    specie=specie)

    return name_database


def add_bigg_names(
        name_database,
        child,
        specie,
        sbml_namespace,
        search_string='notes'):
    """Add names to dictionary to map species ID
    """

    for c in child:
        if c.tag == str(sbml_namespace + 'notes'):
            for z in c[0]:
                _z = z.text
                if 'bigg' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie
                if 'biopath' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie
                if 'kegg' in _z.lower():
                    if ',' in _z:
                        for _z_ in _z.split(','):
                            name_database[_z_.replace(' ', '')] = specie
                    else:
                        name_database[_z.replace(' ', '')] = specie
                if 'metacyc' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie
                if 'mxnref' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie

    name_database[specie] = specie
    return name_database


def add_alternative_names(
        name_database,
        item,
        specie):
    """Add alternative names to name database for mapping
    """

    _remove = item[item.find('('): item.find(')') + 1]
    mod_item = item.replace(_remove, '')
    name_database[mod_item] = specie

    return name_database


def check_chebi(
        item):
    """Some special formatting handling for CHEBI entries
    """

    item_parsed = item.lower().split(chebi_split.lower())[1]
    item_returned = 'CHEBI:' + item_parsed

    return item_returned


def process_item(item):
    """Extract ID and entity type from "is" tag

    Args:
        item (xml tag packet): Extracted "is" tag from xml record

    Returns:
        _id <str>, _type <str>: Extracted ID and entity type
    """
                      
    if 'chebi' in item.lower() and chebi_split.lower() in item.lower():
        _id = str(chebi_split) + str(item.split(chebi_split)[1])
        _type = 'metabolite_component'
    elif 'chebi' in item.lower() and chebi_split.lower() not in item.lower():
        _id = item.split(other_split)[-1]
        if 'chebi' not in _id.lower():
            _id = "CHEBI:" + str(_id)
        _type = 'metabolite_component'
    elif 'uniprot' in item.lower() and uniprot_split.lower() in item.lower():
        _id = item.split(uniprot_split)[1]
        _type = 'protein_component'
    elif 'uniprot' in item.lower() and uniprot_split.lower() not in item.lower():
        _id = item.split(other_split)[-1]
        _type = 'protein_component'
    elif 'gene' in item.lower() and gene_split.lower() in item.lower():
        _id = item.split(gene_split)[1]
        _type = 'gene_component'
    elif 'gene' in item.lower() and gene_split.lower() not in item.lower():
        _id = item.split(other_split)[-1]
        _type = 'gene_component'
    elif 'mirbase' in item.lower() and mirbase_split.lower() in item.lower():
        _id = item.split(mirbase_split)[1]
        _type = 'mirna_component'
    elif 'mirbase' in item.lower() and mirbase_split.lower() not in item.lower():
        _id = item.split(other_split)[-1]
        _type = 'mirna_component'
    else:
        _id = item.split(other_split)[-1]
        _type = 'other'
    
    return _id, _type
                    

def add_species(
        species_database,
        name_database,
        compartment_database,
        components_database,
        pathway_record,
        sbml_namespace,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Add species records for pathway to database
    """

    species = pathway_record.findall(
        str(sbml_namespace + 'listOfSpecies')
    )[0]

    # Collect species information
    for child in species:

        # Initialize specie record and add common name
        specie = child.attrib['id']
        name = child.attrib['name']
        compartment = child.attrib['compartment']

        if '[' in name:
            name = name.split(' [')[0]

        species_database[specie] = name

        # Add compartment membership for specie
        compartment_database[specie] = compartment

        # Add names and ids to name dictionary
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

        for rank in child.iter(str(bqbiol_namespace + 'is')):
            for _rank in rank.iter(str(rdf_namespace + 'li')):
                item = _rank.attrib[str(rdf_namespace + 'resource')]
                if 'reactome' not in item.lower():
                    _id, _type = process_item(item)
                    components_database[specie]['is'] = _id
                    components_database[specie]['type'] = _type
                else:
                    if reactome_split in item.lower():
                        r_id = item.split(reactome_split)[1]
                        components_database[specie]['reactome_id'] = r_id
                    else:
                        r_id = item.split(other_split)[-1]
                        components_database[specie]['reactome_id'] = r_id

        for rank in child.iter(str(bqbiol_namespace + 'hasPart')):
            for _rank in rank.iter(str(rdf_namespace + 'li')):
                item = _rank.attrib[str(rdf_namespace + 'resource')]
                if 'reactome' not in item:
                    components_database[specie]['type'] = 'complex_component'
                    _id, _type = process_item(item)
                    components_database[specie]['hasPart'].append(_id)
                else:
                    _id = item.split(other_split)[-1]
                    components_database[specie]['hasPart'].append(_id)

        # Add source ID
        name_database = add_names(
            name_database=name_database,
            child=child,
            specie=specie,
            search_string='is',
            bqbiol_namespace=bqbiol_namespace,
            rdf_namespace=rdf_namespace)

    return (species_database, name_database, compartment_database,
            components_database)


def process_components(
        output_dir,
        pathways_dir,
        pathways_list,
        species_id,
        args_dict=None,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Process species-specific pathways
    """

    # Initialize databases
    pathway_database = {}
    reaction_database = {}
    species_database = {}
    name_database = {}
    compartment_database = {}
    compartment_dictionary = {}
    components_database = {}

    print('Extracting pathway-level reaction data for: ' + str(species_id))

    counter = 0
    pathway_number = len(pathways_list)

    # Cycle through each pathway database and extract  contents
    for pathway in pathways_list:
        counter = track_progress(args_dict, counter, pathway_number, 7)

        db = get_database(
            pathways_dir,
            pathway)
        sbml_namespace = get_namespace(
            sbml_tree=db
        )

        pathway_record = db.findall(
            str(sbml_namespace + 'model')
        )[0]

        pathway_info = pathway_record.attrib

        id = pathway_info['id']
        pathway_database[pathway] = {
            'id': id,
            'reactome': pathway,
            'name': pathway_info['name'],
            'reactions': set()
        }

        # Parse out reactions
        reactions = pathway_record.findall(
            str(sbml_namespace + 'listOfReactions')
        )[0]

        # Parse out compartment IDs and names
        compartments = pathway_record.findall(
            str(sbml_namespace + 'listOfCompartments')
        )[0]
        for c in range(len(compartments)):
            id = compartments[c].attrib['id']
            name = compartments[c].attrib['name']
            compartment_dictionary[id] = name

        # Extract reactions from pathway
        for reaction in reactions:

            # Get metadata
            compartment, id, reactome, name, reversible, notes = get_metadata(
                reaction=reaction,
                sbml_namespace=sbml_namespace)

            # Get pathway high-level information (reactions, name, compartment)
            pathway_database, reaction_id = add_reaction(
                pathway_database=pathway_database,
                reaction=reaction,
                pathway=pathway,
                bqbiol_namespace=bqbiol_namespace,
                rdf_namespace=rdf_namespace)

            name_database[name] = reaction_id
            reaction_database[id] = {
                'compartment': compartment,
                'id': id,
                'reactome': reactome,
                'name': name,
                'reversible': reversible,
                'notes': notes}

            # Collect reactants for a given reaction by species ID
            reaction_database[reaction_id]['reactants'] = add_reaction_components(
                type='listOfReactants',
                reaction=reaction,
                sbml_namespace=sbml_namespace)

            # Collect products for a given reaction by species ID
            reaction_database[reaction_id]['products'] = add_reaction_components(
                type='listOfProducts',
                reaction=reaction,
                sbml_namespace=sbml_namespace)

            # Collect modifiers for a given reaction by species ID
            reaction_database[reaction_id]['modifiers'] = add_reaction_components(
                type='listOfModifiers',
                reaction=reaction,
                sbml_namespace=sbml_namespace)

        # Convert reaction set for pathway to list
        pathway_database[pathway]['reactions'] = list(
            pathway_database[pathway]['reactions'])

        # Generate species dict
        species_database, name_database, compartment_database, \
            components_database = add_species(
                species_database=species_database,
                name_database=name_database,
                compartment_database=compartment_database,
                components_database=components_database,
                pathway_record=pathway_record,
                sbml_namespace=sbml_namespace,
                bqbiol_namespace=bqbiol_namespace,
                rdf_namespace=rdf_namespace)

    return (args_dict, pathway_database, reaction_database, species_database,
    name_database, compartment_database, compartment_dictionary,
    components_database)


def load_sbml(
        sbml_url):
    """Load supported SBML file for organism network curation
    """

    parsed_url = sbml_url.split(os.path.sep)
    parsed_path = os.path.sep.join(parsed_url[0:-1]) + os.path.sep
    parsed_file = parsed_url[-1].split('.')[0]
    parsed_extension = '.' + parsed_url[-1].split('.')[1]

    return get_database(
        parsed_path,
        parsed_file,
        extension=parsed_extension)

def load_custom_json(
        sbml_url):
    """Load custom SBML file for organism network curation
    """
    
    parsed_url = sbml_url.split(os.path.sep)
    parsed_path = os.path.sep.join(parsed_url[0:-1]) + os.path.sep
    parsed_file = parsed_url[-1].split('.')[0]
    parsed_extension = '.' + parsed_url[-1].split('.')[1]
    
    if not parsed_path.endswith(os.path.sep):
        parsed_path = parsed_path + os.path.sep
    pathway_file = parsed_path + parsed_file + parsed_extension
    
    with open(pathway_file) as json_file:
        data = json.load(json_file)
    
    return data

def update_model_metadata(
        sbml_db,
        args_dict):
    """Get model metadata and update session info
    """

    session_file = args_dict['session_data']
    update_session(
        session_file=session_file,
        key='organism_id',
        value=sbml_db[0].attrib['id'])
    args_dict['organism_id'] = sbml_db[0].attrib['id']

    if 'name' in sbml_db[0].attrib:
        update_session(
            session_file=session_file,
            key='organism',
            value=sbml_db[0].attrib['name'])
    else:
        update_session(
            session_file=session_file,
            key='organism',
            value='unknown')
    if 'metaid' in sbml_db[0].attrib:
        _ver = sbml_db[0].attrib['metaid'] + ' (' + args_dict['database_source'] + ')'
        update_session(
            session_file=session_file,
            key='database_version',
            value=_ver)
        args_dict['database_version'] = _ver
    else:
        update_session(
            session_file=session_file,
            key='database_version',
            value='N/A')
        args_dict['database_version'] = 'N/A'

    return args_dict


def process_manual(
        sbml_db,
        args_dict):
    """Parse network curation elements
    """

    sbml_namespace = get_namespace(
        sbml_tree=sbml_db
    )

    # Initialize databases
    pathway_database = {
        'All': {
            'id': 'All',
            'reactome': 'All',
            'name': 'All',
            'reactions': set()
        }
    }
    reaction_database = {}
    name_database = {}
    compartment_dictionary = {}
    compartment_database = {}
    species_database = {}
    components_database = {}

    # Get model information
    args_dict = update_model_metadata(
        sbml_db=sbml_db,
        args_dict=args_dict
    )

    # Get model categories
    elements = [x for x in sbml_db[0]]

    # Generate compartment dictionary
    for x in elements:
        if x.tag == str(sbml_namespace + 'listOfCompartments'):
            for child in x:
                if child.tag == str(sbml_namespace + 'compartment'):
                    id = child.attrib['id']
                    name = child.attrib['name']
                    compartment_dictionary[id] = name

    # Generate species database
    for x in elements:
        if x.tag == str(sbml_namespace + 'listOfSpecies'):
            for child in x:
                if child.tag == str(sbml_namespace + 'species'):
                    specie = child.attrib['id']
                    if 'name' in child.attrib:
                        name = child.attrib['name']
                    else:
                        name = specie
                    if 'sboTerm' in child.attrib:
                        sboTerm = child.attrib['sboTerm']
                    else:
                        sboTerm = ''
                    compartment = child.attrib['compartment']

                    species_database[specie] = name
                    compartment_database[specie] = compartment
                    name_database[name] = specie
                    components_database[specie] = {
                        'id': specie,
                        'reactome_id': sboTerm,
                        'name': name,
                        'is': specie,
                        'isEncodedBy': '',
                        'hasPart': [],
                        'type': '',
                        'compartment': compartment
                    }

                    for rank in child.iter(str(bqbiol_namespace + 'is')):
                        for _rank in rank.iter(str(rdf_namespace + 'li')):
                            item = _rank.attrib[str(rdf_namespace + 'resource')]
                            if 'reactome' not in item.lower():
                                if 'chebi' in item.lower() \
                                        or 'kegg' in item.lower() \
                                        or 'hmdb' in item.lower() \
                                        or 'bigg' in item.lower():
                                    _id = item.split('/')[-1]
                                    components_database[specie]['is'] = _id
                                    components_database[specie]['type'] = 'metabolite_component'
                                elif 'uniprot' in item.lower():
                                    _id = item.split('/')[-1]
                                    components_database[specie]['is'] = _id
                                    components_database[specie]['type'] = 'protein_component'
                                else:
                                    components_database[specie]['type'] = 'other'
                            else:
                                r_id = item.split('/')[-1]
                                components_database[specie]['reactome_id'] = r_id

                    for rank in child.iter(str(bqbiol_namespace + 'hasPart')):
                        for _rank in rank.iter(str(rdf_namespace + 'li')):
                            item = _rank.attrib[str(rdf_namespace + 'resource')]
                            if 'reactome' not in item:
                                components_database[specie]['type'] = 'complex_component'
                                if 'chebi' in item.lower() \
                                        or 'kegg' in item.lower() \
                                        or 'hmdb' in item.lower() \
                                        or 'bigg' in item.lower():
                                    _id = item.split('/')[-1]
                                    components_database[specie]['hasPart'].append(
                                        _id)
                                elif 'uniprot' in item.lower():
                                    _id = item.split('/')[-1]
                                    components_database[specie]['hasPart'].append(
                                        _id)
                                elif 'mirbase' in item.lower():
                                    _id = item.split('acc=')[1]
                                    components_database[specie]['hasPart'].append(
                                        _id)
                                else:
                                    pass

                    for rank in child.iter(str(bqbiol_namespace + 'isEncodedBy')):
                        for _rank in rank.iter(str(rdf_namespace + 'li')):
                            item = _rank.attrib[str(rdf_namespace + 'resource')]
                            if 'reactome' not in item:
                                if 'kegg.genes' in item.lower():
                                    _id = item.split('/')[-1]
                                    _id_ = _id.split(':')[-1]
                                    components_database[specie]['isEncodedBy'] = _id_
                                else:
                                    _id = item.split('/')[-1]
                                    components_database[specie]['isEncodedBy'] = _id

                    # Add source ID
                    name_database = add_names(
                        name_database=name_database,
                        child=child,
                        specie=specie,
                        search_string='is',
                        bqbiol_namespace=bqbiol_namespace,
                        rdf_namespace=rdf_namespace)

                    name_database = add_bigg_names(
                        name_database=name_database,
                        child=child,
                        specie=specie,
                        sbml_namespace=sbml_namespace,
                        search_string='notes')

    # Generate reaction database
    for x in elements:
        if x.tag == str(sbml_namespace + 'listOfReactions'):
            for child in x:
                if child.tag == str(sbml_namespace + 'reaction'):
                    # Get metadata
                    _id = child.attrib['id']
                    if 'name' in child.attrib:
                        _name = child.attrib['name']
                    else:
                        _name = _id
                    if 'reversible' in child.attrib:
                        _reversible = child.attrib['reversible']
                    else:
                        _reversible = 'false'

                    name_database[_name] = _id
                    pathway_database['All']['reactions'].add(_id)
                    reaction_database[_id] = {
                        'compartment': '',
                        'id': _id,
                        'name': _name,
                        'reversible': _reversible,
                        'notes': ''}
                    reaction_database[_id]['reactants'] = add_reaction_components_manual(
                        type='listOfReactants',
                        reaction=child,
                        sbml_namespace=sbml_namespace)
                    reaction_database[_id]['products'] = add_reaction_components_manual(
                        type='listOfProducts',
                        reaction=child,
                        sbml_namespace=sbml_namespace)
                    reaction_database[_id]['modifiers'] = add_reaction_components_manual(
                        type='listOfModifiers',
                        reaction=child,
                        sbml_namespace=sbml_namespace)

    return (args_dict, pathway_database, reaction_database, species_database,
    name_database, compartment_database, compartment_dictionary,
    components_database)


def update_model_metadata_custom(
        sbml_url,
        sbml_db,
        args_dict):
    """Get custom model metadata and update session info
    """

    session_file = args_dict['session_data']

    args_dict['organism_id'] = args_dict['organism'] = sbml_url.split(os.path.sep)[-1].split(".json")[0]
    update_session(
        session_file=session_file,
        key='organism_id',
        value=args_dict['organism_id'])
    update_session(
        session_file=session_file,
        key='organism',
        value=args_dict['organism'])

    args_dict['database_version'] = 'N/A'
    update_session(
        session_file=session_file,
        key='database_version',
        value=args_dict['database_version'])
    
    return args_dict


def add_names_custom(
        name_database,
        metabolite_dictionary,
        name,
        specie):
    
    name_database[specie] = specie
    
    alt_name = name \
        .replace("_", " ") \
        .replace("?", "") \
        .replace("\u00b1", "") \
        .replace("0.001", "") \
        .replace("0.999", "") \
        .replace("2.0", "") \
        .replace("4.0", "") \
        .replace("6.0", "")
        
    if len(alt_name) > 0:
        if alt_name[0] == "-":
            alt_name = alt_name[1:]
        name_database[alt_name] = specie

        if alt_name[-1] == ")":
            alt_name2 = alt_name[alt_name.find('('): alt_name.find(')') + 1]
            name_database[alt_name2] = specie

    if alt_name in metabolite_dictionary:
        name_database[metabolite_dictionary[alt_name]] = specie
    
    return name_database


def process_custom(
        sbml_db,
        sbml_url,
        args_dict):
    """Parse custom network curation elements
    """

    # Initialize databases
    pathway_database = {
        'All': {
            'id': 'All',
            'reactome': 'All',
            'name': 'All',
            'reactions': set()
        }
    }
    reaction_database = {}
    name_database = {}
    compartment_dictionary = {}
    compartment_database = {}
    species_database = {}
    components_database = {}

    # Get model information
    args_dict = update_model_metadata_custom(
        sbml_url=sbml_url,
        sbml_db=sbml_db,
        args_dict=args_dict
    )

    # Generate compartment dictionary
    compartment_dictionary["N/A"] = "N/A"

    # Generate species database
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
        
        # Add source ID
        name_database = add_names_custom(
            name_database=name_database,
            metabolite_dictionary=sbml_db["synonyms"],
            name=name,
            specie=specie)
        
    # Generate reaction database
    for x in sbml_db["reactions"]:
        _id = x
        _name = sbml_db["reactions"][x]["name"]
        _reversible = "N/A"
        reactants = sbml_db["reactions"][x]['reactants']
        products  = sbml_db["reactions"][x]['products']
        modifiers = sbml_db["reactions"][x]['modifiers']
        
        name_database[_name] = _id
        pathway_database['All']['reactions'].add(_id)
        reaction_database[_id] = {
            'compartment': '',
            'id': _id,
            'name': _name,
            'reversible': _reversible,
            'notes': ''}
        reaction_database[_id]['reactants'] = reactants
        reaction_database[_id]['products']  = products
        reaction_database[_id]['modifiers'] = modifiers
        
    return (args_dict, pathway_database, reaction_database, species_database,
        name_database, compartment_database, compartment_dictionary,
        components_database)
    

def __main__(
        species_id,
        output_dir,
        database_source,
        sbml_url,
        args_dict):
    """Fetch all reactions for a given organism
    """

    # Get pathways files
    if database_source.lower() == 'reactome':
        pathways_dir = unpack_pathways(
            output_dir=output_dir)
        progress_feed(args_dict, "graph", 10)

        pathways_list = get_pathways(
            species_id=species_id,
            pathways_dir=pathways_dir)
        progress_feed(args_dict, "graph", 5)

        # Get list of reaction files to use for populating database
        args_dict, pathway_database, reaction_database, species_database, \
        name_database, compartment_database, compartment_dictionary, \
        components_database = process_components(
            output_dir=output_dir,
            pathways_dir=pathways_dir,
            pathways_list=pathways_list,
            species_id=species_id,
            args_dict=args_dict)

        if 'sbml' in pathways_dir:
            handle_folder_contents(
                dir=pathways_dir)
        else:
            print(
                'Could not find SMBL file directory, skipping removal of this directory...')

    elif database_source.lower() == 'biomodels/bigg' and sbml_url != "None":
        sbml_db = load_sbml(
            sbml_url=sbml_url)
        progress_feed(args_dict, "graph", 10)

        args_dict, pathway_database, reaction_database, species_database, \
        name_database, compartment_database, compartment_dictionary, \
        components_database = process_manual(
                sbml_db=sbml_db,
                args_dict=args_dict)
        progress_feed(args_dict, "graph", 13)

    elif database_source.lower() == 'custom' and sbml_url != "None":
        sbml_db = load_custom_json(
            sbml_url=sbml_url)
        progress_feed(args_dict, "graph", 10)
        
        args_dict, pathway_database, reaction_database, species_database, \
        name_database, compartment_database, compartment_dictionary, \
        components_database = process_custom(
                sbml_db=sbml_db,
                sbml_url=sbml_url,
                args_dict=args_dict)
        progress_feed(args_dict, "graph", 13)
        
    else:
        raise Exception('Input database type not supported by Metaboverse. If you would like the database type included, please submit an issue at <https://github.com/Metaboverse/Metaboverse/issues>.')

    return (args_dict, pathway_database, reaction_database, species_database,
    name_database, compartment_dictionary, components_database)

def test():
    sbml_db = load_sbml(
        sbml_url="C:\\Users\\jorda\\Desktop\\iIS312.xml")
    args_dict = {
        'session_data': '',
        'database_source': ''}

    sbml_db[0].attrib