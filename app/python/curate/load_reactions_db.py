"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
    Copyright (C) 2019  Jordan A. Berg
    jordan <dot> berg <at> biochem <dot> utah <dot> edu

    This program is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License along with
    this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""Import dependencies
"""
import os
import sys
import re
import shutil
import time
import hashlib
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from utils import progress_feed

"""Global variables
"""
smbl_namespace = '{{http://www.sbml.org/sbml/level{0}/version{1}/core}}'
rdf_namespace = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
bqbiol_namespace = '{http://biomodels.net/biology-qualifiers/}'
smbl_level = '3'
smbl_version = '1'

"""Functions
"""
def test():

    output_dir = '/Users/jordan/Desktop/'
    pathways_dir = '/Users/jordan/Desktop/metaboverse_data/records/HSA/'
    species_id = 'HSA'
    args_dict = None


def unpack_pathways(
        output_dir,
        url='https://reactome.org/download/current/all_species.3.1.sbml.tgz'):
    """Load tarballed sbml reactome pathway files from reactome site
    """

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    file = output_dir + url.split('/')[-1]

    os.system('curl -L ' + url + ' -o ' + file)

    pathways_dir = file[:-4] + '/'
    if os.path.exists(pathways_dir):
        shutil.rmtree(pathways_dir)
    os.makedirs(pathways_dir)
    os.system('tar -zxvf ' + file + ' -C ' + pathways_dir)
    try:
        os.remove(file)
    except:
        print('Could not find file: ' + file)

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
    dir = os.path.abspath(pathways_dir) + '/'

    # Get list of files and their reaction name
    file_list = os.listdir(dir)
    pathways_list = [f for f in file_list if species_id in f]
    pathways_list = [f.split('.')[:-1][0] for f in pathways_list]

    return pathways_list

def get_database(
        pathways_dir,
        pathway_name):
    """Import sbml reaction data
    """

    if pathways_dir[-1] != '/':
        pathways_dir = pathways_dir + '/'

    pathway_file = pathways_dir + pathway_name + '.sbml'
    pathway_contents = et.parse(pathway_file)
    contents = pathway_contents.getroot()

    return contents

def get_metadata(
        reaction,
        smbl_level,
        smbl_version,
        smbl_namespace=smbl_namespace):
    """Get basic metadata for a reaction
    """

    compartment = reaction.attrib['compartment']
    id = reaction.attrib['id']
    name = reaction.attrib['name']

    reversible = reaction.attrib['reversible']
    if reversible == 'false':
        if '<' in name and '>' in name:
            reversible = 'true'

    try:
        notes = reaction.findall(
            str(smbl_namespace + 'notes').format(
                smbl_level,
                smbl_version
            )
        )[0][0].text
    except:
        notes = ''
        print('No notes available for', name)

    return compartment, id, name, reversible, notes

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
        smbl_namespace=smbl_namespace,
        smbl_level=smbl_level,
        smbl_version=smbl_version):
    """Add reaction components to reactions database
    For type, options are "listOfReactants", "listOfProducts", or
    "listOfModifiers"
    """

    # Collect modifiers for a given reaction by species ID
    component_list = reaction.findall(
        str(smbl_namespace + type).format(
            smbl_level,
            smbl_version
        )
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

            # If element has parentheses, remove what's in between as
            # additional key
            if '(' in _id and ')' in _id:
                name_database = add_alternative_names(
                    name_database=name_database,
                    item=_id,
                    specie=specie)

    return name_database

def add_alternative_names(
        name_database,
        item,
        specie):
    """Add alternative names to name database for mapping
    """

    _remove = item[item.find('(') : item.find(')') + 1]
    mod_item = item.replace(_remove, '')
    name_database[mod_item] = specie

    return name_database

def check_chebi(
        item):
    """Some special formatting handling for CHEBI entries
    """

    item_parsed = item.lower().split('chebi:')[1]
    item_returned = 'CHEBI:' + item_parsed

    return item_returned

def add_species(
        species_database,
        name_database,
        compartment_database,
        components_database,
        pathway_record,
        smbl_namespace=smbl_namespace,
        smbl_level=smbl_level,
        smbl_version=smbl_version,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Add species records for pathway to database
    """

    species = pathway_record.findall(
        str(smbl_namespace + 'listOfSpecies').format(
            smbl_level,
            smbl_version
        )
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
            'hasPart': [],
            'type': '',
            'compartment': compartment
        }

        for rank in child.iter(str(bqbiol_namespace + 'is')):
            for _rank in rank.iter(str(rdf_namespace + 'li')):
                item = _rank.attrib[str(rdf_namespace + 'resource')]
                if 'reactome' not in item.lower():
                    if 'chebi' in item.lower():
                        _id = item.split('chebiId=')[1]
                        components_database[specie]['is'] = _id
                        components_database[specie]['type'] = 'metabolite_component'

                    elif 'uniprot' in item.lower():
                        _id = item.split('/')[-1]
                        components_database[specie]['is'] = _id
                        components_database[specie]['type'] = 'protein_component'

                    elif 'mirbase' in item.lower():
                        _id = item.split('acc=')[1]
                        components_database[specie]['is'] = _id
                        components_database[specie]['type'] = 'mirna_component'

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
                    if 'chebi' in item.lower():
                        _id = item.split('chebiId=')[1]
                        components_database[specie]['hasPart'].append(_id)

                    elif 'uniprot' in item.lower():
                        _id = item.split('/')[-1]
                        components_database[specie]['hasPart'].append(_id)

                    elif 'mirbase' in item.lower():
                        _id = item.split('acc=')[1]
                        components_database[specie]['hasPart'].append(_id)

                    else:
                        pass

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
        smbl_namespace=smbl_namespace,
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

    # Cycle through each pathway database and extract  contents
    for pathway in pathways_list:

        db = get_database(
            pathways_dir,
            pathway)
        smbl_level = db.attrib['level']
        smbl_version = db.attrib['version']

        pathway_record = db.findall(
            str(smbl_namespace + 'model').format(
                smbl_level,
                smbl_version
            )
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
            str(smbl_namespace + 'listOfReactions').format(
                smbl_level,
                smbl_version
            )
        )[0]

        # Parse out compartment IDs and names
        compartments = pathway_record.findall(
            str(smbl_namespace + 'listOfCompartments').format(
                smbl_level,
                smbl_version
            )
        )[0]
        for c in range(len(compartments)):
            id = compartments[c].attrib['id']
            name = compartments[c].attrib['name']
            compartment_dictionary[id] = name

        # Extract reactions from pathway
        for reaction in reactions:

            # Get metadata
            compartment, id, name, reversible, notes = get_metadata(
                reaction=reaction,
                smbl_level=smbl_level,
                smbl_version=smbl_version,
                smbl_namespace=smbl_namespace)

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
                'name': name,
                'reversible': reversible,
                'notes': notes}

            # Collect reactants for a given reaction by species ID
            reaction_database[reaction_id]['reactants'] = add_reaction_components(
                type='listOfReactants',
                reaction=reaction,
                smbl_namespace=smbl_namespace,
                smbl_level=smbl_level,
                smbl_version=smbl_version)

            # Collect products for a given reaction by species ID
            reaction_database[reaction_id]['products'] = add_reaction_components(
                type='listOfProducts',
                reaction=reaction,
                smbl_namespace=smbl_namespace,
                smbl_level=smbl_level,
                smbl_version=smbl_version)

            # Collect modifiers for a given reaction by species ID
            reaction_database[reaction_id]['modifiers'] = add_reaction_components(
                type='listOfModifiers',
                reaction=reaction,
                smbl_namespace=smbl_namespace,
                smbl_level=smbl_level,
                smbl_version=smbl_version)

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
            smbl_namespace=smbl_namespace,
            smbl_level=smbl_level,
            smbl_version=smbl_version,
            bqbiol_namespace=bqbiol_namespace,
            rdf_namespace=rdf_namespace)

    return (pathway_database, reaction_database, species_database,
        name_database, compartment_database, compartment_dictionary,
        components_database)

def __main__(
        species_id,
        output_dir,
        args_dict):
    """Fetch all reactions for a given organism
    """

    #############
    # Make pathway id and reaction ids non R-HSA-etc
    #############

    # Get pathways files
    pathways_dir = unpack_pathways(
        output_dir=output_dir)
    progress_feed(args_dict, "curate", 10)

    pathways_list = get_pathways(
        species_id=species_id,
        pathways_dir=pathways_dir)
    progress_feed(args_dict, "curate", 7)

    # Get list of reaction files to use for populating database
    pathway_database, reaction_database, species_database, \
    name_database, compartment_database, compartment_dictionary, \
    components_database = process_components(
        output_dir=output_dir,
        pathways_dir=pathways_dir,
        pathways_list=pathways_list,
        species_id=species_id,
        args_dict=args_dict)
    progress_feed(args_dict, "curate", 5)

    if 'sbml' in pathways_dir:
        shutil.rmtree(pathways_dir)
    else:
        print('Could not find SMBL file directory, skipping removal of this directory...')

    return (pathway_database, reaction_database, species_database,
        name_database, compartment_dictionary, components_database)
