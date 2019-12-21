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
import shutil
import time
import hashlib
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from app.python.utils import progress_feed

"""Set globals
"""
analyte_prefix = 'R-ALL-'

"""Get list of reactions to parse
"""
def get_reactions(
        species_id,
        reaction_dir):

    # Check provided path exists
    if not os.path.isdir(reaction_dir):
        raise Exception(reaction_dir, 'does not exist')

    # Clean up path
    dir = os.path.abspath(reaction_dir) + '/'

    # Get list of files and their reaction name
    file_list = os.listdir(dir)
    reactions_list = [f for f in file_list if species_id in f]
    reactions_list = [f.split('.')[:-1][0] for f in reactions_list]

    return reactions_list

"""Get pathway names
"""
def add_pathways(
        database):

    pathways = {}

    # Cycle through processes
    for key in database.keys():

        if database[key]['pathway_name'] in pathways.keys():
            pathways[database[key]['pathway_name']].append(database[key]['reactome_id'])

        else:
            pathways[database[key]['pathway_name']] = [database[key]['reactome_id']]

    # Check for duplicate processes
    for key, value in pathways.items():

        if len(value) > 1:

            reactions = []

            for v in value:

                reactions.append(list(database[v]['reactions'].keys()))

            # Convert lists of reactions to hashes for comparison
            hashes = []

            for r in reactions:

                text_hashed = hashlib.sha256(
                    bytes(
                        str(r),
                        encoding='utf8')
                        ).hexdigest()
                hashes.append(text_hashed)

            # If a list of lists contain identical reactions, keep the first reaction name
            if len(set(hashes)) == 1:
                pathways[key] = [value[0]]

    return pathways

"""Get compartments
"""
def add_compartments(
        database):

    compartments = []

    # Cycle through processes
    for key_x in database.keys():

        # Cycle through reactions
        for key_y in database[key_x]['reactions'].keys():

            compartments.append(database[key_x]['reactions'][key_y]['compartment'])

    return set(compartments)

"""Curate database for all reactions and processes for an organism
"""
def curate_reactions(
        reaction_dir,
        reactions_list,
        species_id,
        species_key='species',
        args_dict=None):

    species_tag = 'R-' + species_id + '-'

    reactions = {}

    counter = 1
    total = len(reactions_list)

    for reaction in reactions_list:

        database = get_database(
            reaction_dir,
            reaction)

        reactions[reaction] = {}

        for child in database:

            if 'model' in child.tag:
                reactions[reaction] = get_process_information(
                    child=child,
                    reactome_id=reaction,
                    reaction_dict=reactions[reaction])

                reactions[reaction]['reactions'] = {}

                for child_2 in child:

                    if 'listOfReactions' in child_2.tag:

                        for child_3 in child_2:

                            reaction_name = species_tag + child_3.attrib['id'].split('_')[1]
                            reactions[reaction]['reactions'] = get_reaction_information(
                                child=child_3,
                                reaction_name=reaction_name,
                                reaction_dict=reactions[reaction]['reactions'],
                                species_tag=species_tag)

                            reactions[reaction]['reactions'][reaction_name]['reactants'] = {}
                            reactions[reaction]['reactions'][reaction_name]['products'] = {}
                            reactions[reaction]['reactions'][reaction_name]['modifiers'] = {}

                            for child_4 in child_3:


                                if 'listOfReactants' in child_4.tag:

                                    for child_5 in child_4:

                                        reactant_name = analyte_prefix + child_5.attrib[species_key].split('_')[1]
                                        reactions[reaction]['reactions'][reaction_name]['reactants'] = populate_reactants(
                                            child=child_5,
                                            species_name=reactant_name,
                                            reactant_dict=reactions[reaction]['reactions'][reaction_name]['reactants'])

                                elif 'listOfProducts' in child_4.tag:

                                    for child_6 in child_4:

                                        product_name = analyte_prefix + child_6.attrib[species_key].split('_')[1]
                                        reactions[reaction]['reactions'][reaction_name]['products'] = populate_products(
                                            child=child_6,
                                            species_name=product_name,
                                            products_dict=reactions[reaction]['reactions'][reaction_name]['products'])

                                elif 'listOfModifiers' in child_4.tag:

                                    for child_7 in child_4:

                                        modifier_name = analyte_prefix + child_7.attrib[species_key].split('_')[1]
                                        reactions[reaction]['reactions'][reaction_name]['modifiers'] = populate_modifiers(
                                            child=child_7,
                                            species_name=modifier_name,
                                            modifiers_dict=reactions[reaction]['reactions'][reaction_name]['modifiers'])

                                else:
                                    pass

        if int(counter % (total / 5)) == 0 and args_dict != None:
            progress_feed(args_dict, "reactions")

        counter += 1

    return reactions

"""Import sbml reaction data
"""
def get_database(
        reaction_dir,
        reaction_name):

    if reaction_dir[-1] != '/':
        reaction_dir = reaction_dir + '/'

    reactome_file = reaction_dir + reaction_name + '.sbml'
    reactome = et.parse(reactome_file)
    root = reactome.getroot()

    return root

"""Add overall process information
"""
def get_process_information(
        child,
        reactome_id,
        reaction_dict,
        name_key='name'):

    reaction_dict['reactome_id'] = reactome_id
    reaction_dict['pathway_name'] = child.attrib[name_key]

    return reaction_dict

"""Add general information for reaction
"""
def get_reaction_information(
        child,
        reaction_dict,
        reaction_name,
        species_tag,
        species_key='species',
        species_splitter='_',
        species_position=1,
        id_key='id',
        id_splitter='_',
        id_position=1,
        name_key='name',
        reversible_key='reversible',
        fast_key='fast',
        compartment_key='compartment'):

    reaction_dict[reaction_name] = {}
    reaction_dict[reaction_name]['id'] = species_tag + child.attrib[id_key].split(id_splitter)[id_position]
    reaction_dict[reaction_name]['name'] = child.attrib[name_key]
    reaction_dict[reaction_name]['reversible'] = child.attrib[reversible_key]
    reaction_dict[reaction_name]['fast'] = child.attrib[fast_key]
    reaction_dict[reaction_name]['compartment'] = species_tag + child.attrib[compartment_key].split(species_splitter)[species_position]

    return reaction_dict

"""Add reactant information for reaction
"""
def populate_reactants(
        child,
        species_name,
        reactant_dict,
        species_key='species',
        species_splitter='_',
        species_position=1,
        id_key='id',
        id_splitter='_',
        id_position=2,
        constant_key='constant',
        stoichiometry_key='stoichiometry'):

    reactant_dict[species_name] = {}
    reactant_dict[species_name]['species_id'] = analyte_prefix + child.attrib[species_key].split(species_splitter)[species_position]
    reactant_dict[species_name]['type'] = child.attrib[id_key].split(id_splitter)[id_position]
    reactant_dict[species_name]['constant'] = child.attrib[constant_key]
    reactant_dict[species_name]['stoichiometry'] = child.attrib[stoichiometry_key]

    return reactant_dict

"""Add product information for reaction
"""
def populate_products(
        child,
        species_name,
        products_dict,
        species_key='species',
        species_splitter='_',
        species_position=1,
        id_key='id',
        id_splitter='_',
        id_position=2,
        constant_key='constant',
        stoichiometry_key='stoichiometry'):

    products_dict[species_name] = {}
    products_dict[species_name]['species_id'] = analyte_prefix + child.attrib[species_key].split(species_splitter)[species_position]
    products_dict[species_name]['type'] = child.attrib[id_key].split(id_splitter)[id_position]
    products_dict[species_name]['constant'] = child.attrib[constant_key]
    products_dict[species_name]['stoichiometry'] = child.attrib[stoichiometry_key]

    return products_dict

"""Add modifier information for reaction
"""
def populate_modifiers(
        child,
        species_name,
        modifiers_dict,
        species_key='species',
        species_splitter='_',
        species_position=1,
        id_key='id',
        id_splitter='_',
        id_position=2):

    modifiers_dict[species_name] = {}
    modifiers_dict[species_name]['species_id'] = analyte_prefix + child.attrib[species_key].split(species_splitter)[species_position]
    modifiers_dict[species_name]['type'] = child.attrib[id_key].split(id_splitter)[id_position]
    modifiers_dict[species_name]['constant'] = None
    modifiers_dict[species_name]['stoichiometry'] = None

    return modifiers_dict

"""Load tarballed sbml reactome reaction files from reactome site
"""
def unpack_reactions(
        output_dir,
        url='https://reactome.org/download/current/all_species.3.1.sbml.tgz'):

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    file = output_dir + url.split('/')[-1]

    os.system('curl -L ' + url + ' -o ' + file)

    reactions_dir = file[:-4] + '/'
    if os.path.exists(reactions_dir):
        shutil.rmtree(reactions_dir)
    os.makedirs(reactions_dir)
    os.system('tar -zxvf ' + file + ' -C ' + reactions_dir)
    try:
        os.remove(file)
    except:
        print('Could not find file: ' + file)

    return reactions_dir

"""Extract global reactions
"""
def add_global_reactions(
        pathways):

    global_reactions = {}

    # Parse through each pathway to extract its reactions
    for pathway in pathways.keys():

        pathway_id = pathways[pathway]['reactome_id']
        pathway_name = pathways[pathway]['pathway_name']

        # Parse out reactions to global reaction dictionary, keep same structure, but add pathway membership
        for reaction in pathways[pathway]['reactions'].keys():

            global_reactions[reaction] = {}
            global_reactions[reaction]['pathway_id'] = pathway_id
            global_reactions[reaction]['pathway_name'] = pathway_name
            global_reactions[reaction]['id'] = pathways[pathway]['reactions'][reaction]['id']
            global_reactions[reaction]['name'] = pathways[pathway]['reactions'][reaction]['name']
            global_reactions[reaction]['reversible'] = pathways[pathway]['reactions'][reaction]['reversible']
            global_reactions[reaction]['fast'] = pathways[pathway]['reactions'][reaction]['fast']
            global_reactions[reaction]['compartment'] = pathways[pathway]['reactions'][reaction]['compartment']
            global_reactions[reaction]['reactants'] = pathways[pathway]['reactions'][reaction]['reactants']
            global_reactions[reaction]['products'] = pathways[pathway]['reactions'][reaction]['products']
            global_reactions[reaction]['modifiers'] = pathways[pathway]['reactions'][reaction]['modifiers']

    return global_reactions

"""Generate pathway and reaction id dictionary for accessing from global network
"""
def add_pathway_dictionary(
        pathways):

    pathway_dictionary = {}

    for x in pathways.keys():

        id = pathways[x]['reactome_id']
        reactions = list(pathways[x]['reactions'].keys())
        pathway_dictionary[id] = reactions

    return pathway_dictionary

"""
"""
def add_reactions_dictionary(
        pathways):

    reactions_dictionary = {}
    for key in pathways.keys():

        for reaction_id in pathways[key]['reactions'].keys():

            components = []

            reactants = list(pathways[key]['reactions'][reaction_id]['reactants'].keys())
            products = list(pathways[key]['reactions'][reaction_id]['products'].keys())
            modifiers = list(pathways[key]['reactions'][reaction_id]['modifiers'].keys())

            components.append(reaction_id)
            for x in [reactants, products, modifiers]:

                for y in x:
                    components.append(y)

            reactions_dictionary[reaction_id] = components

    return reactions_dictionary

"""Fetch all reactions for a given organism
"""
def __main__(
        species_id,
        output_dir,
        args_dict): # Location to output database file

    #species_id='HSA'
    #output_dir='/Users/jordan/Desktop/reactome_test'
    # Get reaction files
    reactions_dir = unpack_reactions(
        output_dir=output_dir)
    progress_feed(args_dict, "reactions")

    #reactions_dir = '/Users/jordan/Desktop/reactome_test/all_species.3.1.sbml/'
    reactome_database = {}

    # Get list of reaction files to use for populating database
    reactions_list = get_reactions(
        species_id=species_id,
        reaction_dir=reactions_dir)
    progress_feed(args_dict, "reactions")

    # Curate reactions database for organism of interest
    reactome_database['pathways'] = curate_reactions(
        reaction_dir=reactions_dir,
        reactions_list=reactions_list,
        species_id=species_id,
        args_dict=args_dict)

    # Add lists of available pathways and compartments found in the database
    reactome_database['pathway_types'] = add_pathways(
        reactome_database['pathways'])

    reactome_database['compartment_types'] = add_compartments(
        reactome_database['pathways'])

    reactome_database['global_reactions'] = add_global_reactions(
        pathways=reactome_database['pathways'])

    reactome_database['pathway_dictionary'] = add_pathway_dictionary(
        pathways=reactome_database['pathways'])

    reactome_database['reactions_dictionary'] = add_reactions_dictionary(
        pathways=reactome_database['pathways'])
    progress_feed(args_dict, "reactions")

    shutil.rmtree(reactions_dir)

    return reactome_database
