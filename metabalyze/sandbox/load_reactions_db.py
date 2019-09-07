"""Import dependencies
"""
import os
import sys
import time
import pickle
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from metabalyze.utils import progress_bar

"""Set globals
"""
analyte_prefix='R-ALL-'

"""Run tests
"""
def test():
    reactions = __main__(
            species_id='HSA',
            reaction_dir='/Users/jordan/Desktop/reactome_test/all_species.3.1.sbml',
            output_dir='/Users/jordan/Desktop/reactome_test')

    reactions['R-HSA-2562578']['reactome_id']
    reactions['R-HSA-2562578']['pathway_name']

    reactions['R-HSA-2562578']['reactions'].keys()
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['name']

    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['id']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['name']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['reversible']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['fast']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['compartment']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['reactants'].keys()
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['products'].keys()
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['products']['R-ALL-2562577']['species_id']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['products']['R-ALL-2562577']['type']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['products']['R-ALL-2562577']['constant']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['products']['R-ALL-2562577']['stoichiometry']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['modifiers'].keys()

    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['reactants']['R-ALL-2562542']
    reactions['R-HSA-2562578']['reactions']['R-HSA-2562564']['modifiers']['R-ALL-2562542']


    id = 'R-HSA-162865'
    for k in reactions.keys():
        if id in k:
            print('1')
            print(reactions[k]['reactome_id'])
            print(reactions[k]['pathway_name'])
        for kk in reactions[k]['reactions'].keys():
            if id in reactions[k]['reactions'][kk]['reactants'].keys() \
            or id in reactions[k]['reactions'][kk]['products'].keys() \
            or id in reactions[k]['reactions'][kk]['modifiers'].keys():
                print('2')
                print(reactions[k]['reactome_id'])
                print(reactions[k]['pathway_name'])





    unpack_reactions(
        output_dir='/Users/jordan/Desktop/')

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

"""Curate database for all reactions and processes for an organism
"""
def curate_reactions(
        reaction_dir,
        reactions_list,
        species_id,
        species_key='species'):

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

        progress_bar(
            counter,
            total,
            status='Processing reactions')
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

"""Write reactions database to pickle file
"""
def write_database(
        output,
        file,
        database):

    # Check provided path exists
    if not os.path.isdir(output):
        os.makedirs(output)

    # Clean up path
    dir = os.path.abspath(output) + '/'

    # Write information to file
    with open(dir + file, 'wb') as file_product:
        pickle.dump(database, file_product)

"""Load tarballed sbml reactome reaction files from reactome site
"""
def unpack_reactions(
        output_dir,
        url='https://reactome.org/download/current/homo_sapiens.3.1.sbml.tgz'):

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    file = output_dir + url.split('/')[-1]

    os.system('curl -L ' + url + ' -o ' + file)

    reactions_dir = file[:-4] + '/'
    os.makedirs(reactions_dir)
    os.system('tar -zxvf ' + file + ' -C ' + reactions_dir)

    return reactions_dir

"""Fetch all reactions for a given organism
"""
def __main__(
        species_id,
        reaction_dir,
        output_dir):

    reactions_list = get_reactions(
        species_id=species_id,
        reaction_dir=reaction_dir)

    reactions_database = curate_reactions(
        reaction_dir=reaction_dir,
        reactions_list=reactions_list,
        species_id=species_id)

    write_database(
        output=output_dir,
        file=species_id + '_reactome_reactions.pickle',
        database=reactions_database)

    return reactions_database
