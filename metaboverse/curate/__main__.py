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
import re
import pickle

"""Import internal dependencies
"""
from metaboverse.curate.load_reactions_db import __main__ as load_reactions
from metaboverse.curate.load_chebi_db import __main__ as load_chebi
from metaboverse.curate.load_uniprot_db import __main__ as load_uniprot
from metaboverse.curate.load_ensembl_db import __main__ as load_ensembl
from metaboverse.curate.load_ncbi_db import __main__ as load_ncbi
from metaboverse.curate.load_mirbase_db import __main__ as load_mirbase
from metaboverse.curate.load_complexes_db import __main__ as load_complexes

"""Plan
"""
# This file will run loading of reactions database, chebi, ensemble, uniprot, complex, etc.
# Will then create interface dictionary for metabolites, proteins, etc relation info, name to id, etc.
# Output total network as pickle

def test():

    __main__(
        {'output':'/Users/jordan/Desktop/',
        'species':'HSA'}
    )

def parse_table(
        reference,
        key):

    if 'source_id' in reference[key].columns.tolist():
        column_names = [
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'reaction_name',
            'source_id']
    else:
        column_names = [
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'reaction_name']

    reference_parsed = reference[key][column_names].copy()
    reference_parsed['analyte'] = reference_parsed['analyte_name'].str.split(' \[').str[0]
    reference_parsed['compartment'] = reference_parsed['analyte_name'].str.split(' \[').str[1].str.split('\]').str[0]

    reference_dictionary = {}

    for index, row in reference_parsed.iterrows():

        reference_dictionary[row[0]] = {}
        reference_dictionary[row[0]]['analyte_id'] = row[0]
        reference_dictionary[row[0]]['reaction_id'] = row[2]
        reference_dictionary[row[0]]['reaction_name'] = row[3]

        if 'source_id' in reference[key].columns.tolist():
            reference_dictionary[row[0]]['source_id']  = row[4]
            reference_dictionary[row[0]]['analyte'] = row[5]
            reference_dictionary[row[0]]['compartment'] = row[6]

        else:
            reference_dictionary[row[0]]['analyte'] = row[4]
            reference_dictionary[row[0]]['compartment'] = row[5]

    return reference_dictionary

def parse_complexes(
        reference):

    pathway_dictionary = {}
    for index, row in reference['complex_pathway'].iterrows():

        pathway_dictionary[row[0]] = {}
        pathway_dictionary[row[0]]['complex'] = row[0]
        pathway_dictionary[row[0]]['pathway'] = row[1]
        pathway_dictionary[row[0]]['top_level_pathway'] = row[2]

    column_names = [
        'identifier',
        'name',
        'participants',
        'participatingComplex']
    complexes_information = reference['complex_participants'][column_names].copy()
    complexes_information['complex'] = complexes_information['name'].str.split(' \[').str[0]
    complexes_information['compartment'] = complexes_information['name'].str.split(' \[').str[1].str.split('\]').str[0]

    complex_dictionary = {}

    for index, row in complexes_information.iterrows():

        complex_dictionary[row[0]] = {}
        complex_dictionary[row[0]]['complex_id'] = row[0]
        complex_dictionary[row[0]]['complex_name'] = row[4]
        complex_dictionary[row[0]]['compartment'] = row[5]

        if row[3] == '-':
            complex_dictionary[row[0]]['participating_complex'] = None
        else:
            complex_dictionary[row[0]]['participating_complex'] = row[3]

        if row[0] in pathway_dictionary.keys():
            complex_dictionary[row[0]]['pathway'] = pathway_dictionary[row[0]]['pathway']
            complex_dictionary[row[0]]['top_level_pathway'] = pathway_dictionary[row[0]]['top_level_pathway']

        complex_dictionary[row[0]]['participants'] = {}
        complex_dictionary[row[0]]['participants']['chebi'] = []
        complex_dictionary[row[0]]['participants']['uniprot'] = []
        complex_dictionary[row[0]]['participants']['ensembl'] = []
        complex_dictionary[row[0]]['participants']['mirbase'] = []
        complex_dictionary[row[0]]['participants']['ncbi'] = []

        participants = row[2].split('|')

        for x in participants:
            if 'chebi' in x:
                complex_dictionary[row[0]]['participants']['chebi'].append(x.split(':')[1])
            if 'uniprot' in x:
                complex_dictionary[row[0]]['participants']['uniprot'].append(x.split(':')[1])
            if 'ensembl' in x:
                complex_dictionary[row[0]]['participants']['ensembl'].append(x.split(':')[1])
            if 'mirbase' in x:
                complex_dictionary[row[0]]['participants']['mirbase'].append(x.split(':')[1])
            if 'ncbi' in x:
                complex_dictionary[row[0]]['participants']['ncbi'].append(x.split(':')[1])
            else:
                pass

    return complex_dictionary

"""
"""
def make_master(
        database):

    master_reference = {}

    reference_list = [
        'chebi_reference',
        'uniprot_reference',
        'ensembl_reference',
        'ncbi_reference',
        'mirbase_reference']

    for x in reference_list:

        for key in database[x].keys():

            master_reference[database[x][key]['analyte_id']] = database[x][key]['analyte']

    for key in database['complexes_reference']['complex_dictionary'].keys():

        master_reference[database['complexes_reference']['complex_dictionary'][key]['complex_id']] = database['complexes_reference']['complex_dictionary'][key]['complex_name']

    # Add pathway id to name mapping
    for key in database['global_reactions'].keys():

        master_reference[database['global_reactions'][key]['pathway_id']] = database['global_reactions'][key]['pathway_name']

    return master_reference

"""
"""
def prepare_mappers(
        databases):

    mappers = {}

    mappers['ensembl'] = {}
    for key in databases['ensembl_reference'].keys():

        mappers['ensembl'][databases['ensembl_reference'][key]['source_id']] = databases['ensembl_reference'][key]['analyte_id']

    mappers['mirbase'] = {}
    for key in databases['mirbase_reference'].keys():

        mappers['mirbase'][databases['mirbase_reference'][key]['source_id']] = databases['mirbase_reference'][key]['analyte_id']

    mappers['uniprot'] = {}
    for key in databases['uniprot_reference'].keys():

        mappers['uniprot'][databases['uniprot_reference'][key]['source_id']] = databases['uniprot_reference'][key]['analyte_id']

    mappers['chebi'] = {}
    for key in databases['chebi_reference'].keys():

        mappers['chebi'][databases['chebi_reference'][key]['source_id']] = databases['chebi_reference'][key]['analyte_id']

    return mappers

"""Map complex names to lists of component genes
"""
def map_complexes_genes(
        complex_participants,
        databases):

    complex_mapper = {}

    mappers = prepare_mappers(
        databases=databases)

    for index, row in complex_participants.iterrows():

        complex_id = row[0]
        participants = row[2].split('|')
        partner_complexes = row[3].split('|')

        participants_ids = []

        for x in participants:

            if x.split(':')[0] == 'ensembl':
                try:
                    id = mappers['ensembl'][x.split(':')[1]]
                    participants_ids.append(id)

                except:
                    print('Could not find', x)

            elif x.split(':')[0] == 'mirbase':
                try:
                    id = mappers['mirbase'][x.split(':')[1]]
                    participants_ids.append(id)

                except:
                    print('Could not find', x)

            elif x.split(':')[0] == 'ncbi':
                print('Skipping NCBI records...')

            elif x.split(':')[0] == 'uniprot':
                try:
                    id = mappers['uniprot'][x.split(':')[1]]
                    participants_ids.append(id)

                except:
                    print('Could not find', x)

            elif x.split(':')[0] == 'chebi':
                try:
                    id = mappers['chebi'][x.split(':')[1]]
                    participants_ids.append(id)

                except:
                    print('Could not find', x)

            else:
                print('Warning: Unable to find', x, 'in database')

        complex_mapper[complex_id] = {
            'id': complex_id,
            'participants': participants_ids,
            'partner_complexes': partner_complexes}

    return complex_mapper

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

"""Curate reactome database
"""
def __main__(
        args_dict,
        log_file):

    #args_dict = {
    #    'output':'/Users/jordan/Desktop/',
    #    'species':'HSA'}

    # Load reactions
    print('Curating Reactome network database. Please be patient, this will take several minutes...')
    print('Loading reactions...')
    reactions_database = load_reactions(
        species_id=args_dict['species'],
        output_dir=args_dict['output'])

    # Add interfacing information to reactions database
    print('\nLoading ChEBI database...')
    chebi_reference = load_chebi(
        output_dir=args_dict['output'])
    print('Parsing ChEBI database...')
    reactions_database['chebi_reference'] = parse_table(
        reference=chebi_reference,
        key='chebi_pe_all_levels')

    print('Loading UniProt database...')
    uniprot_reference = load_uniprot(
        output_dir=args_dict['output'])
    print('Parsing UniProt database...')
    reactions_database['uniprot_reference'] = parse_table(
        reference=uniprot_reference,
        key='uniprot_pe_all_levels')

    print('Loading Ensembl database...')
    ensembl_reference = load_ensembl(
        output_dir=args_dict['output'])
    print('Parsing Ensembl database...')
    reactions_database['ensembl_reference'] = parse_table(
        reference=ensembl_reference,
        key='ensembl_pe_all_levels')

    print('Loading NCBI database...')
    ncbi_reference = load_ncbi(
        output_dir=args_dict['output'])
    print('Parsing NCBI database...')
    reactions_database['ncbi_reference'] = parse_table(
        reference=ncbi_reference,
        key='ncbi_pe_all_levels')

    print('Loading miRBase database...')
    mirbase_reference = load_mirbase(
        output_dir=args_dict['output'])
    print('Parsing miRBase database...')
    reactions_database['mirbase_reference'] = parse_table(
        reference=mirbase_reference,
        key='mirbase_pe_all_levels')

    print('Loading complex database...')
    reactions_database['complexes_reference'] = load_complexes(
        output_dir=args_dict['output'])
    print('Parsing complex database...')
    reactions_database['complexes_reference']['complex_dictionary'] = parse_complexes(
            reactions_database['complexes_reference'])

    reactions_database['complexes_reference']['complex_mapper'] = map_complexes_genes(
        complex_participants=reactions_database['complexes_reference']['complex_participants'],
        databases=reactions_database)

    reactions_database['master_reference'] = make_master(
        reactions_database)

    # Write database to file
    print('Writing metaboverse database to file...')
    write_database(
        output=args_dict['output'],
        file=args_dict['species'] + '_metaboverse_db.pickle',
        database=reactions_database)
    print('Metaboverse database curation complete.')
