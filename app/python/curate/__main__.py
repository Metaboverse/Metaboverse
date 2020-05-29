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
import requests
from datetime import date
import pickle
import pandas as pd

"""Import internal dependencies
"""
from curate.load_reactions_db import __main__ as load_reactions
from curate.load_complexes_db import __main__ as load_complexes
from utils import progress_feed

def test():

    __main__(
        {'output':'/Users/jordan/Desktop/',
        'species_id':'HSA'}
    )

def parse_table(
        reference,
        key,
        args_dict=None):

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

    counter = 0
    total = len(reference_parsed.index.tolist())

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

        if int(counter % (total / 15)) == 0 and args_dict != None:
            progress_feed(args_dict, "reactions")

        counter += 1

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
                complex_dictionary[row[0]]['participants']['chebi'].append(
                    x.split(':')[1])
            if 'uniprot' in x:
                complex_dictionary[row[0]]['participants']['uniprot'].append(
                    x.split(':')[1])
            if 'ensembl' in x:
                complex_dictionary[row[0]]['participants']['ensembl'].append(
                    x.split(':')[1])
            if 'mirbase' in x:
                complex_dictionary[row[0]]['participants']['mirbase'].append(
                    x.split(':')[1])
            if 'ncbi' in x:
                complex_dictionary[row[0]]['participants']['ncbi'].append(
                    x.split(':')[1])
            else:
                pass

    return complex_dictionary

def parse_ensembl_synonyms(
        output_dir,
        species_id,
        url='https://reactome.org/download/current/Ensembl2Reactome_PE_All_Levels.txt',
        file_name='Ensembl2Reactome_PE_All_Levels.txt',
        reactome_location=3,
        name_location=2,
        id_location=0):
    """Retrieve Ensembl gene entity synonyms
    """
    # output_dir='/Users/jordan/Desktop/'

    os.system('curl -L ' + url + ' -o ' + output_dir + file_name)

    ensembl = pd.read_csv(
        output_dir + file_name,
        sep='\t',
        header=None)
    os.remove(output_dir + file_name)

    ensembl[name_location] = ensembl[name_location].str.split(' \[').str[0].tolist()

    ensembl = ensembl[ensembl[reactome_location].str.contains(species_id)]

    ensembl_name_dictionary = pd.Series(
        ensembl[name_location].values,
        index=ensembl[id_location]).to_dict()

    return ensembl_name_dictionary

def parse_uniprot_synonyms(
        output_dir,
        species_id,
        url='https://reactome.org/download/current/UniProt2Reactome_PE_All_Levels.txt',
        file_name='UniProt2Reactome_PE_All_Levels.txt',
        reactome_location=3,
        name_location=2,
        id_location=0):
    """Retrieve UniProt protein entity synonyms
    """

    os.system('curl -L ' + url + ' -o ' + output_dir + file_name)

    uniprot = pd.read_csv(
        output_dir + file_name,
        sep='\t',
        header=None)
    os.remove(output_dir + file_name)

    uniprot[name_location] = uniprot[name_location].str.split(' \[').str[0].tolist()

    uniprot = uniprot[uniprot[reactome_location].str.contains(species_id)]

    uniprot_name_dictionary = pd.Series(
        uniprot[name_location].values,
        index=uniprot[id_location]).to_dict()

    return uniprot_name_dictionary

def parse_chebi_synonyms(
        output_dir,
        url='ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz',
        file_name='names.tsv',
        name_string='NAME',
        id_string='COMPOUND_ID',
        source_string='SOURCE'):
    """Retrieve CHEBI chemical entity synonyms
    """

    os.system('curl -L ' + url + ' -o ' + output_dir + file_name + '.gz')
    os.system('gzip -d ' + output_dir + file_name + '.gz')

    chebi = pd.read_csv(
        output_dir + file_name,
        sep='\t')
    os.remove(output_dir + file_name)

    name_index = None
    id_index = None
    source_index = None

    col_names = chebi.columns.tolist()

    for x in range(len(col_names)):
        if col_names[x].upper() == name_string:
            name_index = x

        if col_names[x].upper() == id_string:
            id_index = x

        if col_names[x].upper() == source_string:
            source_index = x

    chebi_dictionary = {}
    uniprot_metabolites = {}
    if name_index != None and id_index != None:

        for index, row in chebi.iterrows():

            if 'KEGG' in row[source_index].upper() \
            or 'CHEM' in row[source_index].upper() \
            or 'CHEBI' in row[source_index].upper() \
            or 'HMDB' in row[source_index].upper() \
            or 'DRUG' in row[source_index].upper() \
            or 'IUPAC' in row[source_index].upper() \
            or 'LIPID' in row[source_index].upper() \
            or 'METACYC' in row[source_index].upper() \
            or 'SUBMITTER' in row[source_index].upper():

                chebi_dictionary[row[name_index]] = 'CHEBI:' + str(row[id_index])

            else:
                uniprot_metabolites[row[name_index]] = 'CHEBI:' + str(row[id_index])

    else:
        print('Unable to parse CHEBI file as expected...')

    return chebi_dictionary, uniprot_metabolites

def reference_complex_species(
        reference,
        name_database):
    """Correct complex dictionary keys to be searchable by species ID
    """

    new_dict = {}

    for k, v in reference.items():

        if reference[k]['complex_id'] in list(name_database.keys()):

            new_dict[name_database[reference[k]['complex_id']]] = reference[k]

    return new_dict

def get_reactome_version():
    """Get most recent Reactome database version at time of curation
    """
    reactome_url = "https://reactome.org/tag/release"
    f = requests.get(reactome_url)
    matches = re.findall("Version (.*) Released", f.text);
    current_version = max(matches)
    return current_version

def write_database(
        output,
        file,
        database):
    """Write reactions database to pickle file
    """

    # Check provided path exists
    if not os.path.isdir(output):
        os.makedirs(output)

    # Clean up path
    dir = os.path.abspath(output) + '/'

    # Write information to file
    with open(dir + file, 'wb') as file_product:
        pickle.dump(database, file_product)

def add_genes(
        name_database,
        ensembl_reference):
    """Self map all ensembl gene ids for network creation and mapping
    """

    for k, v in ensembl_reference.items():
        name_database[v] = v

    return name_database

def __main__(
        args_dict):
    """Curate reactome database
    """

    # Load reactions
    print('Curating Reactome network database. Please be patient, this will take several minutes...')
    print('Loading reactions...')
    progress_feed(args_dict, "curate", 3)
    pathway_database, reaction_database, species_database, \
    name_database, compartment_dictionary, \
    components_database = load_reactions(
        species_id=args_dict['species_id'],
        output_dir=args_dict['output'],
        args_dict=args_dict)

    print('Loading complex database...')
    complexes_reference = load_complexes(
        output_dir=args_dict['output'])
    progress_feed(args_dict, "curate", 3)

    print('Parsing complex database...')
    complexes_reference['complex_dictionary'] = parse_complexes(
        complexes_reference)
    progress_feed(args_dict, "curate", 2)

    print('Finalizing complex database...')
    complexes_reference['complex_dictionary'] = reference_complex_species(
        reference=complexes_reference['complex_dictionary'],
        name_database=name_database)
    progress_feed(args_dict, "curate", 2)

    print('Parsing Ensembl database...')
    ensembl_reference = parse_ensembl_synonyms(
        output_dir=args_dict['output'],
        species_id=args_dict['species_id'])
    progress_feed(args_dict, "curate", 3)

    print('Adding gene IDs to name database...')
    name_database = add_genes(
        name_database=name_database,
        ensembl_reference=ensembl_reference)
    progress_feed(args_dict, "curate", 2)

    print('Parsing UniProt database...')
    uniprot_reference = parse_uniprot_synonyms(
        output_dir=args_dict['output'],
        species_id=args_dict['species_id'])
    progress_feed(args_dict, "curate", 3)

    print('Parsing ChEBI database...')
    chebi_reference, uniprot_metabolites = parse_chebi_synonyms(
        output_dir=args_dict['output'])
    progress_feed(args_dict, "curate", 5)

    metaboverse_db = {
        'species_id': args_dict['species_id'],
        'pathway_database': pathway_database,
        'reaction_database': reaction_database,
        'species_database': species_database,
        'name_database': name_database,
        'ensembl_synonyms': ensembl_reference,
        'uniprot_synonyms': uniprot_reference,
        'chebi_synonyms': chebi_reference,
        'uniprot_metabolites': uniprot_metabolites,
        'complex_dictionary': complexes_reference['complex_dictionary'],
        'compartment_dictionary': compartment_dictionary,
        'components_database': components_database,
        'curation_date': date.today().strftime('%Y-%m-%d'),
        'reactome_version': get_reactome_version()
    }

    # Write database to file
    print('Writing metaboverse database to file...')
    args_dict['curation'] = args_dict['species_id'] + '_metaboverse_db.pickle'
    write_database(
        output=args_dict['output'],
        file=args_dict['curation'],
        database=metaboverse_db)
    progress_feed(args_dict, "curate", 5)

    print('Metaboverse database curation complete.')

    return args_dict
