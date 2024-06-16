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
import pandas as pd
from datetime import date
import requests
import re
import os

# Import internal dependencies
from cli.metaboverse_cli.curate.load_reactions import load_reactions
from cli.metaboverse_cli.curate.load_complexes import (
    load_complexes, parse_complexes, reference_complex_species, supplement_components
)
from cli.metaboverse_cli.curate.synonyms import (
    parse_ensembl_synonyms, parse_uniprot_synonyms, parse_chebi_synonyms
)
from utils import progress_feed, write_database, write_database_json, get_metaboverse_cli_version


def parse_table(reference, key, args_dict=None):
    """Parse the provided reference table."""
    column_names = ['analyte_id', 'analyte_name', 'reaction_id', 'reaction_name']
    if 'source_id' in reference[key].columns:
        column_names.append('source_id')

    reference_parsed = reference[key][column_names].copy()
    reference_parsed['analyte'] = reference_parsed['analyte_name'].str.split(' \[').str[0]
    reference_parsed['compartment'] = reference_parsed['analyte_name'].str.split(' \[').str[1].str.split('\]').str[0]

    reference_dictionary = {}
    total = len(reference_parsed)

    for counter, (index, row) in enumerate(reference_parsed.iterrows()):
        ref_dict = {
            'analyte_id': row['analyte_id'],
            'reaction_id': row['reaction_id'],
            'reaction_name': row['reaction_name'],
            'analyte': row['analyte'],
            'compartment': row['compartment']
        }
        if 'source_id' in reference[key].columns:
            ref_dict['source_id'] = row['source_id']

        reference_dictionary[row['analyte_id']] = ref_dict

        if args_dict and counter % (total // 15) == 0:
            progress_feed(args_dict, "graph")

    return reference_dictionary


def get_reactome_version():
    """Get the most recent Reactome database version at the time of curation."""
    reactome_url = "https://reactome.org/tag/release"
    response = requests.get(reactome_url)
    matches = re.findall(r"\bv(\d+)\s+released", response.text, re.IGNORECASE)
    current_version = max(map(int, matches)) if matches else "Unknown"
    print("Current version:", current_version)
    return current_version


def add_genes(name_database, ensembl_reference):
    """Self-map all Ensembl gene IDs for network creation and mapping."""
    for gene_id in ensembl_reference.values():
        name_database[gene_id] = gene_id
    return name_database


def main(args_dict):
    """Curate the database."""
    if args_dict.get('organism_curation_file') == 'None':
        args_dict['organism_curation_file'] = os.path.join(args_dict['output'], f"{args_dict['organism_id']}.mvdb")
    args_dict['network'] = args_dict['organism_curation_file']

    print('Curating reaction network database. Please be patient, this will take several minutes...')
    print('Loading reactions...')
    args_dict, pathway_database, reaction_database, species_database, name_database, compartment_dictionary, components_database = load_reactions(
        species_id=args_dict['organism_id'],
        output_dir=args_dict['output'],
        database_source=args_dict['database_source'],
        sbml_url=args_dict['organism_curation_file'],
        args_dict=args_dict
    )

    species_database, name_database, components_database = supplement_components(
        species_database=species_database,
        name_database=name_database,
        components_database=components_database,
        compartment_dictionary=compartment_dictionary,
        species_id=args_dict['organism_id'],
        output_dir=args_dict['output']
    )

    print('Parsing ChEBI database...')
    chebi_mapper, chebi_synonyms, uniprot_metabolites = parse_chebi_synonyms(
        output_dir=args_dict['output'],
        url='ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz',
        file_name='names.tsv'
    )
    progress_feed(args_dict, "graph", 5)

    if args_dict['database_source'].lower() == 'reactome':
        print('Loading complex database...')
        complexes_reference = load_complexes(output_dir=args_dict['output'])
        progress_feed(args_dict, "graph", 2)

        print('Parsing complex database...')
        complexes_reference['complex_dictionary'] = parse_complexes(complexes_reference)
        progress_feed(args_dict, "graph", 1)

        print('Finalizing complex database...')
        complexes_reference['complex_dictionary'] = reference_complex_species(
            complexes_reference['complex_dictionary'], name_database
        )
        progress_feed(args_dict, "graph", 1)

        print('Parsing Ensembl database...')
        ensembl_reference = parse_ensembl_synonyms(
            output_dir=args_dict['output'],
            species_id=args_dict['organism_id'],
            url='https://reactome.org/download/current/Ensembl2Reactome_PE_All_Levels.txt',
            file_name='Ensembl2Reactome_PE_All_Levels.txt'
        )
        progress_feed(args_dict, "graph", 7)

        print('Adding gene IDs to name database...')
        name_database = add_genes(name_database, ensembl_reference)
        progress_feed(args_dict, "graph", 1)

        print('Parsing UniProt database...')
        uniprot_reference = parse_uniprot_synonyms(
            output_dir=args_dict['output'],
            species_id=args_dict['organism_id'],
            url='https://reactome.org/download/current/UniProt2Reactome_PE_All_Levels.txt',
            file_name='UniProt2Reactome_PE_All_Levels.txt'
        )
        progress_feed(args_dict, "graph", 3)

        database_version = f"{get_reactome_version()} (Reactome)"
    else:
        complexes_reference = {'complex_dictionary': {}}
        ensembl_reference = {}
        uniprot_reference = {}
        database_version = args_dict['database_version']

    metaboverse_db = {
        'organism_id': args_dict['organism_id'],
        'pathway_database': pathway_database,
        'reaction_database': reaction_database,
        'species_database': species_database,
        'name_database': name_database,
        'ensembl_synonyms': ensembl_reference,
        'uniprot_synonyms': uniprot_reference,
        'chebi_mapper': chebi_mapper,
        'chebi_synonyms': chebi_synonyms,
        'uniprot_metabolites': uniprot_metabolites,
        'complex_dictionary': complexes_reference['complex_dictionary'],
        'compartment_dictionary': compartment_dictionary,
        'components_database': components_database,
        'curation_date': date.today().strftime('%Y-%m-%d'),
        'metaboverse-curate_version': get_metaboverse_cli_version(),
        'database_version': database_version,
        'database_date': date.today().strftime('%Y-%m-%d')
    }

    print('Writing metaboverse database to file...')
    if args_dict['cmd'] == 'curate':
        args_dict['curation'] = f"{args_dict['organism_id']}.mvdb"
        write_database(output=args_dict['output'], file=args_dict['curation'], database=metaboverse_db)
    elif args_dict['cmd'] == 'electrum':
        args_dict['curation'] = f"{args_dict['organism_id']}.eldb"
        write_database_json(output=args_dict['output'], file=args_dict['curation'], database=metaboverse_db)
    else:
        raise Exception('Unable to output database file.')

    progress_feed(args_dict, "graph", 5)
    print('Metaboverse database curation complete.')

    return args_dict
