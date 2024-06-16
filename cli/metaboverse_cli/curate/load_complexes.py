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

import os
import pandas as pd
from curate.utils import get_table

def load_complexes(output_dir):
    """
    Load complex participants and pathways from Reactome.

    Parameters:
        output_dir (str): Directory to save the downloaded files.

    Returns:
        dict: Dictionary containing complex participants and pathways data.
    """
    complex_participants_url = 'https://reactome.org/download/current/ComplexParticipantsPubMedIdentifiers_human.txt'
    complex_pathway_url = 'https://reactome.org/download/current/Complex_2_Pathway_human.txt'

    complex_participants = get_table(
        output_dir=output_dir,
        url=complex_participants_url,
        column_names=0
    )
    os.remove(os.path.join(output_dir, os.path.basename(complex_participants_url)))

    complex_pathway = get_table(
        output_dir=output_dir,
        url=complex_pathway_url,
        column_names=0
    )
    os.remove(os.path.join(output_dir, os.path.basename(complex_pathway_url)))

    return {
        'complex_participants': complex_participants,
        'complex_pathway': complex_pathway
    }


def parse_complexes(reference):
    """
    Parse complex data from the provided reference.

    Parameters:
        reference (dict): Dictionary containing complex participants and pathways data.

    Returns:
        dict: Parsed complex data.
    """
    def create_pathway_dictionary(complex_pathway_df):
        return {
            row['complex']: {
                'complex': row['complex'],
                'pathway': row['pathway'],
                'top_level_pathway': row['top_level_pathway']
            }
            for _, row in complex_pathway_df.iterrows()
        }

    def prepare_complexes_information(complex_participants_df):
        complexes_information = complex_participants_df[['identifier', 'name', 'participants', 'participatingComplex']].copy()
        complexes_information['complex'] = complexes_information['name'].str.split(' \[').str[0]
        complexes_information['compartment'] = complexes_information['name'].str.extract(r'\[(.*?)\]')
        return complexes_information

    def parse_participants(participants_str):
        participants = {'chebi': [], 'uniprot': [], 'ensembl': [], 'mirbase': [], 'ncbi': []}
        for part in participants_str.split('|'):
            for key in participants:
                if key in part:
                    participants[key].append(part.split(':')[1])
        return participants

    def create_complex_dictionary(complexes_information, pathway_dictionary):
        complex_dictionary = {}
        for _, row in complexes_information.iterrows():
            complex_id = row['identifier']
            complex_dictionary[complex_id] = {
                'complex_id': complex_id,
                'complex_name': row['complex'],
                'compartment': row['compartment'] if pd.notna(row['compartment']) else None,
                'participating_complex': None if row['participatingComplex'] == '-' else row['participatingComplex'],
                'pathway': pathway_dictionary.get(complex_id, {}).get('pathway', None),
                'top_level_pathway': pathway_dictionary.get(complex_id, {}).get('top_level_pathway', None),
                'participants': parse_participants(row['participants'])
            }
        return complex_dictionary

    pathway_dictionary = create_pathway_dictionary(reference['complex_pathway'])
    complexes_information = prepare_complexes_information(reference['complex_participants'])
    return create_complex_dictionary(complexes_information, pathway_dictionary)


def reference_complex_species(reference, name_database):
    """
    Correct complex dictionary keys to be searchable by species ID.

    Parameters:
        reference (dict): Dictionary containing parsed complex data.
        name_database (dict): Dictionary containing name data.

    Returns:
        dict: Updated complex data with species ID as keys.
    """
    return {
        name_database.get(v['complex_id'], v['complex_id']): v
        for v in reference.values()
    }
