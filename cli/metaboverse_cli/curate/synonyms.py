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

import pandas as pd
import os

from utils import download_file_ftp

def parse_ensembl_synonyms(output_dir, species_id, url, file_name):
    """Retrieve Ensembl gene entity synonyms."""
    print(f'Downloading Ensembl synonym database...\n\t{url}')
    download_file(url, output_dir, file_name)

    ensembl = pd.read_csv(f'{output_dir}{file_name}', sep='\t', header=None)
    os.remove(f'{output_dir}{file_name}')

    ensembl[2] = ensembl[2].str.split(' \[').str[0]
    ensembl = ensembl[ensembl[3].str.contains(species_id)]
    return pd.Series(ensembl[2].values, index=ensembl[0]).to_dict()


def parse_uniprot_synonyms(output_dir, species_id, url, file_name):
    """Retrieve UniProt protein entity synonyms."""
    print(f'Downloading UniProt synonym database...\n\t{url}')
    download_file(url, output_dir, file_name)

    uniprot = pd.read_csv(f'{output_dir}{file_name}', sep='\t', header=None)
    os.remove(f'{output_dir}{file_name}')

    uniprot[2] = uniprot[2].str.split(' \[').str[0]
    uniprot = uniprot[uniprot[3].str.contains(species_id)]
    return pd.Series(uniprot[2].values, index=uniprot[0]).to_dict()


def parse_chebi_synonyms(output_dir, url, file_name):
    """Retrieve CHEBI chemical entity synonyms."""
    print(f'Downloading ChEBI synonym database...\n\t{url}')
    output_path = f"{output_dir}{file_name}.gz"
    download_file_ftp(url, output_path)

    chebi = pd.read_csv(output_path, sep='\t', compression='gzip')
    os.remove(output_path)

    chebi_dictionary = {}
    chebi_synonyms = {}
    uniprot_metabolites = {}

    for _, row in chebi.iterrows():
        chebi_id = f'CHEBI:{row["COMPOUND_ID"]}'
        chebi_dictionary[row["NAME"]] = chebi_id

        if chebi_id in chebi_synonyms:
            chebi_synonyms[chebi_id].append(row["NAME"])
        else:
            chebi_synonyms[chebi_id] = [row["NAME"]]

    return chebi_dictionary, chebi_synonyms, uniprot_metabolites


def supplement_components(species_database, name_database, components_database, compartment_dictionary, species_id, output_dir, url, file_name):
    """Supplement components with ChEBI synonyms."""
    print(f'Downloading ChEBI synonym database...\n\t{url}')
    download_file(url, output_dir, file_name)

    chebi = pd.read_csv(f'{output_dir}{file_name}', sep='\t', header=None)
    os.remove(f'{output_dir}{file_name}')

    chebi = chebi[chebi[3].str.contains(species_id)]
    reversed_compartments = {v: k for k, v in compartment_dictionary.items()}

    for _, row in chebi.iterrows():
        species_id = f'species_{row[1].split("-")[-1]}'
        reactome_id = row[1]
        name = row[2].split(' [')[0] if ' [' in row[2] else row[2]
        compartment = row[2].split(' [')[1].split(']')[0] if ' [' in row[2] else None
        compartment_id = reversed_compartments.get(compartment, compartment_dictionary.get(compartment, compartment))

        if species_id not in components_database:
            components_database[species_id] = {
                'id': species_id,
                'reactome_id': reactome_id,
                'name': name,
                'is': reactome_id,
                'isEncodedBy': '',
                'hasPart': [],
                'type': 'metabolite_component',
                'compartment': compartment_id
            }

        species_database[species_id] = name
        name_database[species_id] = species_id
        name_database[name] = species_id

    return species_database, name_database, components_database


def download_file(url, output_dir, file_name):
    """Download a file from a URL and save it to the specified directory."""
    os.system(f'curl -kL {url} -o "{output_dir}{file_name}"')
