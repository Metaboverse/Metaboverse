"""Import dependencies
"""
import os
import sys
import pandas as pd

def get_table(
        output_dir,
        url,
        column_names,
        organism='Homo sapiens',
        organism_key='organism'):

    # chebi_reactom_reactions
    file = unpack_table(
            url=url,
            output_dir=output_dir)
    data = pd.read_csv(
        file,
        sep='\t',
        header=None)
    data.columns = column_names
    data_organism = data.loc[data[organism_key] == organism]

    return data_organism

def unpack_table(
        url,
        output_dir='./'):

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    file = output_dir + url.split('/')[-1]

    os.system('curl -L ' + url + ' -o ' + file)

    return file
