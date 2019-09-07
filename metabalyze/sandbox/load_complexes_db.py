"""Import dependencies
"""
import pandas as pd

"""Import internal dependencies
"""
#from load_utils import get_table
"""Including until modules in package======================================
"""
import os
import sys
import pandas as pd

def get_table(
        output_dir,
        url,
        column_names=None,
        header=None,
        organism='Homo sapiens',
        organism_key='organism'):

    # ncbi_reactom_reactions
    file = unpack_table(
            url=url,
            output_dir=output_dir)
    data = pd.read_csv(
        file,
        sep='\t',
        header=header,
        low_memory=False)

    if column_names != None:
        data.columns = column_names
        data_organism = data.loc[data[organism_key] == organism]
    else:
        data_organism = data

    return data_organism

def unpack_table(
        url,
        output_dir='./'):

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    file = output_dir + url.split('/')[-1]

    os.system('curl -L ' + url + ' -o ' + file)

    return file
"""=======================================================================
"""

def __main__(
        output_dir):

    complex_participants = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/ComplexParticipantsPubMedIdentifiers_human.txt',
        header=0)

    complex_pathway = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/Complex_2_Pathway_human.txt',
        header=0)

    return {
        'complex_participants': complex_participants,
        'complex_pathway': complex_pathway}



output_dir = '/Users/jordan/Desktop/reactome_test/'
complex = __main__(
    output_dir)


complex['complex_participants'].shape
complex['complex_pathway'].shape


complex['complex_participants'].head()

complex['complex_participants'].loc[complex['complex_participants']['name'].str.contains('MPC')]



complex['complex_pathway'].head()
