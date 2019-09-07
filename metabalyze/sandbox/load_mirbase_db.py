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

    all_levels = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/miRBase2Reactome_All_Levels.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'url',
            'reaction_id',
            'go_evidence', #TAS = traceable author statement, IEA = electrong annotation not manually reviewed
            'organism'])

    pe_all_levels = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/miRBase2Reactome_PE_All_Levels.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    pe_pathways = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/miRBase2Reactome_PE_Pathway.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    pe_reactions = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/miRBase2Reactome_PE_Reactions.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    reactome = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/miRBase2Reactome.txt',
        column_names=[
            'source_id',
            'process_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    reactome_reactions = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/miRBase2ReactomeReactions.txt',
        column_names=[
            'source_id',
            'process_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    return {
        'mirbase_all_levels': all_levels,
        'mirbase_pe_all_levels': pe_all_levels,
        'mirbase_pe_pathways': pe_pathways,
        'mirbase_pe_reactions': pe_reactions,
        'mirbase_reactome': reactome,
        'mirbase_reactome_reactions': reactome_reactions}



output_dir = '/Users/jordan/Desktop/reactome_test/'
mirbase = __main__(
    output_dir)


mirbase['mirbase_all_levels'].shape
mirbase['mirbase_pe_all_levels'].shape
mirbase['mirbase_pe_pathways'].shape
mirbase['mirbase_pe_reactions'].shape

mirbase['mirbase_pe_all_levels'].head()
