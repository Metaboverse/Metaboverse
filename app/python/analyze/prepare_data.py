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
import pandas as pd

def read_data(
        url,
        delimiter='\t'):
    """Expected to contain an even number of columns for the data type
    """

    data = pd.read_csv(
        url,
        sep=delimiter,
        index_col=0)

    if len(data.columns.tolist()) % 2 != 0:
        raise Exception('Improperly formatted datatable provided: ', url)

    return data

def format_data(
        data,
        reference):
    """Format data for processing
    0) Transcriptomics
    - Accepts entities mapping to ensembl genes and their synonyms
    - Expects fold change values and their p-values side-by-side
    - Columns should be ordered in intended order (for multi-conditions)
    - gene_dictionary is an Ensembl GTF file
    - Columns should be condition_fc, condition_p, etc.
    1) Proteomics
    - Accepts entities mapping to uniprot proteins and their synonyms
    - Expects fold change values and their p-values side-by-side
    - Columns should be ordered in intended order (for multi-conditions)
    - gene_dictionary is an Ensembl GTF file
    - Columns should be condition_fc, condition_p, etc.
    3) Metabolomics
    - Accepts entities mapping to ensembl genes and their synonyms
    - Expects fold change values and their p-values side-by-side
    - Columns should be ordered in intended order (for multi-conditions)
    - gene_dictionary is an Ensembl GTF file
    - Columns should be condition_fc, condition_p, etc.
    """

    data_output = data.copy()
    data_unmapped = data.copy()

    reference_ids = list(reference.values())

    data_output.index = data_output.index.to_series().map(reference)
    data_output = data_output[data_output.index.isin(reference_ids)]

    data_unmapped.index = data_unmapped.index.to_series().map(reference)
    data_unmapped = data_unmapped[~data_unmapped.index.isin(reference_ids)]

    return data_output, data_unmapped

def output_unmapped(
        data,
        url,
        delimiter='\t'):
    """Output unmapped entities for user information
    """

    if len(data.index.tolist()):
        data.to_csv(
            url[:-4] + '_unmapped.txt',
            sep=delimiter)

def extract_data(
        data):
    """Seperate out fold change and p-values
    - Formatting follows fold change, p-values, etc...
    """
    data_c = data.copy()

    _values = data_c.T[::2].T
    _stats = data_c.T[1::2].T

    return _values, _stats

def broadcast_transcriptomics(
        transcriptomics,
        gene_dictionary,
        protein_dictionary):
    """If transcript levels are provided and no protein levels are given,
    the user can opt to broadcast the gene expression levels to the protein
    nodes. If an enzyme node is constituted of several entities, the user can
    choose to use the average or use the most extreme value (the fc furthest
    from 0)
    """



def catenate_data(
        array):
    """Combine data
    - average duplicates
    - remove NAs
    Return: combined dataframe where indices are species IDs
    """

    combined = pd.concat(array)
    combined = combined.dropna(axis=0)

    combined = combined.sort_index()
    combined = combined.groupby(combined.index).mean()

    return combined

def __main__(
        network,
        transcriptomics_url,
        proteomics_url,
        metabolomics_url):
    """Get user data and preprocess
    """
    #############################
    def read_network(
            network_url):
        """Read in network from previous curation module
        - was provided as a URL to the file and saved to args_dict['network'] in
        "curate" sub-module
        """
        import pickle
        with open(network_url, 'rb') as network_file:
            network = pickle.load(network_file)

        return network

    network = read_network(
        network_url='/Users/jordan/Desktop/HSA_metaboverse_db.pickle')

    transcriptomics_url='/Users/jordan/Desktop/metaboverse/app/python/analyze/test/transcriptomics.txt'
    proteomics_url='/Users/jordan/Desktop/metaboverse/app/python/analyze/test/proteomics.txt'
    metabolomics_url='/Users/jordan/Desktop/metaboverse/app/python/analyze/test/metabolomics.txt'

    #############################

    # Initialize array of filled data
    data_array = []
    stats_array = []

    # Process transcriptomics
    if transcriptomics_url.lower() != 'none':

        transcriptomics = read_data(
            url=transcriptomics_url)
        transcriptomics, transcriptomics_unmapped = format_data(
            data=transcriptomics,
            reference=network['ensembl_synonyms'])
        output_unmapped(
            data=transcriptomics_unmapped,
            url=transcriptomics_url)
        transcriptomics, transcriptomics_stats = extract_data(
            data=transcriptomics)
        data_array.append(transcriptomics)
        stats_array.append(transcriptomics_stats)



transcriptomics_unmapped

    # Process proteomics
    if proteomics_url.lower() != 'none':

        proteomics, proteomics_stats = read_data(
            url=proteomics_url)
        proteomics, proteomics_unmapped = format_proteomics(
            data=proteomics,
            reference=network['uniprot_synonyms'])
        output_unmapped(
            data=proteomics_unmapped,
            url=proteomics_url)
        proteomics, proteomics_stats = extract_data(
            data=proteomics)
        data_array.append(proteomics)
        stats_array.append(proteomics_stats)

    # Process metabolomics
    if metabolomics_url.lower() != 'none':

        metabolomics = read_data(
            url=proteomics_url)
        metabolomics, metabolomics_unmapped = format_metabolomics(
            data=metabolomics,
            reference=network['chebi_synonyms'])
        output_unmapped(
            data=metabolomics_unmapped,
            url=metabolomics_url)
        metabolomics, metabolomics_stats = extract_data(
            data=metabolomics)
        data_array.append(metabolomics)
        stats_array.append(metabolomics_stats)

    # Check for broadcasting
    if proteomics_url.lower() == 'none' \
    and transcriptomics_url.lower() != 'none':

        #


    # Combine data
    data = catenate_data(
        array=data_array)

    stats = catenate_data(
        array=stats_array)

    return data, stats
