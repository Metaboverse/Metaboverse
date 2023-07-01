"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) 2022 Metaboverse

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
import numpy as np
import ast
import sys
import re


def eval_table(table):
    
    return table.applymap(lambda x: ast.literal_eval(str(x)))
    
    
# Source: https://www.geeksforgeeks.org/python-program-to-flatten-a-nested-list-using-recursion/
def flattenList(nestedList):
 
    # check if list is empty
    if not(bool(nestedList)):
        return nestedList
 
     # to check instance of list is empty or not
    if isinstance(nestedList[0], list):
 
        # call function with sublist as argument
        return flattenList(*nestedList[:1]) + flattenList(nestedList[1:])
 
    # call function with sublist as argument
    return nestedList[:1] + flattenList(nestedList[1:])


def read_data(
        url,
        delimiter='\t',
        duplicates=False):
    """Expected to contain an even number of columns for the data type
    """

    data = pd.read_csv(
        url,
        sep=delimiter,
        index_col=0)

    data = data.dropna(axis=1, how="all")
    
    # handle duplicate indices 
    if len(data.loc[data.index.duplicated()].index.tolist()) != 0:
        print(
            'Warning: Input data table contained duplicate row names. ' 
            + 'Duplicates will be removed.'
            + '\n\tRemoving '
                + str(len(data.loc[data.index.duplicated(keep=duplicates)].index.tolist()))
                + ' rows.'
            + '\n\tFile: ' 
                + str(url)
            + '\nEntities removed:' 
            )
        for x in data.loc[data.index.duplicated(keep=duplicates)].index.tolist():
            print("\t- " + str(x))
        data = data.loc[~data.index.duplicated(keep=duplicates)] 

    if len(data.columns.tolist()) % 2 != 0:
        raise Exception('Improperly formatted datatable provided: ', url)

    return data


def check_data(data, data_type="unknown"):

    print("{} data conversion errors:".format(data_type))
    should_exit = False 
    _count = 0

    # Check cell values
    for r in data.index:
        for c in data.columns:
            try:
                float(data.at[r, c])
            except ValueError:
                print("    Formatting error in data cell at coordinates: {0}, {1}  (row {2})".format(r, c, data.index.get_loc(r)))
                _count += 1
                should_exit = True
    
    # Check index values
    for r in data.index:
        try:
            str(r.lstrip().rstrip())
        except AttributeError:
            print("    Formatting error in index: {0}  (row {1})".format(r, data.index.get_loc(r)))
            _count += 1
            should_exit = True

    if _count == 0:
        print("    None")

    return should_exit


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
    """

    data_output = data.copy()
    data_output.index = data_output.index.str.upper()

    reference_ids = list(reference.values())

    data_output.index = data_output.index.to_series().replace(reference)
    data_unmapped = data_output.copy()

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
    data_c = data_c.dropna(axis=0)
    data_c.index = [d.lstrip().rstrip() for d in data_c.index.tolist()]

    _values = data_c.T[::2].T
    _stats = data_c.T[1::2].T

    _values.columns = [x for x in range(len(_values.columns))]
    _stats.columns = [x for x in range(len(_stats.columns))]

    return _values, _stats


def broadcast_transcriptomics(
        transcriptomics,
        transcriptomics_stats,
        gene_dictionary,
        protein_dictionary):
    """If transcript levels are provided and no protein levels are given,
    the user can opt to broadcast the gene expression levels to the protein
    nodes. If an enzyme node is constituted of several entities, the user can
    choose to use the average or use the most extreme value (the fc furthest
    from 0)
    """

    proteomics = transcriptomics.copy()
    proteomics_stats = transcriptomics_stats.copy()

    uniprot_dict = {}
    for x in protein_dictionary.keys():
        uniprot_dict[protein_dictionary[x]] = x

    proteomics.index = proteomics.index.to_series().replace(
        gene_dictionary)
    proteomics_stats.index = proteomics_stats.index.to_series().replace(
        gene_dictionary)

    proteomics.index = proteomics.index.to_series().replace(
        uniprot_dict)
    proteomics_stats.index = proteomics_stats.index.to_series().replace(
        uniprot_dict)

    return proteomics, proteomics_stats


def copy_columns(
        data,
        stats,
        _max):
    """Copy data columns for every time point provided by the user
    """

    data_c = data.copy()
    stats_c = stats.copy()

    counter = 1

    while counter < _max:
        data_c[counter] = data_c[0]
        stats_c[counter] = stats_c[0]
        counter += 1

    return data_c, stats_c


def convert_to_float(
        lists):
    """Recursively force elements of nested list to float
    Source: https://stackoverflow.com/a/46549533

    Args:
        lists (<list>): Nested list

    Returns:
        [<list>]: Nested list with elements converted to floats
    """
    
    return [float(el) if not isinstance(el, list) else convert_to_float(el) for el in lists]
                
                
def catenate_data(
        array):
    """Combine data
    - average duplicates
    - remove NAs
    Return: combined dataframe where indices are species IDs
    """

    tables = [eval_table(x) for x in array]

    combined = pd.concat(tables)
    combined = combined.dropna(axis=0)
    combined = combined.sort_index()

    # Check that types are the same (p-values or confidence intervals)
    if type(combined.iloc[0,0]) == list:
        if len(set(flattenList(combined.applymap(type).values.tolist()))) != 1:
            raise Exception("Input data types do not match. Please check that all fold change and statistical value types match between datasets.")

    ### Need to fix
    removers = []  # Remove non-numbers
    for idx, row in combined.iterrows():
        for x in row:
            if type(x) == list:
                try:
                    x = convert_to_float(x)
                except:
                    removers.append(idx)
            else: 
                try:
                    float(x)
                except:
                    removers.append(idx)

    combined = combined[~combined.index.isin(removers)]
    
    return combined


def __main__(
        network,
        transcriptomics_url,
        proteomics_url,
        metabolomics_url,
        database_source='reactome'):
    """Get user data and preprocess
    """

    should_transcriptomics_exit = False
    should_proteomics_exit = False
    should_metabolomics_exit = False

    # Process transcriptomics
    if transcriptomics_url.lower() != 'none':

        transcriptomics = read_data(
            url=transcriptomics_url)
        should_transcriptomics_exit = check_data(transcriptomics, data_type="Transcriptomics")
        if not should_transcriptomics_exit:
            e_sym = {}
            if database_source.lower() == 'reactome':
                for k, v in network['ensembl_synonyms'].items():
                    if 'phospho-' in v and '-phospho-' not in v:
                        v = v.replace('phospho-', '')
                    if '(' in v and ')' in v:
                        v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
                    if '  ' in v:
                        v = v.replace('  ', ' ')
                    v = v.strip()
                    e_sym[v.upper()] = k
                    e_sym[k.upper()] = k
                transcriptomics, transcriptomics_unmapped = format_data(
                    data=transcriptomics,
                    reference=e_sym)
                output_unmapped(
                    data=transcriptomics_unmapped,
                    url=transcriptomics_url)
            transcriptomics, transcriptomics_stats = extract_data(
                data=transcriptomics)
            transcriptomics_length = len(transcriptomics.columns.tolist())
    else:
        transcriptomics_length = 0

    # Process proteomics
    if proteomics_url.lower() != 'none':

        proteomics = read_data(
            url=proteomics_url)
        should_proteomics_exit = check_data(proteomics, data_type="Proteomics")
        if not should_proteomics_exit:
            u_sym = {}
            if database_source.lower() == 'reactome':
                for k, v in network['uniprot_synonyms'].items():
                    if 'phospho-' in v and '-phospho-' not in v:
                        v = v.replace('phospho-', '')
                    if '(' in v and ')' in v:
                        v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
                    if '  ' in v:
                        v = v.replace('  ', ' ')
                    v = v.strip()
                    u_sym[v.upper()] = k
                    u_sym[k.upper()] = k
                proteomics, proteomics_unmapped = format_data(
                    data=proteomics,
                    reference=u_sym)
                output_unmapped(
                    data=proteomics_unmapped,
                    url=proteomics_url)
            proteomics, proteomics_stats = extract_data(
                data=proteomics)
            proteomics_length = len(proteomics.columns.tolist())
    else:
        proteomics_length = 0

    # Process metabolomics
    if metabolomics_url.lower() != 'none':

        metabolomics = read_data(
            url=metabolomics_url)
        should_metabolomics_exit = check_data(metabolomics, data_type="Metabolomics")
        if not should_metabolomics_exit:
            metabolomics, metabolomics_stats = extract_data(
                data=metabolomics)
            metabolomics_length = len(metabolomics.columns.tolist())
    else:
        metabolomics_length = 0

    # Check if should exit if encountered invalid value in input data
    if should_transcriptomics_exit \
    or should_proteomics_exit \
    or should_metabolomics_exit:
        print("""
        ----------------------------------------------------------------
           Formatting/conversion error(s) in data. Please check logs, 
                   remedy problematic data points, and rerun.            
        ----------------------------------------------------------------
        """)
        sys.exit(1)

    # Check for broadcasting
    if proteomics_url.lower() == 'none' \
            and transcriptomics_url.lower() != 'none':

        proteomics, proteomics_stats = broadcast_transcriptomics(
            transcriptomics=transcriptomics,
            transcriptomics_stats=transcriptomics_stats,
            gene_dictionary=network['ensembl_synonyms'],
            protein_dictionary=network['uniprot_synonyms'])

    # Allow for unequal filling
    lengths = []
    for x in [
        transcriptomics_length,
        proteomics_length,
        metabolomics_length,
    ]:
        if x != 0:
            lengths.append(x)
        else:
            lengths.append(1)

    if sum(lengths) != 3 and len(list(set(lengths))) == 2:

        _max = max(lengths)

        if transcriptomics_url.lower() != 'none' \
        and len(transcriptomics.columns.tolist()) != _max \
        and len(transcriptomics.columns.tolist()) == 1:
            transcriptomics, transcriptomics_stats = copy_columns(
                data=transcriptomics,
                stats=transcriptomics_stats,
                _max=_max)

        if proteomics_url.lower() != 'none' \
        and len(proteomics.columns.tolist()) != _max \
        and len(proteomics.columns.tolist()) == 1:
            proteomics, proteomics_stats = copy_columns(
                data=proteomics,
                stats=proteomics_stats,
                _max=_max)

        if metabolomics_url.lower() != 'none' \
        and len(metabolomics.columns.tolist()) != _max \
        and len(metabolomics.columns.tolist()) == 1:
            metabolomics, metabolomics_stats = copy_columns(
                data=metabolomics,
                stats=metabolomics_stats,
                _max=_max)

    elif len(list(set(lengths))) == 1:
        pass
    else:
        print("When providing multi-omic timecourse data with unequals times, other omics types must match the maximum number of time points or must only provide one time point.")

    # Initialize array of filled data
    data_array = []
    stats_array = []

    if transcriptomics_url.lower() != 'none':
        data_array.append(transcriptomics)
        stats_array.append(transcriptomics_stats)
    if proteomics_url.lower() != 'none':
        data_array.append(proteomics)
        stats_array.append(proteomics_stats)
    if metabolomics_url.lower() != 'none':
        data_array.append(metabolomics)
        stats_array.append(metabolomics_stats)

    # Combine data
    data = catenate_data(
        array=data_array)

    stats = catenate_data(
        array=stats_array)

    unmapped = {}
    try:
        unmapped['transcriptomics_unmapped'] = \
            transcriptomics_unmapped.index.tolist()
    except:
        unmapped['transcriptomics_unmapped'] = []
    try:
        unmapped['proteomics_unmapped'] = \
            proteomics_unmapped.index.tolist()
    except:
        unmapped['proteomics_unmapped'] = []

    return data, stats, unmapped
