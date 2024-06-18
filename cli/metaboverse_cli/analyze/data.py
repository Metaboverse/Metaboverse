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
from scipy.stats import gmean
import statistics
import math
import ast
import sys
import re


def eval_table(table):
    """Evaluate table elements using literal_eval."""
    for col in table.columns:
        table[col] = table[col].apply(lambda x: ast.literal_eval(str(x)))
    return table
    
    
# Source: https://www.geeksforgeeks.org/python-program-to-flatten-a-nested-list-using-recursion/
def flattenList(nestedList):
    """Flatten a nested list using recursion."""
    if not(bool(nestedList)):
        return nestedList

    if isinstance(nestedList[0], list):
        return flattenList(*nestedList[:1]) + flattenList(nestedList[1:])
    
    return nestedList[:1] + flattenList(nestedList[1:])


def read_data(url, delimiter='\t', duplicates=False):
    """Read data from a file with given delimiter."""
    data = pd.read_csv(
        url,
        sep=delimiter,
        index_col=0,
        encoding='utf-8',
        encoding_errors='backslashreplace'
    )
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
    """Check data conversion for any formatting issues."""
    print("Checking {} data conversion:".format(data_type))
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


def format_data(data, reference):
    """Format data for processing."""
    data_output = data.copy()
    data_output.index = data_output.index.str.upper()
    reference_ids = list(reference.values())
    data_output.index = data_output.index.to_series().replace(reference)
    data_unmapped = data_output[~data_output.index.isin(reference_ids)]
    return data_output, data_unmapped


def output_unmapped(data, url, delimiter='\t'):
    """Output unmapped entities for user information."""
    if len(data.index.tolist()):
        data.to_csv(
            url[:-4] + '_unmapped.txt',
            sep=delimiter)


def extract_data(data):
    """Separate out fold change and p-values."""
    data_c = data.copy()
    data_c = data_c.dropna(axis=0)
    data_c.index = [d.lstrip().rstrip() for d in data_c.index.tolist()]
    _values = data_c.T[::2].T
    _stats = data_c.T[1::2].T
    _values.columns = [x for x in range(len(_values.columns))]
    _stats.columns = [x for x in range(len(_stats.columns))]
    return _values, _stats


def broadcast_transcriptomics(transcriptomics, transcriptomics_stats, gene_dictionary, protein_dictionary):
    """Broadcast transcriptomics data to proteomics."""
    proteomics = transcriptomics.copy()
    proteomics_stats = transcriptomics_stats.copy()

    uniprot_dict = {}
    for x in protein_dictionary.keys():
        uniprot_dict[protein_dictionary[x]] = x

    proteomics.index = proteomics.index.to_series().replace(gene_dictionary)
    proteomics_stats.index = proteomics_stats.index.to_series().replace(gene_dictionary)
    proteomics.index = proteomics.index.to_series().replace(uniprot_dict)
    proteomics_stats.index = proteomics_stats.index.to_series().replace(uniprot_dict)

    return proteomics, proteomics_stats


def copy_columns(data, stats, _max):
    """Copy data columns for every time point provided by the user."""
    data_c = data.copy()
    stats_c = stats.copy()

    counter = 1
    while counter < _max:
        data_c[counter] = data_c[0]
        stats_c[counter] = stats_c[0]
        counter += 1

    return data_c, stats_c


def convert_to_float(lists):
    """Recursively force elements of nested list to float."""
    return [float(el) if not isinstance(el, list) else convert_to_float(el) for el in lists]
                
                
def catenate_data(array):
    """Combine data - average duplicates, remove NAs."""
    tables = [eval_table(x) for x in array]
    combined = pd.concat(tables)
    combined = combined.dropna(axis=0)
    combined = combined.sort_index()

    # Check that types are the same (p-values or confidence intervals)
    if type(combined.iloc[0, 0]) == list:
        if len(set(flattenList(combined.applymap(type).values.tolist()))) != 1:
            raise Exception("Input data types do not match. Please check that all fold change and statistical value types match between datasets.")

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


def prepare_data(network, transcriptomics_url, proteomics_url, metabolomics_url, database_source='reactome'):
    """Get user data and preprocess."""
    should_transcriptomics_exit = False
    should_proteomics_exit = False
    should_metabolomics_exit = False

    # Process transcriptomics
    if transcriptomics_url.lower() != 'none':
        transcriptomics = read_data(url=transcriptomics_url)
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
        proteomics = read_data(url=proteomics_url)
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
        metabolomics = read_data(url=metabolomics_url)
        should_metabolomics_exit = check_data(metabolomics, data_type="Metabolomics")
        if not should_metabolomics_exit:
            metabolomics, metabolomics_stats = extract_data(
                data=metabolomics)
            metabolomics_length = len(metabolomics.columns.tolist())
    else:
        metabolomics_length = 0

    # Check if should exit if encountered invalid value in input data
    if should_transcriptomics_exit or should_proteomics_exit or should_metabolomics_exit:
        print("""
----------------------------------------------------------------
    Formatting/conversion error(s) in data. Please check logs, 
            remedy problematic data points, and rerun.            
----------------------------------------------------------------
""")
        sys.exit(1)

    # Check for broadcasting
    if proteomics_url.lower() == 'none' and transcriptomics_url.lower() != 'none':
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
        if transcriptomics_url.lower() != 'none' and len(transcriptomics.columns.tolist()) != _max and len(transcriptomics.columns.tolist()) == 1:
            transcriptomics, transcriptomics_stats = copy_columns(
                data=transcriptomics,
                stats=transcriptomics_stats,
                _max=_max)
        if proteomics_url.lower() != 'none' and len(proteomics.columns.tolist()) != _max and len(proteomics.columns.tolist()) == 1:
            proteomics, proteomics_stats = copy_columns(
                data=proteomics,
                stats=proteomics_stats,
                _max=_max)
        if metabolomics_url.lower() != 'none' and len(metabolomics.columns.tolist()) != _max and len(metabolomics.columns.tolist()) == 1:
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
    data = catenate_data(array=data_array)
    stats = catenate_data(array=stats_array)

    unmapped = {}
    try:
        unmapped['transcriptomics_unmapped'] = transcriptomics_unmapped.index.tolist()
    except:
        unmapped['transcriptomics_unmapped'] = []
    try:
        unmapped['proteomics_unmapped'] = proteomics_unmapped.index.tolist()
    except:
        unmapped['proteomics_unmapped'] = []

    return data, stats, unmapped


def remove_nulls(values):
    if all(None in v for v in values):
        return []
    else:
        for v in values:
            if None in v:
                values.remove(v)
        return values


def infer_protein_values(values, length):
    protein_vals = []
    for i in range(length):
        pos = []
        for j in range(len(values)):
            if values[j][i] is not None:
                pos.append(values[j][i])
        protein_vals.append(statistics.median(pos))
    return protein_vals


def infer_protein_stats(stats, length, stat_type="float"):
    if stat_type == "array":
        protein_stats = [None for x in range(length)]
    else:
        protein_stats = []
        for i in range(length):
            pos = []
            for j in range(len(stats)):
                if stats[j][i] is not None:
                    pos.append(stats[j][i])
            if len(pos) == 1:
                this_stat = pos[0]
            else:
                this_stat = (math.e * gmean(pos))
            if this_stat > 1.0:
                this_stat = 1.0
            protein_stats.append(this_stat)
    return protein_stats


def prepare_mapping_data(graph, data, stats):
    n = len(data.columns.tolist())
    data_renamed, stats_renamed = reindex_data(data, stats)
    data_max = abs(data_renamed).max().max()
    
    if type(stats_renamed.iloc[0,0]) == list:
        stats_logged = -1 
        stats_max = -1 
    else:
        stats_logged = -1 * np.log10(stats_renamed + 1e-100)
        stats_max = abs(stats_logged).max().max()

    temp_idx = [''.join(c.lower() for c in str(i) if c.isalnum()) for i in data_renamed.index.tolist()]
    temp_idx_set = set(temp_idx)

    chebi_mapping = {}
    for current_id in list(graph.nodes()):
        map_id = graph.nodes()[current_id]['map_id']
        name = graph.nodes()[current_id]['name']
        if 'chebi' in map_id.lower():
            if name in chebi_mapping:
                chebi_mapping[name].add(map_id)
            else:
                chebi_mapping[name] = set()
                chebi_mapping[name].add(map_id)

    return data_renamed, stats_renamed, data_max, stats_max, n, temp_idx, temp_idx_set, chebi_mapping


def reindex_data(data, stats):
    data_renamed = data.copy()
    data_renamed = data_renamed.loc[data_renamed.dropna(axis=0).index.drop_duplicates(keep=False)]
    d_cols = data_renamed.columns
    data_renamed[d_cols] = data_renamed[d_cols].apply(pd.to_numeric, errors='coerce')

    stats_renamed = stats.copy()
    stats_renamed = stats_renamed.loc[stats_renamed.dropna(axis=0).index.drop_duplicates(keep=False)]
    s_cols = stats_renamed.columns
    
    if type(stats_renamed.iloc[0,0]) != list:
        stats_renamed[s_cols] = stats_renamed[s_cols].apply(pd.to_numeric, errors='coerce')

    if len(data_renamed.index.tolist()) != len(data.index.tolist()) or len(stats_renamed.index.tolist()) != len(stats.index.tolist()):
        print('Warning: Duplicate row names were found. All duplicate data were \nremoved from downstream processing. Choose one row per name to \nmaintain the duplicated entity in downstream processing.')
        merge = list(set(data.index.tolist() + stats.index.tolist()))
        merge_after = list(set(data_renamed.index.tolist() + stats_renamed.index.tolist()))
        print("\n Unmapped entities:")
        for x in [y for y in merge if y not in merge_after]:
            print("\t-> " + str(x))

    return data_renamed, stats_renamed
