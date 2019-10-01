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
import pickle
import pandas as pd
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

cmap = matplotlib.cm.get_cmap('RdYlBu')

"""Import internal dependencies
"""
from metaboverse.metaboverse_curate.__main__ import __main__ as curate

from metaboverse.metaboverse_analyze.utils import file_path
from metaboverse.metaboverse_analyze.utils import check_suffix
from metaboverse.metaboverse_analyze.utils import add_data
from metaboverse.metaboverse_analyze.utils import format_data
from metaboverse.metaboverse_analyze.utils import format_times
from metaboverse.metaboverse_analyze.utils import sort_columns
from metaboverse.metaboverse_analyze.utils import even_spacing
from metaboverse.metaboverse_analyze.utils import ratio_spacing
from metaboverse.metaboverse_analyze.utils import retrieve_pathways
from metaboverse.metaboverse_analyze.utils import map_ids

from metaboverse.metaboverse_analyze.graph import create_graph_output
from metaboverse.metaboverse_analyze.graph import add_expression_average
from metaboverse.metaboverse_analyze.graph import add_expression_value
from metaboverse.metaboverse_analyze.graph import add_component_node_edge
from metaboverse.metaboverse_analyze.graph import fetch_analyte_info
from metaboverse.metaboverse_analyze.graph import add_node_edge
from metaboverse.metaboverse_analyze.graph import add_reaction_node
from metaboverse.metaboverse_analyze.graph import process_reactions
from metaboverse.metaboverse_analyze.graph import layout_graph
from metaboverse.metaboverse_analyze.graph import plot_graph
from metaboverse.metaboverse_analyze.graph import build_graph
from metaboverse.metaboverse_analyze.graph import output_graph

"""Set globals
"""
__path__ = os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/Desktop/Metaboverse/tests/analysis_tests/'

def test_file_path():

    path = '.'
    full_path = file_path(
        input=path)

    assert '/' in full_path, 'Unable to obtain full path'

def test_check_suffix():

    file_0 = '/path/file.csv'
    file_1 = '/path/file.tsv'
    file_2 = '/path/file.txt'
    file_3 = '/path/file.xxx'

    suf_0 = check_suffix(
        file=file_0)
    assert suf_0 == ',', 'Unable to parse suffix type'

    suf_1 = check_suffix(
        file=file_1)
    assert suf_1 == '\t', 'Unable to parse suffix type'

    suf_2 = check_suffix(
        file=file_2)
    assert suf_2 == '\t', 'Unable to parse suffix type'

    try:
        check_suffix(
            file=file_3)
    except Exception:
        pass
    else:
        raise Exception('Unable to catch exception case')

def test_add_data():

    cols = [
        '10m_treat_rna_rep1',
        '10m_treat_rna_rep2',
        '00m_treat_rna_rep1',
        '00m_treat_rna_rep2']
    truth = [
        ['Q0010',	75,	56,	73,	85],
        ['Q0017',	18,	18,	30,	17],
        ['Q0032',	158,	174,	185,	148],
        ['Q0045',	19739,	20298,	20750,	21270],
        ['Q0050',	57969,	58509,	67234,	67894]]
    data_truth = pd.DataFrame(truth)
    data_truth.set_index(0, inplace=True)
    del data_truth.index.name
    data_truth.columns = cols

    data = add_data(
        file=__path__ + 'rnaseq_test.txt')
    assert data.head().equals(data_truth)

    return data

def test_format_data(
        data):

    metadata = add_data(
        file=__path__ + 'rnaseq_test_meta.txt')

    formatted = format_data(
        data=data,
        metadata=metadata)

    assert formatted.at['Q0045', '10m'] == ((data.at['Q0045', '10m_treat_rna_rep1'] + data.at['Q0045', '10m_treat_rna_rep2']) / 2), 'Failure in formatting table'

    assert formatted.at['Q0050', '00m'] == ((data.at['Q0050', '00m_treat_rna_rep1'] + data.at['Q0050', '00m_treat_rna_rep2']) / 2), 'Failure in formatting table'

    return formatted

def test_format_times():

    name_0 = 't00'
    name_1 = 't10'
    name_2 = '10t'
    name_3 = 'hello'

    out_0 = format_times(
        name=name_0)
    assert out_0 == 0, 'Failed formatting timepoint name'

    out_1 = format_times(
        name=name_1)
    assert out_1 == 10, 'Failed formatting timepoint name'

    out_2 = format_times(
        name=name_2)
    assert out_2 == 10, 'Failed formatting timepoint name'

    try:
        out_3 = format_times(
            name=name_3)
    except Exception:
        pass
    else:
        raise Exception('Failed to catch a name without numerical values')

def test_sort_columns(
        data):

    truth = [
        't00m_treat_rna_rep1',
        't00m_treat_rna_rep2',
        't10m_treat_rna_rep1',
        't10m_treat_rna_rep2']

    data = sort_columns(data)

    assert truth == data.columns.tolist(), 'Unable to sort columns'

def test_even_spacing(
        data):

    data['30m'] = data['10m'] + 3
    data_c0 = data.copy()
    data_c1 = data.copy()

    total_timepoints = 8
    cols_numeric = [format_times(name) for name in data_c0.columns.tolist()]
    tp_numbers = len(data_c0.columns)
    tp_total = total_timepoints - tp_numbers
    data_c0.columns = cols_numeric

    data_spaced_0 = even_spacing(
        data=data_c0,
        cols_numeric=cols_numeric,
        tp_numbers=tp_numbers,
        tp_total=tp_total)

    total_timepoints = 7
    cols_numeric = [format_times(name) for name in data_c1.columns.tolist()]
    tp_numbers = len(data_c1.columns)
    tp_total = total_timepoints - tp_numbers
    data_c1.columns = cols_numeric

    data_spaced_1 = even_spacing(
        data=data_c1,
        cols_numeric=cols_numeric,
        tp_numbers=tp_numbers,
        tp_total=tp_total)

    assert data_spaced_0.shape == data_spaced_1.shape, 'Failed to properly round number of columns to be created'

    cols_0 = data_spaced_0.columns.tolist()
    cols_0 = [round(x, 1) for x in cols_0]

    cols_1 = data_spaced_1.columns.tolist()
    cols_1 = [round(x, 1) for x in cols_1]

    assert set(cols_0) == set(cols_1) == {0.0, 3.3, 6.7, 10.0, 16.7, 23.3, 30.0}, 'Unable to space '

    data_sub = data_spaced_0.head()
    data_sub.columns = [str(round(x, 1)) for x in data_sub.columns.tolist()]
    assert data_sub.at['Q0045', '23.3'] - data_sub.at['Q0045', '16.7'] == 1, 'Failed calculating smoothened datapoints for even placed times'
    assert data_sub.at['Q0045', '6.7'] - data_sub.at['Q0045', '3.3'] != 1, 'Failed calculating smoothened datapoints for even placed times'

def test_ratio_spacing(
        data):

    data['30m'] = data['10m'] + 3
    data_c0 = data.copy()
    data_c1 = data.copy()

    total_timepoints = 14
    cols_numeric = [format_times(name) for name in data_c0.columns.tolist()]
    tp_numbers = len(data_c0.columns)
    tp_total = total_timepoints - tp_numbers
    data_c0.columns = cols_numeric

    data_ratio_0 = ratio_spacing(
        data=data_c0,
        cols_numeric=cols_numeric,
        tp_numbers=tp_numbers,
        tp_total=tp_total)

    total_timepoints = 7
    cols_numeric = [format_times(name) for name in data_c0.columns.tolist()]
    tp_numbers = len(data_c0.columns)
    tp_total = total_timepoints - tp_numbers
    data_c0.columns = cols_numeric

    data_ratio_1 = even_spacing(
        data=data_c1,
        cols_numeric=cols_numeric,
        tp_numbers=tp_numbers,
        tp_total=tp_total)

def test_create_graph_output():

    graph_path = create_graph_output(
        output=__path__)

    assert graph_path == __path__ + 'graphs/'
    os.rmdir(__path__ + 'graphs/')

def test_add_expression_average():

    data = {
        'Hello': 15,
        'World': None,
        'What':10}

    G = nx.Graph()
    G.add_nodes_from(['Hello', 'World', 'Is'])

    data

    for x in G.nodes():
        if x in data.keys():
            G, value = add_expression_average(
                G,
                100,
                data[x],
                x)

    assert G.nodes()['World']['color'] == 'white', 'Failed to add average expression value for None'
    assert G.nodes()['Hello']['color'] != 'white', 'Failed to add average expression value'
    assert 'Is' in G.nodes().keys(), 'Failed at average expression value addition'
    assert G.nodes()['Is'] == {}, 'Failed at average expression value addition'
    assert 'What' not in G.nodes().keys(), 'Failed at average expression value addition'

    G.nodes()['World']

def test_add_expression_value():

    data = pd.DataFrame()
    data[0] = [15, None, 10]
    data.index = ['Hello','World','What']

    G = nx.Graph()
    G.add_nodes_from(['Hello', 'World', 'Is'])

    for x in G.nodes():
        if x in data.index.tolist():
            G, value = add_expression_value(
                graph=G,
                data=data,
                analyte=x,
                analyte_id=x)

    assert G.nodes()['World']['color'] == 'white', 'Failed to add expression value for None'
    assert G.nodes()['Hello']['color'] != 'white', 'Failed to add expression value'

def test_add_component_node_edge(
        network):

    network = reference

    data = pd.DataFrame()
    data[0] = [15, None, 10]
    data.index = ['Hello','World','What']

    G = nx.Graph()


def test_add_reaction_node():

    G = nx.Graph()
    G.add_nodes_from(['reactions', 'other'])



def __main__():

    # Utils tests
    test_file_path()
    test_check_suffix()
    data = test_add_data()
    formatted_data = test_format_data(
        data=data)
    test_format_times()
    test_sort_columns(
        data=data)
    #test_even_spacing( # Run more tests with this
    #    data=formatted_data)
    #test_ratio_spacing( # Finish these tests
    #    data=formatted_data)

    # Graph tests
    test_create_graph_output()
    test_add_expression_average()
    test_add_expression_value()

    args_dict = {
        'species': 'HSA',
        'output': __path__}
    curate(
        args_dict=args_dict)
    with open(__path__ + 'HSA_metaboverse_db.pickle', 'rb') as network_file:
        reference = pickle.load(network_file)

    test_add_component_node_edge(
        network=reference)








    os.remove(__path__ + 'HSA_metaboverse_db.pickle')
