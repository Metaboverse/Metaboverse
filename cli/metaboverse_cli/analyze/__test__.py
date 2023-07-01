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
import zipfile
import importlib.util
import numpy as np
import networkx as nx
import xml.etree.ElementTree as et
import pandas as pd
import pickle
import os

"""prepare_data.py
"""
print("Testing prepare_data.py")
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
prepare_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare_data)
read_data = prepare_data.read_data
format_data = prepare_data.format_data
output_unmapped = prepare_data.output_unmapped
extract_data = prepare_data.extract_data
broadcast_transcriptomics = prepare_data.broadcast_transcriptomics
copy_columns = prepare_data.copy_columns
catenate_data = prepare_data.catenate_data
__main__ = prepare_data.__main__

zipped_net = os.path.abspath(
    './metaboverse_cli/analyze/test/HSA.zip')
with zipfile.ZipFile(zipped_net, 'r') as zip_file:
    zip_file.extractall(
        os.path.abspath(
            './metaboverse_cli/analyze/test'))

network_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/HSA.mvdb")
with open(network_url, 'rb') as network_file:
    network = pickle.load(network_file)

transcriptomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/transcriptomics.txt")
proteomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/proteomics.txt")
metabolomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/metabolomics.txt")

# read_data()
print("Testing read_data()")
transcriptomics_df = read_data(
    url=transcriptomics_url)
assert len(transcriptomics_df.columns.tolist()) == 2, "read_data() failed"

# format_data()
print("Testing format_data()")
e_sym = {}
for k, v in network['ensembl_synonyms'].items():
    e_sym[v] = k
    e_sym[k] = k
t_mapped, t_unmapped = format_data(
    data=transcriptomics_df,
    reference=e_sym)
assert t_mapped.shape[0] >= t_unmapped.shape[0], "format_data()"
assert t_mapped.shape[1] == 2, "format_data()"
assert t_unmapped.shape[1] == 2, "format_data()"
assert transcriptomics_df.loc['GCN1'][0] == t_mapped.loc['ENST00000300648'][0], 'format_data() failed'
try:
    t_mapped.loc['GCN1']
except KeyError:
    pass
else:
    raise Exception('format_data() failed')

# output_unmapped()
print("Testing output_unmapped()")
output_unmapped(
    data=t_unmapped,
    url=transcriptomics_url)

# extract_data()
print("Testing extract_data()")
_v, _s = extract_data(
    data=t_mapped)
assert _v.shape[1] == 1, "extract_data()"
assert _s.shape[1] == 1, "extract_data()"

# broadcast_transcriptomics()
# How is this mapped gene vs protein in graphing?
print("Testing broadcast_transcriptomics()")
_p, _p_stats = broadcast_transcriptomics(
    transcriptomics=_v,
    transcriptomics_stats=_s,
    gene_dictionary=network['ensembl_synonyms'],
    protein_dictionary=network['uniprot_synonyms'])
assert _v.loc['ENST00000300648'][0] == _p.loc['Q92616'][0], 'broadcast_transcriptomics() failed'

# copy_columns()
print("Testing copy_columns()")
data_col, stat_col = copy_columns(_v, _s, 6)
assert len(_v.columns.tolist()) == 1, 'copy_columns() failed'
assert len(data_col.columns.tolist()) == 6, 'copy_columns() failed'
assert data_col[0].values.all(
) == data_col[5].values.all(), 'copy_columns() failed'
assert stat_col[0].values.all(
) == stat_col[3].values.all(), 'copy_columns() failed'

# catenate_data()
print("Testing catenate_data()")
concat_df = catenate_data([_v, _p])
assert concat_df.shape > _v.shape and concat_df.shape > _p.shape, 'catenate_data() failed'

# Test main()
print("Testing main()")
data, stats, unmapped = __main__(
    network=network,
    transcriptomics_url=transcriptomics_url,
    proteomics_url=proteomics_url,
    metabolomics_url=metabolomics_url)
assert data.shape > _v.shape, '__main__() from prepare_data.py failed'

try:
    stats.loc['UBE2N']
except KeyError:
    pass
else:
    raise Exception('__main__() from prepare_data.py failed')

try:
    data.loc['P61088']
except KeyError:
    raise Exception('__main__() from prepare_data.py failed')

try:
    data.loc['ENST00000318066']
except KeyError:
    raise Exception('__main__() from prepare_data.py failed')

try:
    data.loc['WRAP73']
except KeyError:
    raise Exception('__main__() from prepare_data.py failed')


"""utils.py
"""
print("Testing utils.py")
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/utils.py"))
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)
file_path = utils.file_path
check_suffix = utils.check_suffix
add_data = utils.add_data
convert_rgba = utils.convert_rgba

# file_path()
print("Testing file_path()")
file = "./metaboverse_cli/analyze/utils.py"
assert file_path(file) == os.path.abspath(file), 'file_path() failed'

# check_suffix()
print("Testing check_suffix()")
assert check_suffix('file.txt') == '\t', 'check_suffix() failed'

# add_data()
print("Testing add_data()")
df = add_data(transcriptomics_url)
assert type(df) == pd.DataFrame, 'add_data() failed'

# convert_rgba()
print("Testing convert_rgba()")
n = 2
reaction_color = (0.75, 0.75, 0.75, 1)
missing_color = (1, 1, 1, 1)
color1 = [reaction_color for x in range(n)]
assert convert_rgba(color1) == [(191, 191, 191, 1),
                                (191, 191, 191, 1)], 'convert_rgba() failed'
n = 1
color2 = [missing_color for x in range(n)]
assert convert_rgba(color2) == [(255, 255, 255, 1)], 'convert_rgba() failed'

"""model.py
"""
print("Testing model.py")
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/model.py"))
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
name_graph = model.name_graph
build_graph = model.build_graph
process_reactions = model.process_reactions
add_node_edge = model.add_node_edge
check_complexes = model.check_complexes
uniprot_ensembl_reference = model.uniprot_ensembl_reference
map_attributes = model.map_attributes
extract_value = model.extract_value
output_graph = model.output_graph
compile_pathway_degree = model.compile_pathway_degree
compile_node_degrees = model.compile_node_degrees
remove_nulls = model.remove_nulls
infer_protein_values = model.infer_protein_values
infer_protein_stats = model.infer_protein_stats
broadcast_values = model.broadcast_values
make_motif_reaction_dictionary = model.make_motif_reaction_dictionary
load_metabolite_synonym_dictionary = model.load_metabolite_synonym_dictionary

G = nx.DiGraph()
G.add_node('Alpha')
G.nodes()['Alpha']['name'] = 'Alpha'
G.nodes()['Alpha']['map_id'] = 'alpha'
G.nodes()['Alpha']['type'] = 'reactant'
G.nodes()['Alpha']['sub_type'] = 'gene'
G.nodes()['Alpha']['complex'] = 'false'
G.nodes()['Alpha']['values'] = [5]
G.nodes()['Alpha']['stats'] = [0.5]
G.add_node('Beta')
G.nodes()['Beta']['name'] = 'Beta'
G.nodes()['Beta']['map_id'] = 'beta'
G.nodes()['Beta']['type'] = 'reaction'
G.nodes()['Beta']['sub_type'] = 'reaction'
G.nodes()['Beta']['complex'] = 'false'
G.nodes()['Beta']['values'] = [None]
G.nodes()['Beta']['stats'] = [None]
G.add_node('Gamma')
G.nodes()['Gamma']['name'] = 'Gamma'
G.nodes()['Gamma']['map_id'] = 'gamma'
G.nodes()['Gamma']['type'] = 'product'
G.nodes()['Gamma']['sub_type'] = 'product'
G.nodes()['Gamma']['complex'] = 'true'
G.nodes()['Gamma']['values'] = [None]
G.nodes()['Gamma']['stats'] = [None]
G.add_node('Delta')
G.nodes()['Delta']['name'] = 'Delta'
G.nodes()['Delta']['map_id'] = 'delta'
G.nodes()['Delta']['type'] = 'complex_component'
G.nodes()['Delta']['sub_type'] = 'protein_component'
G.nodes()['Delta']['complex'] = 'false'
G.nodes()['Delta']['values'] = [None]
G.nodes()['Delta']['stats'] = [None]
G.add_node('Epsilon')
G.nodes()['Epsilon']['name'] = 'Epsilon'
G.nodes()['Epsilon']['map_id'] = 'epsilon'
G.nodes()['Epsilon']['type'] = 'complex_component'
G.nodes()['Epsilon']['sub_type'] = 'protein_component'
G.nodes()['Epsilon']['complex'] = 'false'
G.nodes()['Epsilon']['values'] = [9]
G.nodes()['Epsilon']['stats'] = [.9]
G.add_edges_from([
    ('Alpha', 'Beta')])
G.add_edges_from([
    ('Beta', 'Gamma')])
G.add_edges_from([
    ('Gamma', 'Alpha')])
G.add_edges_from([
    ('Alpha', 'Delta')])
G.add_edges_from([
    ('Gamma', 'Epsilon')])

net_test = {
    'reaction_database': {
        'reaction_1': {
            'compartment': 'compartment_1',
            'id': 'reaction_1',
            'name': 'test reaction',
            'reversible': 'true',
            'notes': 'This is a test',
            'reactants': ['species_1', 'species_2'],
            'products': ['species_3'],
            'modifiers': [['species_4', 'catalyst']]
        },
        'reaction_2': {
            'compartment': 'compartment_1',
            'id': 'reaction_2',
            'name': 'test reaction 2',
            'reversible': 'true',
            'notes': 'This is a test 2',
            'reactants': ['species_5', 'species_6'],
            'products': ['species_7'],
            'modifiers': [['species_8', 'inhibitor']]
        }
    },
    'pathway_database': {
        'R-HSA-1': {
            'id': 'pathway_1',
            'reactome': 'R-HSA-1',
            'name': 'Test pathway',
            'reactions': ['reaction_1', 'reaction_2']
        },
        'R-HSA-2': {
            'id': 'pathway_2',
            'reactome': 'R-HSA-2',
            'name': 'Test pathway',
            'reactions': ['reaction_1', 'reaction_2', 'reaction_3', 'reaction_4']
        }
    }
}

test_args = {
    'output': os.path.abspath("./metaboverse_cli/analyze/test"),
    'output_file': "test.mvrs",
    'bad_output_file': os.path.abspath("./metaboverse_cli/analyze/test/"),
    'organism_id': "HSA",
    'network': network_url
}

# name_graph()
print("Testing name_graph()")
name = name_graph(
    output_file=test_args['output_file'],
    species_id=test_args['organism_id']
)
assert name == test_args['output_file'], 'name_graph() failed'

name = name_graph(
    output_file=test_args['bad_output_file'],
    species_id=test_args['organism_id']
)
assert name == test_args['organism_id'] + '_global_reactions.mvrs', 'name_graph() failed'

# build_graph()
print("Testing build_graph()")
net_copy = net_test.copy()
gg1, net_copy, pathway_database = build_graph(
    args_dict=test_args,
    network=net_copy['reaction_database'],
    pathway_database=net_copy['pathway_database'],
    species_reference={},
    name_reference={},
    protein_reference={},
    chebi_dictionary={},
    uniprot_reference={},
    complexes={},
    species_id={},
    gene_reference={},
    compartment_reference={'compartment_1': 'hello'},
    component_database={
        'species_1': {'is': 'gene0', 'name': 'gene0', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_2': {'is': 'gene0', 'name': 'gene0', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_3': {'is': 'gene0', 'name': 'gene0', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_4': {'is': 'gene0', 'name': 'gene0', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_5': {'is': 'gene1', 'name': 'geneA', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_6': {'is': 'gene2', 'name': 'geneB', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_7': {'is': 'gene3', 'name': 'geneC', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_8': {'is': 'gene4', 'name': 'geneD', 'compartment': 'none', 'type': 'gene', 'hasPart': []}
    }
)
assert list(gg1.nodes()) == [
    'reaction_1',
    'species_1',
    'species_2',
    'species_3',
    'species_4',
    'reaction_2',
    'species_5',
    'species_6',
    'species_7',
    'species_8'], 'build_graph() failed'

# process_reactions()
print("Testing process_reactions()")
net_copy = net_test['reaction_database'].copy()
gg2 = nx.DiGraph()
key_hash = set()
remove_keys = []
gg2, net_copy, key_hash, remove_keys = process_reactions(
    graph=gg2,
    reactome_id='reaction_2',
    network=net_copy,
    species_reference={},
    name_reference={},
    protein_reference={},
    chebi_dictionary={},
    uniprot_reference={},
    complex_reference={},
    species_id={},
    gene_reference={},
    compartment_reference={'compartment_1': 'hello'},
    component_database={
        'species_5': {'is': 'gene1', 'name': 'geneA', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_6': {'is': 'gene2', 'name': 'geneB', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_7': {'is': 'gene3', 'name': 'geneC', 'compartment': 'none', 'type': 'gene', 'hasPart': []},
        'species_8': {'is': 'gene4', 'name': 'geneD', 'compartment': 'none', 'type': 'gene', 'hasPart': []}
    },
    key_hash=key_hash,
    remove_keys=remove_keys)
try:
    gg2.nodes()['reaction_2']
    gg2.nodes()['species_5']
    gg2.edges()[('species_5', 'reaction_2')]
    gg2.edges()[('species_6', 'reaction_2')]
    gg2.edges()[('reaction_2', 'species_7')]
except:
    raise Exception('process_reactions() failed')

# add_node_edge()
print("Testing add_node_edge()")
G_add = G.copy()
G_add = add_node_edge(
    graph=G_add,
    id='thisisnew',
    map_id='mapper',
    name='test',
    compartment='none',
    reaction_membership='Gamma',
    type='complex_component',
    sub_type='none',
    reversible='false',
    complex_reference={},
    species_reference={},
    name_reference={},
    protein_reference={},
    compartment_reference={})
try:
    G_add.nodes()['thisisnew']
    G_add.edges()[('thisisnew', 'Gamma')]
except:
    raise Exception('add_node_edge() failed')

# check_complexes()
print("Testing check_complexes()")
G_complex = G.copy()
G_complex, add_components = check_complexes(
    species_id='HSA',
    graph=G_complex,
    id='Epsilon',
    complex_reference={'Epsilon': ['x', 'y', 'z']},
    species_reference={''},
    name_reference={},
    protein_reference={},
    chebi_dictionary={},
    uniprot_reference={},
    gene_reference={},
    component_database={'Epsilon': {'hasPart': ['x', 'y', 'z']}},
    compartment_reference={})
try:
    G_complex.nodes()['x']
    G_complex.edges()[('x', 'Epsilon')]
    G_complex.edges()[('y', 'Epsilon')]
    G_complex.edges()[('z', 'Epsilon')]
except:
    raise Exception('check_complexes() failed')

# uniprot_ensembl_reference()
print("Testing uniprot_ensembl_reference()")
uni_ref = {'Alpha': 'Beta', 'Epsilon': 'F'}
ens_ref = {'Beta': 'Delta'}
ref_test = uniprot_ensembl_reference(
    uniprot_reference=uni_ref,
    ensembl_reference=ens_ref)
assert ref_test == {'Alpha': 'Delta'}, 'uniprot_ensembl_reference() failed'

# compile_pathway_degree()
print("Testing compile_pathway_degree()")
s_p = compile_pathway_degree(
    pathways=net_test['pathway_database'],
    scale_factor=2
)
assert len(list(s_p.keys())) == 1, 'compile_pathway_degree() failed'

# compile_node_degrees()
print("Testing compile_node_degrees()")
d_d = {
    'Alpha': 3,
    'Beta': 2,
    'Gamma': 3,
    'Delta': 1,
    'Epsilon': 1
}
degree_dictionary = compile_node_degrees(
    graph=G)
assert degree_dictionary == d_d, 'compile_node_degrees() failed'

# map_attributes()
print("Testing map_attributes()")
G_map = G.copy()
data = pd.DataFrame()
data[0] = [1, 3, 5]
data.index = ['Alpha', 'Gamma', 'epsilon']
stats = pd.DataFrame()
stats[0] = [.1, .3, .5]
stats.index = ['Alpha', 'Gamma', 'epsilon']

G_mapped, data_max, stats_max, non_mappers = map_attributes(
    args_dict={},
    graph=G_map,
    data=data,
    stats=stats,
    name_reference={
        'Alpha': 'Alpha',
        'Beta': 'Beta',
        'Gamma': 'Gamma',
        'Delta': 'Delta',
        'epsilon': 'Epsilon',
    },
    degree_dictionary=degree_dictionary,
    chebi_dictionary={
        'Alpha': 'Alpha',
        'Beta': 'Beta',
        'Gamma': 'Gamma',
        'Delta': 'Delta',
        'epsilon': 'epsilon',
    },
    chebi_synonyms={},
    metabolite_mapper={
        'mapping_dictionary': {
            'Alpha': 'Alpha',
            'Beta': 'Beta',
            'Gamma': 'Gamma',
            'Delta': 'Delta',
            'epsilon': 'epsilon',
        },


    },
    uniprot_mapper={
        'Alpha': 'Alpha',
        'Beta': 'Beta',
        'Gamma': 'Gamma',
        'Delta': 'Delta',
        'epsilon': 'epsilon',
    })

assert data_max == 5, 'map_attributes() failed'
assert stats_max == 1.0, 'map_attributes() failed'
assert non_mappers == [], 'map_attributes() failed'
assert G_mapped.nodes()['Alpha']['values'] == [1], 'map_attributes() failed'
assert G_mapped.nodes()['Alpha']['stats'] == [0.1], 'map_attributes() failed'
assert G_mapped.nodes()['Beta']['values'] == [None], 'map_attributes() failed'
assert G_mapped.nodes()['Beta']['stats'] == [None], 'map_attributes() failed'
assert G_mapped.nodes()['Gamma']['values'] == [3], 'map_attributes() failed'
assert G_mapped.nodes()['Gamma']['stats'] == [0.3], 'map_attributes() failed'
assert G_mapped.nodes()['Delta']['values'] == [None], 'map_attributes() failed'
assert G_mapped.nodes()['Delta']['stats'] == [None], 'map_attributes() failed'
assert G_mapped.nodes()['Epsilon']['values'] == [5], 'map_attributes() failed'
assert G_mapped.nodes()['Epsilon']['stats'] == [0.5], 'map_attributes() failed'

# extract_value()
"""
v_e = [(0.0, 0.0, 0.3, 1.0),
 (1.0, 0.9921568627450981, 0.9921568627450981, 1.0),
 (0.5, 0.0, 0.0, 1.0)]
_v = extract_value(
        value_array=[-5,0,5],
        max_value=5,
        type="value")
assert v_e == _v, 'extract_value() failed'
s_e = [(0.403921568627451, 0.0, 0.05098039215686274, 1.0),
 (0.9882352941176471, 0.6664975009611688, 0.5547558631295655, 1.0),
 (1.0, 0.9607843137254902, 0.9411764705882353, 1.0)]
_s = extract_value(
        value_array=[0,0.5,1],
        max_value=1,
        type="stat")
assert s_e == _s, 'extract_value() failed'
"""

# output_graph()
print("Testing output_graph()")
output_graph(
    graph=G,
    output_name=test_args['output_file'],
    pathway_dictionary={},
    collapsed_pathway_dictionary={},
    super_pathways={},
    reaction_dictionary={},
    collapsed_reaction_dictionary={},
    motif_reaction_dictionary={},
    mod_collapsed_pathways={},
    degree_dictionary={},
    max_value=5,
    max_stat=1,
    categories=[0],
    labels=['OO'],
    blocklist=[],
    species_blocklist=[],
    metadata={},
    unmapped={})
if os.path.exists(test_args['output_file']):
    os.remove(str(test_args['output_file']))
else:
    raise Exception('output_graph() failed')

# remove_nulls()
print("Testing remove_nulls()")
vals1 = [[None, 1, 6, None]]
assert remove_nulls(vals1) == [], 'remove_nulls() failed'
vals2 = [[None, 1, 6, None], [1, 2, 3]]
assert remove_nulls(vals2) == [[1, 2, 3]], 'remove_nulls() failed'

# infer_protein_values()
print("Testing infer_protein_values()")
vals = [[1], [2], [3], [3], [4]]
length = 1
assert infer_protein_values(vals, length) == [3], 'infer_protein_values() failed'

# infer_protein_stats()
print("Testing infer_protein_stats()")
assert infer_protein_stats(vals, length) == [1.0], 'infer_protein_stats() failed'

# broadcast_values()
print("Testing broadcast_values()")
G_update = G.copy()
G_update = broadcast_values(
    args_dict={},
    graph=G_update,
    categories=[0],
    max_value=100,
    max_stat=1,
    broadcast_genes=True,
    broadcast_metabolites=False)
assert G_update.nodes()['Alpha'] == G.nodes()[
    'Alpha'], 'broadcast_values() failed'
assert G_update.nodes()['Beta'] == G.nodes()[
    'Beta'], 'broadcast_values() failed'
assert G_update.nodes()['Epsilon'] == G.nodes()[
    'Epsilon'], 'broadcast_values() failed'
assert G_update.nodes()['Gamma']['values'] == [
    9.0], 'broadcast_values() failed to broadcast proteins to complexes'
assert G_update.nodes()['Delta']['stats'] == [
    0.5], 'broadcast_values() failed to broadcast genes to proteins'

# Test casting metabolite values to protein complexes
G_update = G.copy()
G_update.add_node('Zeta')
G_update.nodes()['Zeta']['name'] = 'Zeta'
G_update.nodes()['Zeta']['map_id'] = 'zeta'
G_update.nodes()['Zeta']['type'] = 'metabolite_component'
G_update.nodes()['Zeta']['sub_type'] = 'metabolite_component'
G_update.nodes()['Zeta']['complex'] = 'false'
G_update.nodes()['Zeta']['values'] = [7]
G_update.nodes()['Zeta']['stats'] = [.7]
G_update.add_edges_from([
    ('Zeta', 'Gamma')])

G_update = broadcast_values(
    args_dict={},
    graph=G_update,
    categories=[0],
    max_value=100,
    max_stat=1,
    broadcast_genes=True,
    broadcast_metabolites=True)
assert G_update.nodes()['Alpha'] == G.nodes()[
    'Alpha'], 'broadcast_values() failed'
assert G_update.nodes()['Beta'] == G.nodes()[
    'Beta'], 'broadcast_values() failed'
assert G_update.nodes()['Epsilon'] == G.nodes()[
    'Epsilon'], 'broadcast_values() failed'
assert G_update.nodes()['Gamma']['values'] == [
    8.0], 'broadcast_values() failed to broadcast metabolites to complexes'
assert G_update.nodes()['Delta']['stats'] == [
    0.5], 'broadcast_values() failed to broadcast genes to proteins'

# make_motif_reaction_dictionary()
print("Testing load_motif_reaction_dictionary()")
updated_reactions = {
    'reaction_1_reaction_3': {
        'collapsed': 'true',
        'collapsed_reactions': ['reaction_1', 'reaction_3'],
        'compartment': 'compartment_1',
        'id': 'reaction_1_reaction_3',
        'name': 'test reaction',
        'reversible': 'true',
        'notes': 'This is a test',
        'reactants': ['species_1', 'species_2'],
        'products': ['species_3'],
        'modifiers': [['species_4', 'catalyst']],
        'additional_components': ['ENST1']
    },
    'reaction_2': {
        'compartment': 'compartment_1',
        'id': 'reaction_2',
        'name': 'test reaction 2',
        'reversible': 'true',
        'notes': 'This is a test 2',
        'reactants': ['species_5', 'species_6'],
        'products': ['species_7'],
        'modifiers': [['species_8', 'inhibitor']]
    }
}
updated_pathway_dictionary = {
    'R-HSA-1': {
        'id': 'pathway_1',
        'reactome': 'R-HSA-1',
        'name': 'Test pathway',
        'reactions': ['reaction_1_reaction_3', 'reaction_2']
    }
}
mot_dic = make_motif_reaction_dictionary(
    network=net_test,
    updated_reactions=updated_reactions,
    updated_pathway_dictionary=updated_pathway_dictionary)
assert mot_dic == {'reaction_1_reaction_3': ['pathway_1'], 'reaction_2': [
    'pathway_1']}, 'make_motif_reaction_dictionary() failed'

# load_metabolite_synonym_dictionary()
print("Testing load_metabolite_synonym_dictionary()")
mapper = load_metabolite_synonym_dictionary()
assert list(mapper.keys()) == ['hmdb_dictionary', 'display_dictionary',
                               'mapping_dictionary'], 'load_metabolite_synonym_dictionary() failed'

"""collapse.py
"""
print("Testing collapse.py")
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/collapse.py"))
collapse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(collapse)
generate_updated_dictionary = collapse.generate_updated_dictionary
collapse_nodes = collapse.collapse_nodes

G_collapse = nx.DiGraph()
G_collapse.add_node('N1')
G_collapse.nodes()['N1']['values'] = [5]
G_collapse.nodes()['N1']['stats'] = [0.5]
G_collapse.add_node('N2')
G_collapse.nodes()['N2']['values'] = [5]
G_collapse.nodes()['N2']['stats'] = [0.5]
G_collapse.add_node('N3')
G_collapse.nodes()['N3']['values'] = [5]
G_collapse.nodes()['N3']['stats'] = [0.5]
G_collapse.add_node('N4')
G_collapse.nodes()['N4']['values'] = [None]
G_collapse.nodes()['N4']['stats'] = [None]
G_collapse.add_node('N5')
G_collapse.nodes()['N5']['values'] = [None]
G_collapse.nodes()['N5']['stats'] = [None]
G_collapse.add_node('N6')
G_collapse.nodes()['N6']['values'] = [5]
G_collapse.nodes()['N6']['stats'] = [0.5]
G_collapse.add_node('N7')
G_collapse.nodes()['N7']['values'] = [None]
G_collapse.nodes()['N7']['stats'] = [None]
G_collapse.add_node('N8')
G_collapse.nodes()['N8']['values'] = [None]
G_collapse.nodes()['N8']['stats'] = [None]
G_collapse.add_node('N29')
G_collapse.nodes()['N29']['values'] = [None]
G_collapse.nodes()['N29']['stats'] = [None]
G_collapse.add_node('N30')
G_collapse.nodes()['N30']['values'] = [None]
G_collapse.nodes()['N30']['stats'] = [None]
G_collapse.add_node('N9')
G_collapse.nodes()['N9']['values'] = [5]
G_collapse.nodes()['N9']['stats'] = [0.5]
G_collapse.add_node('N10')
G_collapse.nodes()['N10']['values'] = [None]
G_collapse.nodes()['N10']['stats'] = [None]
G_collapse.add_node('N11')
G_collapse.nodes()['N11']['values'] = [None]
G_collapse.nodes()['N11']['stats'] = [None]
G_collapse.add_node('N12')
G_collapse.nodes()['N12']['values'] = [5]
G_collapse.nodes()['N12']['stats'] = [0.5]
G_collapse.add_node('N13')
G_collapse.nodes()['N13']['values'] = [5]
G_collapse.nodes()['N13']['stats'] = [0.5]
G_collapse.add_node('N14')
G_collapse.nodes()['N14']['values'] = [None]
G_collapse.nodes()['N14']['stats'] = [None]
G_collapse.add_node('N15')
G_collapse.nodes()['N15']['values'] = [None]
G_collapse.nodes()['N15']['stats'] = [None]
G_collapse.add_node('N16')
G_collapse.nodes()['N16']['values'] = [None]
G_collapse.nodes()['N16']['stats'] = [None]
G_collapse.add_node('N17')
G_collapse.nodes()['N17']['values'] = [None]
G_collapse.nodes()['N17']['stats'] = [None]
G_collapse.add_node('N18')
G_collapse.nodes()['N18']['values'] = [None]
G_collapse.nodes()['N18']['stats'] = [None]
G_collapse.add_node('N19')
G_collapse.nodes()['N19']['values'] = [5]
G_collapse.nodes()['N19']['stats'] = [0.5]
G_collapse.add_node('N20')
G_collapse.nodes()['N20']['values'] = [None]
G_collapse.nodes()['N20']['stats'] = [None]
G_collapse.add_node('N21')
G_collapse.nodes()['N21']['values'] = [5]
G_collapse.nodes()['N21']['stats'] = [0.5]
G_collapse.add_node('N22')
G_collapse.nodes()['N22']['values'] = [None]
G_collapse.nodes()['N22']['stats'] = [None]
G_collapse.add_node('N23')
G_collapse.nodes()['N23']['values'] = [None]
G_collapse.nodes()['N23']['stats'] = [None]
G_collapse.add_node('N24')
G_collapse.nodes()['N24']['values'] = [None]
G_collapse.nodes()['N24']['stats'] = [None]
G_collapse.add_node('N25')
G_collapse.nodes()['N25']['values'] = [5]
G_collapse.nodes()['N25']['stats'] = [0.5]
G_collapse.add_node('N26')
G_collapse.nodes()['N26']['values'] = [None]
G_collapse.nodes()['N26']['stats'] = [None]
G_collapse.add_node('N27')
G_collapse.nodes()['N27']['values'] = [None]
G_collapse.nodes()['N27']['stats'] = [None]
G_collapse.add_node('N28')
G_collapse.nodes()['N28']['values'] = [None]
G_collapse.nodes()['N28']['stats'] = [None]
G_collapse.add_node('N100')
G_collapse.nodes()['N100']['values'] = [5]
G_collapse.nodes()['N100']['stats'] = [0.5]

collapser_database = {
    'R1': {
        'compartment': 'none',
        'id': 'R1',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N1'],
        'products': ['N2'],
        'modifiers': [],
        'additional_components': []
    },
    'R2': {
        'compartment': 'none',
        'id': 'R2',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N3'],
        'products': ['N4', 'N5'],
        'modifiers': [],
        'additional_components': []
    },
    'R3': {
        'compartment': 'none',
        'id': 'R3',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N4', 'N5'],
        'products': ['N6'],
        'modifiers': [],
        'additional_components': []
    },
    'R4': {
        'compartment': 'none',
        'id': 'R4',
        'reactome': 'Re4',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N7'],
        'products': ['N8'],
        'modifiers': [],
        'additional_components': []
    },
    'R5': {
        'compartment': 'none',
        'id': 'R5',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N9'],
        'products': ['N10'],
        'modifiers': [],
        'additional_components': []
    },
    'R6': {
        'compartment': 'none',
        'id': 'R6',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N10'],
        'products': ['N11'],
        'modifiers': [],
        'additional_components': []
    },
    'R7': {
        'compartment': 'none',
        'id': 'R7',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N11'],
        'products': ['N12'],
        'modifiers': [],
        'additional_components': []
    },
    'R8': {
        'compartment': 'none',
        'id': 'R8',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N13'],
        'products': ['N14', 'N15'],
        'modifiers': [],
        'additional_components': []
    },
    'R9': {
        'compartment': 'none',
        'id': 'R9',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N14', 'N15'],
        'products': ['N16', 'N17', 'N18'],
        'modifiers': [],
        'additional_components': []
    },
    'R10': {
        'compartment': 'none',
        'id': 'R10',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N16', 'N17', 'N18'],
        'products': ['N19', 'N20'],
        'modifiers': [],
        'additional_components': []
    },
    'R11': {
        'compartment': 'none',
        'id': 'R11',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N21'],
        'products': ['N22'],
        'modifiers': [],
        'additional_components': []
    },
    'R12': {
        'compartment': 'none',
        'id': 'R12',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N22'],
        'products': ['N23'],
        'modifiers': [],
        'additional_components': []
    },
    'R13': {
        'compartment': 'none',
        'id': 'R13',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N23'],
        'products': ['N24'],
        'modifiers': [],
        'additional_components': []
    },
    'R14': {
        'compartment': 'none',
        'id': 'R14',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N25'],
        'products': ['N26', 'N27'],
        'modifiers': [],
        'additional_components': []
    },
    'R15': {
        'compartment': 'none',
        'id': 'R15',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N26', 'N27'],
        'products': ['N28'],
        'modifiers': [['N100', 'inhibitor'], ['N100', 'catalyst']],
        'additional_components': []
    },
    'R16': {
        'compartment': 'none',
        'id': 'R16',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N8'],
        'products': ['N29'],
        'modifiers': [],
        'additional_components': []
    },
    'R17': {
        'compartment': 'none',
        'id': 'R17',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N29'],
        'products': ['N30'],
        'modifiers': [],
        'additional_components': []
    }
}

collapse_degree_dictionary = {
    'R1': 4,
    'R2': 3,
    'R3': 6,
    'R4': 5,
    'R5': 2,
    'R6': 4,
    'R7': 4,
    'R8': 4,
    'R9': 5,
    'R10': 7,
    'R11': 1,
    'R12': 3,
    'R13': 4,
    'R14': 2,
    'R15': 5,
    'R16': 3,
    'R17': 4,
    'N1': 7,
    'N2': 8,
    'N3': 5,
    'N4': 3,
    'N5': 8,
    'N6': 6,
    'N7': 4,
    'N8': 3,
    'N9': 6,
    'N10': 7,
    'N11': 8,
    'N12': 4,
    'N13': 5,
    'N14': 2,
    'N15': 4,
    'N16': 7,
    'N17': 8,
    'N18': 5,
    'N19': 3,
    'N20': 6,
    'N21': 8,
    'N22': 9,
    'N23': 3,
    'N24': 43,
    'N25': 7,
    'N26': 4,
    'N27': 7,
    'N28': 3,
    'N29': 3,
    'N30': 7,
    'N100': 3
}

# collapse_nodes()
print("Testing collapse_nodes()")
G_coll1 = G_collapse.copy()
G_coll1, updated_rxns1, changed_rxns1, removed_rxn1 = collapse_nodes(
    args_dict={},
    graph=G_coll1,
    reaction_dictionary=collapser_database,
    neighbors_dictionary={},
    degree_dictionary=collapse_degree_dictionary,
    samples=1,
    collapse_with_modifiers=False,
    blocklist=[])

final_reactions1 = [
    'R1',
    'R2_R3',
    'R4',
    'R5_R6_R7',
    'R8_R9_R10',
    'R11',
    'R12',
    'R13',
    'R14',
    'R15',
    'R16',
    'R17'
]
assert list(updated_rxns1.keys()
            ) == final_reactions1, 'collapse_nodes() failed'

G_coll2 = G_collapse.copy()
G_coll2, updated_rxns2, changed_rxns2, removed_rxn2 = collapse_nodes(
    args_dict={},
    graph=G_coll2,
    reaction_dictionary=collapser_database,
    neighbors_dictionary={},
    degree_dictionary=collapse_degree_dictionary,
    samples=1,
    collapse_with_modifiers=True,
    blocklist=[])
final_reactions2 = [
    'R1',
    'R2_R3',
    'R4',
    'R5_R6_R7',
    'R8_R9_R10',
    'R11',
    'R12',
    'R13',
    'R14_R15',
    'R16',
    'R17'
]
assert list(updated_rxns2.keys()
            ) == final_reactions2, 'collapse_nodes() failed'


# collapse_nodes() for partial collapse
print('Testing collapse_nodes() for partial collapse...')
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/collapse.py"))
collapse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(collapse)
generate_updated_dictionary = collapse.generate_updated_dictionary
collapse_nodes = collapse.collapse_nodes

G_collapse_3 = nx.DiGraph()
G_collapse_3.add_node('N1')
G_collapse_3.nodes()['N1']['name'] = ['n_1']
G_collapse_3.nodes()['N1']['values'] = [1]
G_collapse_3.nodes()['N1']['stats'] = [0.1]
G_collapse_3.add_node('N2')
G_collapse_3.nodes()['N2']['name'] = ['n_2']
G_collapse_3.nodes()['N2']['values'] = [None]
G_collapse_3.nodes()['N2']['stats'] = [None]
G_collapse_3.add_node('N3')
G_collapse_3.nodes()['N3']['name'] = ['n_3']
G_collapse_3.nodes()['N3']['values'] = [3]
G_collapse_3.nodes()['N3']['stats'] = [0.3]
G_collapse_3.add_node('N4')
G_collapse_3.nodes()['N4']['name'] = ['n_4']
G_collapse_3.nodes()['N4']['values'] = [4]
G_collapse_3.nodes()['N4']['stats'] = [0.4]
G_collapse_3.add_node('N5')
G_collapse_3.nodes()['N5']['name'] = ['n_5']
G_collapse_3.nodes()['N5']['values'] = [None]
G_collapse_3.nodes()['N5']['stats'] = [None]
G_collapse_3.add_node('N6')
G_collapse_3.nodes()['N6']['name'] = ['n_6']
G_collapse_3.nodes()['N6']['values'] = [6]
G_collapse_3.nodes()['N6']['stats'] = [0.6]
G_collapse_3.add_node('N7')
G_collapse_3.nodes()['N7']['name'] = ['n_7']
G_collapse_3.nodes()['N7']['values'] = [None]
G_collapse_3.nodes()['N7']['stats'] = [None]
G_collapse_3.add_node('N8')
G_collapse_3.nodes()['N8']['name'] = ['n_8']
G_collapse_3.nodes()['N8']['values'] = [None]
G_collapse_3.nodes()['N8']['stats'] = [None]

G_collapse_3.add_node('N9')
G_collapse_3.nodes()['N9']['name'] = ['n_9']
G_collapse_3.nodes()['N9']['values'] = [1]
G_collapse_3.nodes()['N9']['stats'] = [0.1]
G_collapse_3.add_node('N10')
G_collapse_3.nodes()['N10']['name'] = ['n_10']
G_collapse_3.nodes()['N10']['values'] = [None]
G_collapse_3.nodes()['N10']['stats'] = [None]
G_collapse_3.add_node('N11')
G_collapse_3.nodes()['N11']['name'] = ['n_11']
G_collapse_3.nodes()['N11']['values'] = [None]
G_collapse_3.nodes()['N11']['stats'] = [None]
G_collapse_3.add_node('N12')
G_collapse_3.nodes()['N12']['name'] = ['n_12']
G_collapse_3.nodes()['N12']['values'] = [None]
G_collapse_3.nodes()['N12']['stats'] = [None]
G_collapse_3.add_node('N13')
G_collapse_3.nodes()['N13']['name'] = ['n_13']
G_collapse_3.nodes()['N13']['values'] = [None]
G_collapse_3.nodes()['N13']['stats'] = [None]
G_collapse_3.add_node('N14')
G_collapse_3.nodes()['N14']['name'] = ['n_14']
G_collapse_3.nodes()['N14']['values'] = [None]
G_collapse_3.nodes()['N14']['stats'] = [None]
G_collapse_3.add_node('N15')
G_collapse_3.nodes()['N15']['name'] = ['n_15']
G_collapse_3.nodes()['N15']['values'] = [None]
G_collapse_3.nodes()['N15']['stats'] = [None]
G_collapse_3.add_node('N16')
G_collapse_3.nodes()['N16']['name'] = ['n_16']
G_collapse_3.nodes()['N16']['values'] = [16]
G_collapse_3.nodes()['N16']['stats'] = [0.16]

G_collapse_3.add_node('N17')
G_collapse_3.nodes()['N17']['name'] = ['n_17']
G_collapse_3.nodes()['N17']['values'] = [None]
G_collapse_3.nodes()['N17']['stats'] = [None]

G_collapse_3.add_node('N18')
G_collapse_3.nodes()['N18']['name'] = ['n_18']
G_collapse_3.nodes()['N18']['values'] = [1]
G_collapse_3.nodes()['N18']['stats'] = [0.1]
G_collapse_3.add_node('N19')
G_collapse_3.nodes()['N19']['name'] = ['n_19']
G_collapse_3.nodes()['N19']['values'] = [2]
G_collapse_3.nodes()['N19']['stats'] = [0.2]
G_collapse_3.add_node('N20')
G_collapse_3.nodes()['N20']['name'] = ['n_20']
G_collapse_3.nodes()['N20']['values'] = [3]
G_collapse_3.nodes()['N20']['stats'] = [0.3]
G_collapse_3.add_node('N21')
G_collapse_3.nodes()['N21']['name'] = ['n_21']
G_collapse_3.nodes()['N21']['values'] = [None]
G_collapse_3.nodes()['N21']['stats'] = [None]
G_collapse_3.add_node('N22')
G_collapse_3.nodes()['N22']['name'] = ['n_22']
G_collapse_3.nodes()['N22']['values'] = [None]
G_collapse_3.nodes()['N22']['stats'] = [None]
G_collapse_3.add_node('N23')
G_collapse_3.nodes()['N23']['name'] = ['n_23']
G_collapse_3.nodes()['N23']['values'] = [None]
G_collapse_3.nodes()['N23']['stats'] = [None]
G_collapse_3.add_node('N24')
G_collapse_3.nodes()['N24']['name'] = ['n_24']
G_collapse_3.nodes()['N24']['values'] = [24]
G_collapse_3.nodes()['N24']['stats'] = [0.24]
G_collapse_3.add_node('N25')
G_collapse_3.nodes()['N25']['name'] = ['n_25']
G_collapse_3.nodes()['N25']['values'] = [None]
G_collapse_3.nodes()['N25']['stats'] = [None]
G_collapse_3.add_node('N26')
G_collapse_3.nodes()['N26']['name'] = ['n_26']
G_collapse_3.nodes()['N26']['values'] = [17]
G_collapse_3.nodes()['N26']['stats'] = [0.17]

collapser_database_3 = {
    'R1': {
        'compartment': 'none',
        'id': 'R1',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N1', 'N2'],
        'products': ['N3', 'N4', 'N5'],
        'modifiers': [],
        'additional_components': []
    },
    'R2': {
        'compartment': 'none',
        'id': 'R2',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N3', 'N4', 'N5'],
        'products': ['N6', 'N7', 'N8'],
        'modifiers': [],
        'additional_components': []
    },
    'R3': {
        'compartment': 'none',
        'id': 'R3',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N9', 'N10'],
        'products': ['N11', 'N12', 'N13'],
        'modifiers': [],
        'additional_components': []
    },
    'R4': {
        'compartment': 'none',
        'id': 'R4',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N13', 'N14'],
        'products': ['N15', 'N16'],
        'modifiers': [],
        'additional_components': []
    },
    'R5': {
        'compartment': 'none',
        'id': 'R5',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N13', 'N14'],
        'products': ['N15', 'N17'],
        'modifiers': [],
        'additional_components': []
    },
    'R6': {
        'compartment': 'none',
        'id': 'R6',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N18', 'N19'],
        'products': ['N20', 'N21', 'N22', 'N23'],
        'modifiers': [],
        'additional_components': []
    },
    'R7': {
        'compartment': 'none',
        'id': 'R7',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N23', 'N24'],
        'products': ['N25', 'N26'],
        'modifiers': [],
        'additional_components': []
    },
    'R8': {
        'compartment': 'none',
        'id': 'R8',
        'name': 'none',
        'reversible': 'none',
        'notes': 'none',
        'reactants': ['N1'],
        'products': ['N26'],
        'modifiers': [],
        'additional_components': []
    }
}

degree_dict = {
    'N1': 1,
    'N2': 1,
    'N3': 100,
    'N4': 200,
    'N5': 2,
    'N6': 1,
    'N7': 1,
    'N8': 1,
    'N9': 1,
    'N10': 1,
    'N11': 1,
    'N12': 1,
    'N13': 2,
    'N14': 1,
    'N15': 1,
    'N16': 1,
    'N17': 1,
    'N18': 1,
    'N19': 1,
    'N20': 2,
    'N21': 1,
    'N22': 2,
    'N23': 1,
    'N24': 1,
    'N25': 1,
    'N26': 1,
}

blocklist_3 = ['N20']

G_coll3 = G_collapse_3.copy()
G_coll3, updated_rxns3, changed_rxns3, removed_rxn3 = collapse_nodes(
    args_dict={},
    graph=G_coll3,
    reaction_dictionary=collapser_database_3,
    neighbors_dictionary={},
    degree_dictionary=degree_dict,
    samples=1,
    collapse_with_modifiers=False,
    blocklist=blocklist_3)

final_reactions3 = [
    'R1_R2',
    'R2_R1',
    'R3_R4',
    'R4_R3',
    'R5',
    'R6_R7',
    'R7_R6',
    'R8'
]

test_var = True
for k in list(updated_rxns3.keys()):
    if k in final_reactions3:
        pass
    else:
        test_var = False
assert test_var == True, 'collapse_nodes() for partial collapse failed'

# generate_updated_dictionary()
print("Testing generate_updated_dictionary()")
pathway_database = {
    'P1': {
        'reactome': 'Re1',
        'id': 'P1',
        'name': 'React1',
        'reactions': ['R1', 'R2', 'R3', 'R4']
    },
    'P2': {
        'reactome': 'Re2',
        'id': 'P2',
        'name': 'React2',
        'reactions': ['R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12', 'R13', ]
    },
    'P3': {
        'reactome': 'Re3',
        'id': 'P3',
        'name': 'React3',
        'reactions': ['R14', 'R15', 'R16', 'R17']
    }
}
updated_pathway_dictionary = generate_updated_dictionary(
    original_database=pathway_database,
    update_dictionary=changed_rxns1,
    removed_reaction=removed_rxn1)
assert updated_pathway_dictionary['Re1']['reactions'] == [
    'R1', 'R2_R3', 'R4'], 'generate_updated_dictionary() failed'
assert updated_pathway_dictionary['Re2']['reactions'] == [
    'R5_R6_R7', 'R8_R9_R10', 'R11', 'R12', 'R13'], 'generate_updated_dictionary() failed'
assert updated_pathway_dictionary['Re3']['reactions'] == [
    'R14', 'R15', 'R16', 'R17'], 'generate_updated_dictionary() failed'

os.remove(network_url)

print('Tests completed')
