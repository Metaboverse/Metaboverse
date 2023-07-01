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
import importlib.util
import xml.etree.ElementTree as et
import pandas as pd
from shutil import copyfile
import pickle
import os

"""Curation/Utils
"""
spec = importlib.util.spec_from_file_location(
    "__main__", os.path.abspath("./metaboverse_cli/curate/__main__.py"))
curate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(curate)

spec = importlib.util.spec_from_file_location(
    "__main__", os.path.abspath("./metaboverse_cli/curate/load_reactions_db.py"))
load_reactions = importlib.util.module_from_spec(spec)
spec.loader.exec_module(load_reactions)

spec = importlib.util.spec_from_file_location(
    "__main__", os.path.abspath("./metaboverse_cli/curate/load_complexes_db.py"))
load_complexes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(load_complexes)

spec = importlib.util.spec_from_file_location(
    "__main__", os.path.abspath("./metaboverse_cli/curate/fetch_species.py"))
fetch_species = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fetch_species)

spec = importlib.util.spec_from_file_location(
    "get_table", os.path.abspath("./metaboverse_cli/curate/utils.py"))
get_table = importlib.util.module_from_spec(spec)
spec.loader.exec_module(get_table)

spec = importlib.util.spec_from_file_location(
    "unpack_table", os.path.abspath("./metaboverse_cli/curate/utils.py"))
unpack_table = importlib.util.module_from_spec(spec)
spec.loader.exec_module(unpack_table)

spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/curate/load_reactions_db.py"))
load_reactions_db = importlib.util.module_from_spec(spec)
spec.loader.exec_module(load_reactions_db)

# test __main__() -- functional test
args_dict = {
    'organism_id': 'SCE',
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "curate", "test")) + os.path.sep}

# load_reactions_db.py
args_dict, pathway_database, reaction_database, species_database, \
name_database, compartment_dictionary, \
components_database = load_reactions.__main__(
        species_id=args_dict['organism_id'],
        output_dir=args_dict['output'],
        database_source='reactome',
        sbml_url='None',
        args_dict=args_dict)

assert type(pathway_database) == dict, "load_reactions_db.py failed"
assert type(pathway_database[list(pathway_database.keys())[
            0]]) == dict, "load_reactions_db.py failed"

assert type(reaction_database) == dict, "load_reactions_db.py failed"
assert type(reaction_database[list(reaction_database.keys())[
            0]]) == dict, "load_reactions_db.py failed"

assert type(species_database) == dict, "load_reactions_db.py failed"
assert type(species_database[list(species_database.keys())[
            0]]) == str, "load_reactions_db.py failed"

assert type(name_database) == dict, "load_reactions_db.py failed"
assert type(name_database[list(name_database.keys())[0]]
            ) == str, "load_reactions_db.py failed"

assert type(compartment_dictionary) == dict, "load_reactions_db.py failed"
assert type(compartment_dictionary[list(compartment_dictionary.keys())[
            0]]) == str, "load_reactions_db.py failed"

assert type(components_database) == dict, "load_reactions_db.py failed"
assert type(components_database[list(components_database.keys())[
            0]]) == dict, "load_reactions_db.py failed"

# load_complexes_db.py
d = load_complexes.__main__(
    output_dir=args_dict['output'])
assert type(d['complex_participants']
            ) == pd.DataFrame, "load_complexes_db.py failed"
assert type(d['complex_pathway']
            ) == pd.DataFrame, "load_complexes_db.py failed"

# fetch_species.py
organisms = fetch_species.__main__()
assert type(organisms) == list, "fetch_species.py failed"

# write_database()
args_dict = {
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "curate", "test")) + os.path.sep,
    'file': 'test.mvdb'}
database = {
    'chebi_reference': {
        'R-HSA-0000': {
            'analyte_id': 'R-ANA-1017',
            'source_id': '0017',
            'other': None
        }
    },
    'uniprot_reference': {
        'R-HSA-0000': {
            'analyte_id': 'R-ANA-1015',
            'source_id': '0015',
            'other': None
        }
    },
    'ensembl_reference': {
        'R-HSA-0001': {
            'analyte_id': 'R-ANA-1014',
            'source_id': '0014',
            'other': None
        }
    },
    'ncbi_reference': {
        'R-HSA-0000': {
            'analyte_id': 'R-ANA-1013',
            'source_id': '0013',
            'other': None
        }
    },
    'mirbase_reference': {
        'R-HSA-0000': {
            'analyte_id': 'R-ANA-1012',
            'source_id': '0012',
            'other': None
        }
    }
}
curate.write_database(
    output=args_dict['output'],
    file=args_dict['file'],
    database=database)

assert os.path.exists(args_dict['output'] + args_dict['file']
                      ) == True, 'Failed to write database to pickle file'
os.remove(args_dict['output'] + args_dict['file'])

# unpack_table()
url = 'https://reactome.org/download/current/'
test_file = 'models2pathways.tsv'
file = unpack_table.unpack_table(
    url + test_file,
    output_dir=args_dict['output'])
assert file == args_dict['output'] + test_file, 'unpack_table() failed'
os.remove(args_dict['output'] + test_file)

# get_table()
url = 'https://reactome.org/download/current/'
test_file = 'models2pathways.tsv'
table = get_table.get_table(
    args_dict['output'],
    url + test_file,
    column_names=[
        '1',
        '2',
        '3',
        '4',
        '5',
        '6',
        '7'],
    organism='Homo sapiens',
    organism_key='7')  # Key for testing purposes only
assert 'Homo sapiens' in table['7'].tolist(), 'Problem getting table download'
os.remove(args_dict['output'] + test_file)

# parse_table()


def run_checks(
        ref_dict):

    if ref_dict['R-ALL-00000']['analyte_id'] != 'R-ALL-00000' \
            or ref_dict['R-ALL-00000']['reaction_id'] != 'R-HSA-00000' \
            or ref_dict['R-ALL-00000']['reaction_name'] != 'fake_reaction' \
            or ref_dict['R-ALL-00000']['analyte'] != 'fake_molecule' \
            or ref_dict['R-ALL-00000']['compartment'] != 'organelle':
        return False

    else:
        return True


data_1 = [
    ['analyte_id', 'analyte_name', 'reaction_id',
        'reaction_name', 'source_id', 'extra_column'],
    ['R-ALL-00000', 'fake_molecule [organelle]',
        'R-HSA-00000', 'fake_reaction', '0101', '????']
]
table = pd.DataFrame(data_1, columns=data_1[0])
table = table.drop(0, axis=0)
reference = {
    'this_one': table,
    'not_this_one': False}
ref_dict = curate.parse_table(
    reference=reference,
    key='this_one')
assert run_checks(ref_dict) == True, 'Problem parsing Reactome table'

# Test unpacking of a reaction file
species_id = 'HSA'
path = args_dict['output']
file = 'R-HSA-realtest.sbml'
test_db = args_dict['output'] + file
sbml_namespace = '{{http://www.sbml.org/sbml/level{0}/version{1}/core}}'
rdf_namespace = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
bqbiol_namespace = '{http://biomodels.net/biology-qualifiers/}'
sbml_level = '3'
sbml_version = '1'

path_list = load_reactions_db.get_pathways(
    species_id=species_id,
    pathways_dir=path)
assert path_list == ['R-HSA-realtest'], 'get_pathways() failed'

contents = load_reactions_db.get_database(
    pathways_dir=path,
    pathway_name=path_list[0])
assert type(contents) == et.Element, 'get_database() failed'

args_dict, pathway_database, reaction_database, species_database, \
name_database, compartment_database, compartment_dictionary, \
components_database = load_reactions_db.process_components(
        output_dir=args_dict['output'],
        pathways_dir=args_dict['output'],
        pathways_list=path_list,
        species_id=species_id,
        args_dict=args_dict,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace)

assert 'R-HSA-realtest' in pathway_database.keys(), 'process_components() failed'
assert pathway_database['R-HSA-realtest']['id'] == 'pathway_168268', 'process_components() failed'
assert pathway_database['R-HSA-realtest']['reactome'] == 'R-HSA-realtest', 'process_components() failed'
assert pathway_database['R-HSA-realtest']['name'] == 'Virus Assembly and Release', 'process_components() failed'
assert len(pathway_database['R-HSA-realtest']
           ['reactions']) == 21, 'process_components() failed'
truth_d = set([
    'reaction_168884',
    'reaction_168875',
    'reaction_168882',
    'reaction_195733',
    'reaction_169919',
    'reaction_195739',
    'reaction_168894',
    'reaction_168870',
    'reaction_168871',
    'reaction_195734',
    'reaction_168858',
    'reaction_168862',
    'reaction_168895',
    'reaction_168860',
    'reaction_195726',
    'reaction_169921',
    'reaction_195926',
    'reaction_195730',
    'reaction_168869',
    'reaction_188544',
    'reaction_169847'])
assert set(pathway_database['R-HSA-realtest']['reactions']
           ) == truth_d, 'process_components() failed'

assert reaction_database['reaction_195733']['compartment'] == 'compartment_17957', 'process_components() failed'
assert reaction_database['reaction_195733']['id'] == 'reaction_195733', 'process_components() failed'
assert reaction_database['reaction_195733']['name'] == 'M2 protein synthesis', 'process_components() failed'
assert reaction_database['reaction_195733']['reversible'] == 'false', 'process_components() failed'
assert type(reaction_database['reaction_195733']
            ['notes']) == str, 'process_components() failed'
assert reaction_database['reaction_195733']['reactants'] == [
], 'process_components() failed'
assert reaction_database['reaction_195733']['products'] == [
    'species_195756'], 'process_components() failed'
assert reaction_database['reaction_195733']['modifiers'] == [
], 'process_components() failed'
assert reaction_database['reaction_168870']['reactants'] == [
    'species_195945'], 'process_components() failed'
assert reaction_database['reaction_168870']['products'] == [
    'species_195941'], 'process_components() failed'
assert reaction_database['reaction_168870']['modifiers'] == [
    ['species_195915', 'catalyst']], 'process_components() failed'

assert species_database['species_195752'] == 'NP', 'process_components() failed'
assert species_database['species_195915'] == 'NA', 'process_components() failed'

assert name_database['NP'] == 'species_195752', 'process_components() failed'
assert name_database['NA'] == 'species_195915', 'process_components() failed'

assert compartment_database['species_195752'] == 'compartment_876', 'process_components() failed'
assert compartment_database['species_195915'] == 'compartment_984', 'process_components() failed'

assert compartment_dictionary['compartment_876'] == 'plasma membrane', 'process_components() failed'
assert compartment_dictionary['compartment_984'] == 'extracellular region', 'process_components() failed'

assert components_database['species_195752']['id'] == 'species_195752', 'process_components() failed'
assert components_database['species_195752']['reactome_id'] == 'R-FLU-195752', 'process_components() failed'
assert components_database['species_195752']['name'] == 'NP', 'process_components() failed'
assert components_database['species_195752']['is'] == 'P03466', 'process_components() failed'
assert components_database['species_195752']['hasPart'] == [
], 'process_components() failed'
print(components_database['species_195752'])
assert components_database['species_195752']['type'] == 'protein_component', 'process_components() failed'
assert components_database['species_195752']['compartment'] == 'compartment_876', 'process_components() failed'

print(components_database['species_195941'])

assert set([
    'P03466',
    'P03428',
    'P03433',
    'P03431',
    'P03485',
    'P03508',
    'P03468',
    'P03452',
    'P06821']).issubset(set(components_database['species_195941']['hasPart'])), 'process_components() failed'

# BioModels tests
src = os.path.abspath(os.path.join(".", "metaboverse_cli",
                                   "curate", "test", "test_session_data.json"))
dst = os.path.abspath(os.path.join(".", "metaboverse_cli",
                                   "curate", "test", "test_session_data_copy.json"))
s_out = copyfile(src, dst)

args_dict = {
    'organism_id': 'find',
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "curate", "test")) + os.path.sep,
    'session_data': os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data_copy.json"))}

args_dict, pathway_database, reaction_database, species_database, \
name_database, compartment_dictionary, \
components_database = load_reactions.__main__(
    species_id=args_dict['organism_id'],
    output_dir=args_dict['output'],
    database_source='biomodels/bigg',
    sbml_url=os.path.abspath(os.path.join(
        ".", "metaboverse_cli", "curate", "test", "iIS312.xml")),
    args_dict=args_dict)

if 'All' in pathway_database:
    print("Pass")
else:
    raise Exception("BioModels/BiGG database curation failed")
assert reaction_database['R_PROt2r_copy2']['reactants'] == [
    'M_h_e', 'M_pro__L_e'], "BioModels/BiGG database curation failed"
assert species_database['M_pro__L_e'] == 'L-Proline', "BioModels/BiGG database curation failed"
assert name_database['M_pro__L_e'] == 'M_pro__L_e', "BioModels/BiGG database curation failed"
assert name_database['L-Proline'] == 'M_pro__L_e', "BioModels/BiGG database curation failed"
assert compartment_dictionary['c'] == 'cytosol', "BioModels/BiGG database curation failed"
assert components_database['M_pro__L_e']['is'] == 'D00035', "BioModels/BiGG database curation failed"

os.remove(dst)


# BiGG tests
src = os.path.abspath(os.path.join(".", "metaboverse_cli",
                                   "curate", "test", "test_session_data.json"))
dst = os.path.abspath(os.path.join(".", "metaboverse_cli",
                                   "curate", "test", "test_session_data_copy.json"))
s_out = copyfile(src, dst)

args_dict = {
    'organism_id': 'find',
    'database_source': 'biomodels/bigg',
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "curate", "test")) + os.path.sep,
    'session_data': os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data_copy.json"))}

args_dict, pathway_database, reaction_database, species_database, \
name_database, compartment_dictionary, \
components_database = load_reactions.__main__(
    species_id=args_dict['organism_id'],
    output_dir=args_dict['output'],
    database_source='biomodels/bigg',
    sbml_url=os.path.abspath(os.path.join(
        ".", "metaboverse_cli", "curate", "test", "BMID000000141967_url.xml")),
    args_dict=args_dict)

if 'All' in pathway_database:
    print("Pass")
else:
    raise Exception("BioModels/BiGG database curation failed")
assert reaction_database['MNXR42496_i']['reactants'] == [
    'MNXM2885_i', 'bigg_gthrd_i'], "BioModels/BiGG database curation failed"
assert species_database['MNXM2885_i'] == '(1S,2R)-Naphthalene 1,2-oxide', "BioModels/BiGG database curation failed"
assert species_database['MNXM2885_e'] == '(1S,2R)-Naphthalene 1,2-oxide', "BioModels/BiGG database curation failed"
assert name_database['(1S,2R)-Naphthalene 1,2-oxide'] == 'MNXM2885_e', "BioModels/BiGG database curation failed"
assert compartment_dictionary['i'] == 'intracellular', "BioModels/BiGG database curation failed"
assert components_database['MNXM2885_e']['is'] == 'MNXM2885_e', "BioModels/BiGG database curation failed"

os.remove(dst)

# test all
"""
src = os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data.json"))
dst = os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data_copy.json"))
s_out = copyfile(src, dst)

args_dict = {
    'organism_id': 'SCE',
    'database_source': 'reactome',
    'organism_curation': "None",
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "curate", "test")) + os.path.sep,
    'session_data':os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data_copy.json"))}
curate.__main__(
    args_dict=args_dict
)
os.remove(dst)
os.remove(args_dict['output'] + "SCE.mvdb")
"""

# test all bigg
#src = os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data.json"))
#dst = os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data_copy.json"))
#s_out = copyfile(src, dst)
# args_dict = {
#    'organism_id': 'find',
#    'database_source': 'biomodels/bigg',
#    'organism_curation': os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "BMID000000141967_url.xml")),
#    'output': os.path.abspath(
#        os.path.join(".", "metaboverse_cli", "curate", "test")) + os.path.sep,
#    'session_data':os.path.abspath(os.path.join(".", "metaboverse_cli", "curate", "test", "test_session_data_copy.json"))}
# curate.__main__(
#    args_dict=args_dict
# )
# os.remove(dst)
#os.remove(args_dict['output'] + "BMID000000141967.mvdb")

print('Tests completed')
