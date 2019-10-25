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
import pandas as pd

"""Import internal dependencies
"""
from metaboverse.curate.fetch_species import fetch_species
from metaboverse.curate.fetch_version import check_reactome_version

from metaboverse.curate.utils import get_table
from metaboverse.curate.utils import unpack_table

from metaboverse.curate.__main__ import parse_table
from metaboverse.curate.__main__ import parse_complexes
from metaboverse.curate.__main__ import make_master
from metaboverse.curate.__main__ import prepare_mappers
from metaboverse.curate.__main__ import map_complexes_genes
from metaboverse.curate.__main__ import write_database

from metaboverse.curate.load_reactions_db import get_reactions
from metaboverse.curate.load_reactions_db import unpack_reactions
from metaboverse.curate.load_reactions_db import get_database
from metaboverse.curate.load_reactions_db import curate_reactions
from metaboverse.curate.load_reactions_db import add_pathways
from metaboverse.curate.load_reactions_db import add_compartments

"""Set globals
"""
__path__ = os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/Desktop/Metaboverse/tests/curation_tests/'

"""Test that staple species were grabbed from species list
"""
def test_species_list():

    staple_species = [
        'Danio rerio',
        'Drosophila melanogaster',
        'Escherichia coli',
        'Homo sapiens',
        'Mus musculus',
        'Rattus norvegicus',
        'Saccharomyces cerevisiae',
        'Schizosaccharomyces pombe',
        'Xenopus laevis',
        'Xenopus tropicalis']

    organisms = fetch_species()

    for x in staple_species:

        if x not in organisms:
            raise Exception('Could not find', x, 'in Reactome species list')

"""Test Reactome version checking
"""
def test_reactome_version():

    check_reactome_version()

"""Test ability to download Reactome file
"""
def test_util_unpack_table():

    url = 'https://reactome.org/download/current/'
    test_file = 'models2pathways.tsv'

    file = unpack_table(
            url + test_file,
            output_dir=__path__)

    assert file == __path__ + test_file, 'Problem getting table download'
    os.remove(file)

"""Test ability to download file and unpack table
"""
def test_util_get_table():

    url = 'https://reactome.org/download/current/'
    test_file = 'models2pathways.tsv'

    table = get_table(
        __path__,
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
        organism_key='7') # Key for testing purposes only

    assert 'Homo sapiens' in table['7'].tolist(), 'Problem getting table download'
    os.remove(__path__ + test_file)

"""
"""
def test_parse_table(
        data):

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

    table = pd.DataFrame(data_1, columns=data_1[0])
    table = table.drop(0, axis=0)

    reference = {
        'this_one': table,
        'not_this_one': False}

    ref_dict = parse_table(
            reference=reference,
            key='this_one')

    assert run_checks(ref_dict) == True, 'Problem parsing Reactome table'

    return ref_dict

def test_parse_complexes(
        complex_ref):

    truth = {
        'R-ALL-00000': {
            'complex_id': 'R-ALL-00000',
            'complex_name': 'fake_complex',
            'compartment': 'organelle',
            'participating_complex': 'R-HSA-00000',
            'pathway': 'R-PATH-00000',
            'top_level_pathway': 'R-TOP-0000',
            'participants': {
                'chebi': ['00000'],
                'uniprot': ['0001'],
                'ensembl': [],
                'mirbase': ['00000'],
                'ncbi': []}}}

    complexes = parse_complexes(
        reference=complex_ref)

    assert complexes == truth, 'Unable to parse complex data'

def test_make_master():

    test = {
        'chebi_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'analyte': 'example',
                'other': None
            }
        },
        'uniprot_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'analyte': 'example',
                'other': None
            }
        },
        'ensembl_reference': {
            'R-HSA-0001': {
                'analyte_id': 'R-ANA-0001',
                'analyte': 'example1',
                'other': None
            }
        },
        'ncbi_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'analyte': 'example',
                'other': None
            }
        },
        'mirbase_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'analyte': 'example',
                'other': None
            }
        },
        'complexes_reference': {
            'R-COM-0000': {
                'complex_id': 'R-COMP-0000',
                'complex_name': 'complex',
                'other': None
            }
        }
    }

    truth = {
        'R-ANA-0000': 'example',
        'R-ANA-0001': 'example1',
        'R-COMP-0000': 'complex'}

    test_output = make_master(test)

    assert test_output == truth, 'Unable to generate master reference'

def test_prepare_mappers():

    test = {
        'chebi_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'source_id': 'example',
                'other': None
            }
        },
        'uniprot_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'source_id': 'example',
                'other': None
            }
        },
        'ensembl_reference': {
            'R-HSA-0001': {
                'analyte_id': 'R-ANA-0001',
                'source_id': 'example1',
                'other': None
            }
        },
        'ncbi_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'source_id': 'example',
                'other': None
            }
        },
        'mirbase_reference': {
            'R-HSA-0000': {
                'analyte_id': 'R-ANA-0000',
                'source_id': 'example',
                'other': None
            }
        }
    }

    # Not expecting anything from complexes or ncbi
    truth = {
        'ensembl': {
            'example1': 'R-ANA-0001'},
        'mirbase': {
            'example': 'R-ANA-0000'},
        'uniprot': {
            'example': 'R-ANA-0000'},
        'chebi': {
            'example': 'R-ANA-0000'}}

    test_output = prepare_mappers(test)

    assert test_output == truth, 'Unable to generate master mapper reference'

def test_map_complexes_genes():

    test = {
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

    complex = [
        ['identifier', 'name', 'participants', 'participatingComplex'],
        ['R-ALL-00000', 'fake_complex [organelle]', 'chebi:0017|uniprot:0015|mirbase:0012', 'R-HSA-00000']]
    complex = pd.DataFrame(
        complex,
        columns=complex[0])
    complex = complex.drop(
        0,
        axis=0)

    truth = {
        'R-ALL-00000': {
            'id': 'R-ALL-00000',
            'participants': [
                'R-ANA-1017',
                'R-ANA-1015',
                'R-ANA-1012'],
            'partner_complexes': [
                'R-HSA-00000']}}

    mapped = map_complexes_genes(
        complex_participants=complex,
        databases=test)

    assert mapped == truth, 'Unable to generate complex to gene mapper'

def test_write_database():

    output = __path__
    file = 'test.pickle'
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

    write_database(
        output=output,
        file=file,
        database=database)

    assert os.path.exists(output + file) == True, 'Failed to write database to pickle file'
    os.remove(output + file)

def test_get_reactions():

    reactions = get_reactions(
        species_id='HSA',
        reaction_dir=__path__)

    reactions

    assert reactions == ['R-HSA-test'], 'Unable to fetch reaction file list for organism'

"""Run tests
"""
def __main__():

    # Utils
    test_species_list()
    test_reactome_version()
    test_util_get_table()
    test_util_unpack_table()

    # Parse table
    data_1 = [
        ['analyte_id', 'analyte_name', 'reaction_id', 'reaction_name', 'source_id', 'extra_column'],
        ['R-ALL-00000', 'fake_molecule [organelle]', 'R-HSA-00000', 'fake_reaction', '0101', '????']]
    ref_dict = test_parse_table(data_1)

    data_2 = [
        ['analyte_id', 'analyte_name', 'reaction_id', 'reaction_name'],
        ['R-ALL-00000', 'fake_molecule [organelle]', 'R-HSA-00000', 'fake_reaction']]
    ref_dict = test_parse_table(data_2)

    complex_pathway = [
        ['R-ALL-00000', 'R-PATH-00000', 'R-TOP-0000', 'something else']]
    complex_pathway = pd.DataFrame(complex_pathway)
    complex = [
        ['identifier', 'name', 'participants', 'participatingComplex'],
        ['R-ALL-00000', 'fake_complex [organelle]', 'chebi:00000|uniprot:0001|mirbase:00000', 'R-HSA-00000']]
    complex = pd.DataFrame(
        complex,
        columns=complex[0])
    complex = complex.drop(
        0,
        axis=0)
    complex_ref = {
        'complex_pathway': complex_pathway,
        'complex_participants': complex}

    # Parse complexes
    test_parse_complexes(complex_ref)

    # Other database modifications
    test_make_master()
    test_prepare_mappers()
    test_map_complexes_genes()
    test_write_database()

    # Reaction loading
    test_get_reactions()

    truth_reference = {}

    file = 'R-HSA-realtest'
    sbml = get_database(
        reaction_dir=__path__,
        reaction_name=file)
    tag = '{http://www.sbml.org/sbml/level3/version1/core}sbml'
    assert sbml.tag == tag, 'Failed to import SBML file'

    reference = {}
    reference['pathways'] = curate_reactions(
        reaction_dir=__path__,
        reactions_list=[file],
        species_id='HSA')
    reference['pathway_types'] = add_pathways(
        reference['pathways'])
    reference['compartment_types'] = add_compartments(
        reference['pathways'])

    assert list(reference['pathway_types'].keys()) == ['Virus Assembly and Release'], 'Failed to gather pathway name'

    compartments = {
        'R-HSA-12045',
        'R-HSA-17957',
        'R-HSA-20699',
        'R-HSA-24337',
        'R-HSA-70101',
        'R-HSA-876',
        'R-HSA-984'}
    assert compartments == reference['compartment_types'], 'Unable to extract correct compartments'

    assert reference['pathways']['R-HSA-realtest']['reactome_id'] == 'R-HSA-realtest', 'Unable to extract reactome id'
    assert reference['pathways']['R-HSA-realtest']['pathway_name'] == 'Virus Assembly and Release', 'Unable to extract pathway name'
    assert reference['pathways']['R-HSA-realtest']['reactions']['R-HSA-168875']['id'] == 'R-HSA-168875'
    assert reference['pathways']['R-HSA-realtest']['reactions']['R-HSA-168858']['name'] == 'Palmitoylation of cysteine residues on HA in the cis-Golgi network', 'Unable to extract reaction name'
    assert reference['pathways']['R-HSA-realtest']['reactions']['R-HSA-168858']['reactants']['R-ALL-195819']['species_id'] == 'R-ALL-195819', 'Unable to extract sub-reaction level information'
