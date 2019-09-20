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
from metaboverse.metaboverse_curate.fetch_species import fetch_species
from metaboverse.metaboverse_curate.fetch_version import check_reactome_version
from metaboverse.metaboverse_curate.utils import get_table
from metaboverse.metaboverse_curate.utils import unpack_table
from metaboverse.metaboverse_curate.__main__ import parse_table
from metaboverse.metaboverse_curate.__main__ import parse_complexes
from metaboverse.metaboverse_curate.__main__ import make_master
from metaboverse.metaboverse_curate.__main__ import prepare_mappers
from metaboverse.metaboverse_curate.__main__ import map_complexes_genes
from metaboverse.metaboverse_curate.__main__ import write_database
from metaboverse.metaboverse_curate.load_reactions_db import * # add others <---------------------------

"""Set globals
"""
__path__ = os.path.dirname(os.path.realpath(__file__)) + '/'
__path__ = '/Users/jordan/Desktop/Metaboverse/tests/curation_tests/'

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
        ref_dict):

    complexes = parse_complexes(
        reference=ref_dict)

    # <--- continue writing this test


"""Run tests
"""
def __main__():

    test_species_list()
    test_reactome_version()
    test_util_get_table()
    test_util_unpack_table()

    data_1 = [
        ['analyte_id', 'analyte_name', 'reaction_id', 'reaction_name', 'source_id', 'extra_column'],
        ['R-ALL-00000', 'fake_molecule [organelle]', 'R-HSA-00000', 'fake_reaction', '0101', '????']]
    ref_dict = test_parse_table(data_1)

    data_2 = [
        ['analyte_id', 'analyte_name', 'reaction_id', 'reaction_name'],
        ['R-ALL-00000', 'fake_molecule [organelle]', 'R-HSA-00000', 'fake_reaction']]
    ref_dict = test_parse_table(data_2)

    complex_ref = {

    }
    test_parse_complexes(complex_ref)
