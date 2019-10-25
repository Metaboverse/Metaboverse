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

"""Import internal dependencies
"""
from metaboverse.curate.__main__ import __main__ as curate

"""Set globals
"""
__path__ = os.path.dirname(os.path.realpath(__file__)) + '/'
#__path__ = '/Users/jordan/Desktop/Metaboverse/tests/curation_tests/'

args_dict = {
    'species': 'SCE',
    'output': __path__}

curate(
    args_dict=args_dict)

with open(__path__ + 'SCE_metaboverse_db.pickle', 'rb') as network_file:
    reactome_database = pickle.load(network_file)

assert reactome_database['master_reference']['R-ALL-389536'] == 'CO2', 'Unable to extract element from master reference'
assert reactome_database['pathways']['R-HSA-2562578']['reactions']['R-HSA-2562541']['name'] == 'TLR4-induced ripoptosome assembly', 'Unable to extract reaction name'

os.remove(__path__ + 'SCE_metaboverse_db.pickle')
