"""License Information
MetaboNet:

    Thomas Cameron Waller
    tcameronwaller@gmail.com
    Department of Biochemistry
    University of Utah
    Room 4100, Emma Eccles Jones Medical Research Building
    15 North Medical Drive East
    Salt Lake City, Utah 84112
    United States of America

    Portions of this code are modified from MetaboNet
    (https://github.com/tcameronwaller/metabonet/).

    MetaboNet supports definition and analysis of custom metabolic networks.
    Copyright (C) 2019 Thomas Cameron Waller

    MetaboNet is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    MetaboNet is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with MetaboNet. If not, see <http://www.gnu.org/licenses/>.

MetaboNet-Analyzer:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metabalyze
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

"""Import internal dependencies
"""
from metabalyze.utils import progress_bar
from metabalyze.metabonet_network.utils import remove_file, remove_directory

"""Clean intermediate data files
"""
def clean_collect(
        args_dict):

    # Collection.
    remove_file(args_dict['collect'] + 'compartments.pickle')
    remove_file(args_dict['collect'] + 'processes.pickle')
    remove_file(args_dict['collect'] + 'reactions.pickle')
    remove_file(args_dict['collect'] + 'metabolites.pickle')
    remove_file(args_dict['collect'] + 'reactions.tsv')
    remove_file(args_dict['collect'] + 'metabolites.tsv')
    remove_directory(args_dict['collect'])

def clean_extract(
        args_dict):

    # Extraction.
    remove_file(args_dict['extract'] + 'hmdb_summary.pickle')
    remove_file(args_dict['extract'] + 'hmdb_summary.tsv')
    remove_directory(args_dict['extract'])

def clean_enhance(
        args_dict):

    # Enhancement.
    remove_file(args_dict['enhance'] + 'compartments.pickle'))
    remove_file(args_dict['enhance'] + 'processes.pickle'))
    remove_file(args_dict['enhance'] + 'reactions.pickle'))
    remove_file(args_dict['enhance'] + 'metabolites.pickle'))
    remove_file(args_dict['enhance'] + 'reactions.tsv'))
    remove_file(args_dict['enhance'] + 'metabolites.tsv'))
    remove_file(args_dict['enhance'] + 'reactions_filter.tsv'))
    remove_directory(args_dict['enhance'])

def clean_curate(
        args_dict):

    # Curation.
    remove_file(args_dict['curate'] + 'compartments.pickle')
    remove_file(args_dict['curate'] + 'processes.pickle')
    remove_file(args_dict['curate'] + 'reactions.pickle')
    remove_file(args_dict['curate'] + 'metabolites.pickle')
    remove_file(args_dict['curate'] + 'reactions.tsv')
    remove_file(args_dict['curate'] + 'metabolites.tsv')
    remove_directory(args_dict['curate'])

def clean_convert(
        args_dict):

    # Conversion.
    remove_file(args_dict['convert'] + 'compartments.pickle')
    remove_file(args_dict['convert'] + 'processes.pickle')
    remove_file(args_dict['convert'] + 'reactions.pickle')
    remove_file(args_dict['convert'] + 'metabolites.pickle')
    remove_file(args_dict['convert'] + 'compartments.tsv')
    remove_file(args_dict['convert'] + 'processes.tsv')
    remove_file(args_dict['convert'] + 'reactions.tsv')
    remove_file(args_dict['convert'] + 'metabolites.tsv')
    remove_file(args_dict['convert'] + 'dymetabonet.json')
    remove_directory(args_dict['convert'])
