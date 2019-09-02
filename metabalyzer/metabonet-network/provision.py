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
    alias: metabalyzer
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
import shutil
import csv
import copy
import pickle
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from metabalyzer.metabonet-network.utils import confirm_path_directory
from metabalyzer.metabonet-network.utils import write_file_table
from metabalyzer.metabonet-network.utils import read_file_table
from metabalyzer.metabonet-network.enhance import filter_hmdb_entries_by_identifiers


"""Reads and organizes source information from file
arguments:
    directory (str): directory of source files
returns:
    (object): source information
"""
def read_source(
        directory):

    # Specify directories and files.
    path_hmdb = os.path.join(directory, 'hmdb_metabolites.xml')
    path = os.path.join(directory, 'provision')
    path_midas = os.path.join(path, 'midas_original.tsv')

    # Read information from file.
    #hmdb = et.parse(path_hmdb)
    midas_original = read_file_table(
        path_file=path_midas,
        names=None,
        delimiter="\t")

    # Compile and return information.
    return {
        "hmdb": path_hmdb,
        "midas_original": midas_original}

"""Transfers information from Human Metabolome Database (HMDB) to MIDAS
library.
arguments:
    summary_hmdb (dict<dict>): information from HMDB
    midas_original (list<dict<str>>): information in MIDAS
returns:
    (list<dict<str>>): information in MIDAS
"""
def transfer_summary_midas(
        summary_hmdb=None,
        midas_original=None):

    midas_novel = []

    for record_midas in midas_original:

        # Interpretation.
        reference_midas = record_midas['reference_midas']
        name_original = record_midas['name']
        reference_hmdb_original = record_midas['reference_hmdb']
        reference_kegg_original = record_midas['reference_kegg']

        # Determine whether MIDAS record matches a HMDB record.
        if len(reference_hmdb_original) > 0:

            # Match MIDAS record to HMDB record.
            keys_hmdb = filter_hmdb_entries_by_identifiers(
                identifiers=[reference_hmdb_original],
                metabolites_references=summary_hmdb)

            if len(keys_hmdb) > 1:

                print('found multiple hmdb keys!')

            key_hmdb = keys_hmdb[0]
            record_hmdb = summary_hmdb[key_hmdb]
            record_novel = copy.deepcopy(record_hmdb)
            record_novel['reference_midas'] = reference_midas
            record_novel['name_original'] = name_original
            record_novel['reference_kegg_original'] = reference_kegg_original
            midas_novel.append(record_novel)

        else:

            keys_hmdb = list(summary_hmdb.values())[0].keys()
            record_hmdb = {}

            for key in keys_hmdb:

                record_hmdb[key] = 'null'

            record_novel = record_hmdb
            record_novel['reference_midas'] = identifier_midas
            record_novel['name_original'] = name_original
            record_novel['reference_kegg_original'] = reference_kegg_original
            midas_novel.append(record_novel)

    return midas_novel

"""Determine neighbors
Not currently in use
"""
def determine_reactant_product_neighbors():

    # Possibly even iterate only on the reactions already known to involve the
    # metabolite.
    for reaction in reactions.values():
        # Get reaction's reactants and products.
        #reactants =
        #products =
        # Determine whether query metabolite matches any reactants of products.
        # If query metabolite matches, then collect all reactants and products.
        pass

"""Writes product information to file
arguments:
    directory (str): directory for product files
    information (object): information to write to file
"""
def write_product(
        directory,
        information=None):

    # Specify directories and files.
    confirm_path_directory(directory)
    path_pickle = os.path.join(directory, 'hmdb_summary.pickle')
    path_text = os.path.join(directory, 'hmdb_summary.tsv')
    path_midas = os.path.join(directory, 'midas_novel.tsv')

    # Write information to file.
    with open(path_pickle, 'wb') as file_product:
        pickle.dump(information['summary_object'], file_product)

    write_file_table(
        information=information['summary_list'],
        path_file=path_text,
        names=information['summary_list'][0].keys(),
        delimiter='\t')
    write_file_table(
        information=information['midas_novel'],
        path_file=path_midas,
        names=information['midas_novel'][0].keys(),
        delimiter='\t')

"""Function to execute module's main behavior.
The purpose of this procedure is to extract relevant information from the
Human Metabolome Database.
arguments:
    directory (str): path to directory for source and product files
"""
def execute_procedure(
        args_dict):

    # TODO: read in MIDAS, hmdb summary, and information from

    # Read source information from file.
    source = read_source(
        args_dict['source'])

    # Transfer information to MIDAS library.
    midas_novel = transfer_summary_midas(
        summary_hmdb=summary_hmdb,
        midas_original=source["midas_original"])

    # Compile information.
    information = {
        "midas_novel": midas_novel}

    #Write product information to file
    write_product(
    args_dict['provision'],
    information=information)
