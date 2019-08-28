"""License Information
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

"""Remove file if it exists
"""
def remove_file(
        path=None):

    if os.path.exists(path):
        os.remove(path)

"""Remove directory if it exists
"""
def remove_directory(
        path=None):

    if (os.path.exists(path)) and (len(os.listdir(path)) < 1):
        os.rmdir(path)

"""Confirms that a path to a directory exists.
Creates a directory if it does not already exist.
arguments:
    path (str): path to directory
"""
def confirm_path_directory(
        path=None):

    if not os.path.exists(path):

        os.makedirs(path)

"""Prepares a summary report on curation of metabolic sets and entities.
arguments:
    compartments (dict<dict>): information about compartments
    processes (dict<dict>): information about processes
    reactions (dict<dict>): information about reactions
    metabolites (dict<dict>): information about metabolites
returns:
    (str): report of summary information
"""
def prepare_curation_report(
        compartments=None,
        processes=None,
        reactions=None,
        metabolites=None):

    # Count compartments.
    count_compartments = len(compartments)

    # Count processes.
    count_processes = len(processes)

    # Count reactions.
    count_reactions = len(reactions)

    # Count metabolites.
    count_metabolites = len(metabolites)

    # Count reactions with references to MetaNetX.
    count_one = count_entities_with_references(
        references=['metanetx'],
        entities=reactions)
    proportion_one = count_one / count_reactions
    percentage_one = round((proportion_one * 100), 2)

    # Count reactions with references either to genes or enzyme commission.
    count_two = count_entities_with_references(
        references=[
            'gene',
            'enzyme'],
        entities=reactions)
    proportion_two = count_two / count_reactions
    percentage_two = round((proportion_two * 100), 2)

    # Count metabolites with references to MetaNetX.
    count_three = count_entities_with_references(
        references=['metanetx'],
        entities=metabolites)
    proportion_three = count_three / count_metabolites
    percentage_three = round((proportion_three * 100), 2)

    # Count metabolites with references to Human Metabolome Database (HMDB) and
    # PubChem.
    count_four = count_entities_with_references(
        references=[
            'hmdb',
            'pubchem'],
        entities=metabolites)
    proportion_four = count_four / count_metabolites
    percentage_four = round((proportion_four * 100), 2)

    # Compile information.
    report = textwrap.dedent("""\
        --------------------------------------------------
        curation report
        compartments: {count_compartments}
        processes: {count_processes}
        reactions: {count_reactions}
        metabolites: {count_metabolites}
        reactions in MetaNetX: {count_one} ({percentage_one} %)
        reactions with gene or enzyme: {count_two} ({percentage_two} %)
        metabolites in MetaNetX: {count_three} ({percentage_three} %)
        metabolites with HMDB or PubChem: {count_four} ({percentage_four} %)
        --------------------------------------------------
        """
        ).format(
            count_compartments=count_compartments,
            count_processes=count_processes,
            count_reactions=count_reactions,
            count_metabolites=count_metabolites,
            count_one=count_one,
            percentage_one=percentage_one,
            count_two=count_two,
            percentage_two=percentage_two,
            count_three=count_three,
            percentage_three=percentage_three,
            count_four=count_four,
            percentage_four=percentage_four)

    return report
