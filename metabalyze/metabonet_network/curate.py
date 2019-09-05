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
import sys
import csv
import copy
import pickle
import numpy as np
from functools import partial

"""Import internal dependencies
"""
from metabalyze.metabonet_network.model import convert_metabolites_text
from metabalyze.metabonet_network.model import convert_reactions_text
from metabalyze.metabonet_network.utils import confirm_path_directory
from metabalyze.metabonet_network.utils import read_file_table
from metabalyze.metabonet_network.utils import collect_unique_elements
from metabalyze.metabonet_network.utils import collect_value_from_records
from metabalyze.metabonet_network.utils import prepare_curation_report
from metabalyze.utils import progress_bar
from metabalyze.utils import get_cores
from metabalyze.utils import run_chunks

"""Set globals
"""
__customization__  =  os.path.dirname(os.path.realpath(__file__)) + '/customization/'

custom_compartments = 'curation_compartments.tsv'
custom_processes = 'curation_processes.tsv'
custom_reactions = 'curation_reactions.tsv'
custom_metabolites = 'curation_metabolites.tsv'

hmdb_pickle = 'hmdb_summary.pickle'
compartment_pickle = 'compartments.pickle'
process_pickle = 'processes.pickle'
reaction_pickle = 'reactions.pickle'
metabolite_pickle = 'metabolites.pickle'

"""Reads and organizes source information from file
arguments:
    directory (str): directory of source files
returns:
    (object): source information
"""
def read_source(
        directory):

    # Specify directories and files
    path_compartments_curation = __customization__ + custom_compartments
    path_processes_curation = __customization__ + custom_processes
    path_reactions_curation = __customization__ + custom_reactions
    path_metabolites_curation = __customization__ + custom_metabolites

    path_compartments = directory + compartment_pickle
    path_processes = directory + process_pickle
    path_reactions = directory + reaction_pickle
    path_metabolites = directory + metabolite_pickle

    # Read information from file.
    with open(path_compartments, 'rb') as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, 'rb') as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, 'rb') as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, 'rb') as file_source:
        metabolites = pickle.load(file_source)

    compartments_curation = read_file_table(
        path_file=path_compartments_curation,
        names=None,
        delimiter='\t')
    processes_curation = read_file_table(
        path_file=path_processes_curation,
        names=None,
        delimiter='\t')
    reactions_curation = read_file_table(
        path_file=path_reactions_curation,
        names=None,
        delimiter='\t')
    metabolites_curation = read_file_table(
        path_file=path_metabolites_curation,
        names=None,
        delimiter='\t')

    return {
        'compartments': compartments,
        'processes': processes,
        'reactions': reactions,
        'metabolites': metabolites,
        'compartments_curation': compartments_curation,
        'processes_curation': processes_curation,
        'reactions_curation': reactions_curation,
        'metabolites_curation': metabolites_curation}

"""Curates information about specific compartments and relevant reactions.
arguments:
    compartments_curation (list<dict<str>>): information to change about
        specific compartments
    compartments_original (dict<dict>): information about compartments
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict<dict>>): information about compartments and reactions
"""
def curate_compartments(
        compartments_curation,
        compartments_original,
        reactions_original):

    # Copy information
    compartments_novel = copy.deepcopy(compartments_original)
    reactions_novel = copy.deepcopy(reactions_original)

    for record in compartments_curation:

        # Interpretation.
        identifier_original = record['identifier_original']
        identifier_novel = record['identifier_novel']
        name_original = record['name_original']
        name_novel = record['name_novel']

        # Determine method to change information.
        match_identifiers = identifier_original == identifier_novel
        match_names = name_original == name_novel

        if identifier_novel == 'null':

            if identifier_original in compartments_novel:

                # Remove compartment.
                del compartments_novel[identifier_original]

                # Removal of a compartment justifies removal of any reactions
                # within that compartment.
                # Remove relevant reactions.
                removals = []

                for key, record_reaction in reactions_novel.items():

                    # Determine whether any of reaction's participants are in
                    # the compartment.
                    match = determine_reaction_compartment(
                        compartment=identifier_original,
                        reaction=record_reaction)

                    if match:
                        removals.append(key)

                for removal in removals:

                    del reactions_novel[removal]

        elif not match_names:

            # Change name.
            if identifier_original in compartments_novel:
                compartments_novel[identifier_original]['name'] = name_novel

        else:
            pass

    return {
        'compartments': compartments_novel,
        'reactions': reactions_novel}

"""Determines whether any of reaction's participants are in a compartment
arguments:
    compartment (str): identifier of a compartment
    reaction (dict): information about a reaction
returns:
    (bool): whether any of reaction's participants are in the compartment
"""
def determine_reaction_compartment(
        compartment,
        reaction):

    participants = reaction['participants']

    for participant in participants:

        if participant['compartment'] == compartment:
            return True

    return False

"""Curates information about specific processes and relevant reactions.
arguments:
    processes_curation (list<dict<str>>): information to change about
        specific processes
    processes_original (dict<dict>): information about processes
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict<dict>>): information about processes and reactions
"""
def curate_processes(
        processes_curation,
        processes_original,
        reactions_original):

    # Copy information
    processes_novel = copy.deepcopy(processes_original)
    reactions_novel = copy.deepcopy(reactions_original)

    for record in processes_curation:

        # Interpretation
        identifier_original = record['identifier_original']
        identifier_novel = record['identifier_novel']
        name_original = record['name_original']
        name_novel = record['name_novel']

        # Determine method to change information
        match_identifiers = identifier_original == identifier_novel
        match_names = name_original == name_novel

        if identifier_novel == 'null':

            if identifier_original in processes_novel:
                del processes_novel[identifier_original] # remove process

                # carry-over from original:
                # TODO: also remove the process from any reactions' processes
                # Removal of a process does not justify removal of any
                # reactions that participate in that process.

        else:

            if not match_identifiers:

                # Change identifier
                # Remove original
                if identifier_original in processes_novel:
                    del processes_novel[identifier_original]

                # Replace with novel
                if identifier_novel in processes_novel:

                    for reaction in reactions_novel.values():

                        processes = reaction['processes']

                        if identifier_original in processes:

                            for index, process in enumerate(processes):

                                if process == identifier_original:
                                    processes[index] = identifier_novel

                            # Collect unique values.
                            processes_unique = collect_unique_elements(
                                processes)
                            reaction['processes'] = processes_unique
                            reactions_novel[reaction['identifier']] = reaction

            if not match_names:

                # Change name
                if identifier_novel in processes_novel:
                    processes_novel[identifier_novel]['name'] = name_novel

    return {
        'processes': processes_novel,
        'reactions': reactions_novel}

"""Curates information about specific metabolites and relevant reactions.
arguments:
    metabolites_curation (list<dict<str>>): information to change about
        specific metabolites
    metabolites_original (dict<dict>): information about metabolites
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict<dict>>): information about metabolites and reactions
"""
def curate_metabolites(
        metabolites_curation,
        metabolites_original,
        reactions_original):

    # Copy information
    metabolites_novel = copy.deepcopy(metabolites_original)
    reactions_novel = copy.deepcopy(reactions_original)

    for record in metabolites_curation:

        # Interpretation
        identifier_original = record['identifier_original']
        identifier_novel = record['identifier_novel']
        name_original = record['name_original']
        name_novel = record['name_novel']

        # Determine method to change information
        if identifier_novel == 'null':

            if identifier_original in metabolites_novel:
                del metabolites_novel[identifier_original] # Remove metabolite

                # Remove metabolite from relevant reactions
                reactions_novel = change_reactions_participants_metabolite(
                    reactions_original=reactions_novel,
                    metabolite_original=identifier_original,
                    metabolite_novel='null',
                    remove=True,
                    replace=False)

        else:

            if not (identifier_original == identifier_novel):

                # Change identifier.
                # Change identifier in reactions' participants.
                if identifier_original in metabolites_novel \
                and identifier_novel not in metabolites_novel:

                    # Copy original record.
                    metabolite_novel = copy.deepcopy(metabolites_novel[identifier_original])

                    # Change identifier.
                    metabolite_novel['identifier'] = identifier_novel

                    # Replace original record with novel record.
                    del metabolites_novel[identifier_original]
                    metabolites_novel[identifier_novel] = metabolite_novel

                elif identifier_original in metabolites_novel \
                and identifier_novel in metabolites_novel:

                    # Remove original record.
                    del metabolites_novel[identifier_original]

                # Replace metabolite in relevant reactions.
                reactions_novel = change_reactions_participants_metabolite(
                    reactions_original=reactions_novel,
                    metabolite_original=identifier_original,
                    metabolite_novel=identifier_novel,
                    remove=False,
                    replace=True)

            if not (name_original == name_novel):

                # Change name.
                if identifier_novel in metabolites_novel:
                    metabolites_novel[identifier_novel]['name'] = name_novel

            # Curate metabolite's references.
            metabolites_novel = curate_metabolites_references(
                metabolite_curation=record,
                metabolites=metabolites_novel)

    return {
        'metabolites': metabolites_novel,
        'reactions': reactions_novel}

"""Curates references for a metabolite.
arguments:
    metabolite_curation (dict<str>): information to change about a
        metabolite
    metabolites (dict<dict>): information about metabolites
returns:
    (dict<dict>): information about metabolites
"""
def curate_metabolites_references(
        metabolite_curation,
        metabolites):

    def match_items(
            match_list):

        try:
            return match_list[0] != match_list[1]
        except:
            return True

    # Interpretation
    identifier = metabolite_curation['identifier_novel']
    hmdb_novel = metabolite_curation['hmdb_novel']
    hmdb_error = metabolite_curation['hmdb_error']
    pubchem_novel = metabolite_curation['pubchem_novel']
    pubchem_error = metabolite_curation['pubchem_error']

    if identifier in metabolites.keys():

        if len(hmdb_novel) > 0:
            references_original = metabolites[identifier]['references']['hmdb']
            references_original.append(hmdb_novel)
            references_novel = collect_unique_elements(
                references_original)
            metabolites[identifier]['references']['hmdb'] = references_novel

        if len(hmdb_error) > 0:
            references_original = metabolites[identifier]['references']['hmdb']
            references_novel = list(filter(match_items, [references_original, hmdb_error]))
            metabolites[identifier]['references']['hmdb'] = references_novel

        if len(pubchem_novel) > 0:
            references_original = metabolites[identifier]['references']['pubchem']
            references_original.append(pubchem_novel)
            references_novel = collect_unique_elements(
                references_original)
            metabolites[identifier]['references']['pubchem'] = references_novel

        if len(pubchem_error) > 0:
            references_original = metabolites[identifier]['references']['pubchem']
            references_novel = list(filter(match_items, [references_original, pubchem_error]))
            metabolites[identifier]['references']['pubchem'] = references_novel

    return metabolites

"""Changes metabolite in reactions' participants.
arguments:
    reactions_original (dict<dict>): information about reactions
    metabolite_original (str): identifier of a metabolite
    metabolite_novel (str): identifier of a metabolite
    remove (bool): whether to remove the metabolite from participants
    replace (bool): whether to replace the metabolite in participants
returns:
    (dict<dict>): information about reactions
"""
def change_reactions_participants_metabolite(
        reactions_original,
        metabolite_original,
        metabolite_novel,
        remove,
        replace):

    reactions_novel = {}

    for reaction in reactions_original.values():

        participants_original = reaction['participants']
        participants_novel = []

        for party in participants_original:

            if party['metabolite'] != metabolite_original:
                participants_novel.append(party)

            else:

                if remove:
                    pass # Omit participant

                elif replace:
                    party['metabolite'] = metabolite_novel
                    participants_novel.append(party)

        reaction['participants'] = participants_novel
        reactions_novel[reaction['identifier']] = reaction

    return reactions_novel

"""Curates information about specific reactions.
arguments:
    reactions_curation (list<dict<str>>): information to change about
        specific reactions
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict>): information about reactions
"""
def run_curate_reactions(
        reactions_curation_chunk,
        reactions_original,
        reactions_novel):

    counter = 1
    total = len(reactions_curation_chunk) + 1

    for record in reactions_curation_chunk:

        # Interpretation.
        identifier_original = record['identifier_original']
        identifier_novel = record['identifier_novel']
        name_original = record['name_original']
        name_novel = record['name_novel']

        # Determine method to change information.
        match_identifiers = identifier_original == identifier_novel
        match_names = name_original == name_novel

        if identifier_novel == 'null':

            if identifier_original in reactions_novel:
                del reactions_novel[identifier_original] # Remove reaction

        elif not match_names:

            # Change name
            if identifier_original in reactions_novel:
                reactions_novel[identifier_original]['name'] = name_novel

        # Filter references to replicate reactions.
        # Ensure that all references to reactions are valid.
        reactions_replicates = filter_reaction_replicates(
            reactions_original=reactions_novel)

        progress_bar(
            counter,
            total,
            status='Curate reactions')
        counter += 1

    return reactions_replicates

def curate_reactions(
        reactions_curation,
        reactions_original,
        args_dict):

    # Copy information
    reactions_novel = copy.deepcopy(reactions_original)

    # Get cores and chunks
    cores = get_cores(args_dict)

    # Chunk data based on core #
    chunks = np.array_split(
        reactions_curation,
        cores)
    func = partial(
        run_curate_reactions,
        reactions_original=reactions_original,
        reactions_novel=reactions_novel)
    reactions_replicates = run_chunks(
        func,
        chunks,
        cores)

    return reactions_replicates


"""Filters references to replicate reactions.
arguments:
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict>): information about reactions
"""
def filter_reaction_replicates(
        reactions_original):

    reactions_novel = copy.deepcopy(reactions_original)

    def match(
            identifier):

        return identifier in reactions_novel.keys()

    for key in reactions_novel.keys():

        reaction = reactions_novel[key]
        replicates_original = reaction["replicates"]
        replicates_novel = list(filter(match, replicates_original))
        reactions_novel[key]["replicates"] = replicates_novel

    return reactions_novel

"""Accesses summary information about reactions of interest.
arguments:
    reactions_interest (list<dict<str>>): identifiers of reactions of
        interest
    reactions (dict<dict>): information about reactions
    directory (str): path to directory for source and product files
returns:
    (list<dict<str>>): information about measurements and signals for all
        samples
"""
def access_reactions_summary(
        reactions_interest,
        reactions,
        directory):

    # Collect information about reactions of interest
    identifiers = collect_value_from_records(
        key='identifier',
        records=reactions_interest)

    identifiers_unique = collect_unique_elements(identifiers)
    reactions_summary = []

    for identifier in identifiers_unique:

        reaction = reactions[identifier]
        name = reaction['name']
        metanetx = ';'.join(reaction['references']['metanetx'])
        gene = ';'.join(reaction['references']['gene'])
        enzyme = ';'.join(reaction['references']['enzyme'])
        record_product = {
            'identifier': identifier,
            'name': name,
            'metanetx': metanetx,
            'gene': gene,
            'enzyme': enzyme}
        reactions_summary.append(record_product)

    return reactions_summary

"""Writes product information to file
arguments:
    directory (str): directory for product files
    information (object): information to write to file
"""
def write_product(
        directory,
        information):

    # Specify directories and files
    confirm_path_directory(directory)
    path_compartments = directory + "compartments.pickle"
    path_processes = directory + "processes.pickle"
    path_reactions = directory + "reactions.pickle"
    path_reactions_summary = directory + "reactions_summary.tsv"
    path_metabolites = directory + "metabolites.pickle"
    path_metabolites_report = directory + "metabolites.tsv"
    path_reactions_report = directory + "reactions.tsv"

    # Write information to file
    with open(path_compartments, "wb") as file_product:
        pickle.dump(information["compartments"], file_product)
    with open(path_processes, "wb") as file_product:
        pickle.dump(information["processes"], file_product)
    with open(path_reactions, "wb") as file_product:
        pickle.dump(information["reactions"], file_product)
    with open(path_metabolites, "wb") as file_product:
        pickle.dump(information["metabolites"], file_product)

    try:
        write_file_table(
            information=information["metabolites_report"],
            path_file=path_metabolites_report,
            names=information["metabolites_report"][0].keys(),
            delimiter="\t")
        write_file_table(
            information=information["reactions_report"],
            path_file=path_reactions_report,
            names=information["reactions_report"][0].keys(),
            delimiter="\t")
    except:
        write_file_table(
            information=information["reactions_summary"],
            path_file=path_reactions_summary,
            names=information["reactions_summary"][0].keys(),
            delimiter="\t")

"""Function to execute module's main behavior.
The purpose of this procedure is to curate information about metabolic
entities and sets.
arguments:
    directory (str): path to directory for source and product files
"""
def __main__(
        args_dict):

    # Read source information from file
    print('Step 0/6: Reading in source data...')
    source = read_source(
        directory=args_dict['enhance'])

    # Change procedures allow custom changes to metabolites and reactions
    # Curate information about compartments
    print('Step 1/6: Curating compartments...')
    compartments_reactions = curate_compartments(
        compartments_curation=source['compartments_curation'],
        compartments_original=source['compartments'],
        reactions_original=source['reactions'])

    # Curate information about processes
    print('Step 2/6: Curating processes...')
    processes_reactions = curate_processes(
        processes_curation=source['processes_curation'],
        processes_original=source['processes'],
        reactions_original=compartments_reactions['reactions'])

    # Curate information about metabolites
    print('Step 3/6: Curating metabolites...')
    metabolites_reactions = curate_metabolites(
        metabolites_curation=source['metabolites_curation'],
        metabolites_original=source['metabolites'],
        reactions_original=processes_reactions['reactions'])

    # Curate information about reactions
    print('Step 4/6: Curating reactions...')
    reactions = curate_reactions(
        reactions_curation=source['reactions_curation'],
        reactions_original=metabolites_reactions['reactions'],
        args_dict=args_dict)

    # Prepare reports of information for review
    print('\nStep 5/6: Converting text...')
    metabolites_report = convert_metabolites_text(
        metabolites=metabolites_reactions['metabolites'])
    reactions_report = convert_reactions_text(
        reactions=reactions)

    # Compile information
    information = {
        'compartments': compartments_reactions['compartments'],
        'processes': processes_reactions['processes'],
        'metabolites': metabolites_reactions['metabolites'],
        'reactions': reactions,
        'metabolites_report': metabolites_report,
        'reactions_report': reactions_report}

    #Write product information to file
    print('Step 6/6: Writing output...\n')
    write_product(
        directory=args_dict['curate'],
        information=information)

    # Report
    report = prepare_curation_report(
        compartments=compartments_reactions['compartments'],
        processes=processes_reactions['processes'],
        reactions=reactions,
        metabolites=metabolites_reactions['metabolites'])
    print(report)
