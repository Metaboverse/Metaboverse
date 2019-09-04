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
import time
import shutil
import csv
import copy
import pickle
import numpy

"""Import internal dependencies
"""
from metabalyze.utils import progress_bar
from metabalyze.metabonet_network.convert import convert_metabolites_text
from metabalyze.metabonet_network.convert import convert_reactions_text
from metabalyze.metabonet_network.utils import match_hmdb_entries_by_identifiers_names
from metabalyze.metabonet_network.utils import collect_unique_elements
from metabalyze.metabonet_network.utils import collect_reaction_participants_value
from metabalyze.metabonet_network.utils import compare_lists_by_mutual_inclusion
from metabalyze.metabonet_network.utils import collect_value_from_records
from metabalyze.metabonet_network.utils import filter_common_elements
from metabalyze.metabonet_network.utils import find_index
from metabalyze.metabonet_network.utils import confirm_path_directory
from metabalyze.metabonet_network.utils import write_file_table
from metabalyze.metabonet_network.utils import prepare_curation_report

"""Reads and organizes source information from file
arguments:
    directory (str): directory of source files
returns:
    (object): source information
"""
def read_source(
        directory):

    # Specify directories and files.
    path_hmdb = directory + 'extraction/hmdb_summary.pickle'
    path = directory + 'collection/'
    path_compartments = path + 'compartments.pickle'
    path_processes = path + 'processes.pickle'
    path_reactions = path + 'reactions.pickle'
    path_metabolites = path + 'metabolites.pickle'

    # Read information from file.
    with open(path_hmdb, 'rb') as file_source:
        summary_hmdb = pickle.load(file_source)
    with open(path_compartments, 'rb') as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, 'rb') as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, 'rb') as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, 'rb') as file_source:
        metabolites = pickle.load(file_source)

    return {
        'summary_hmdb': summary_hmdb,
        'compartments': compartments,
        'processes': processes,
        'reactions': reactions,
        'metabolites': metabolites}

"""Enhances information about metabolites
arguments:
    metabolites_original (dict<dict>): information about metabolites
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (dict<dict>): information about metabolites
"""
def enhance_metabolites(
        metabolites_original=None,
        summary_hmdb=None):

    metabolites_novel = {}
    counter = 1
    total = len(metabolites_original)

    # Iterate on records for metabolites.
    for metabolite in metabolites_original.values():

        # Enhance information about metabolite.
        metabolite_novel = enhance_metabolite(
            metabolite_original=metabolite,
            summary_hmdb=summary_hmdb)

        # Compile information
        metabolites_novel[metabolite_novel['identifier']] = metabolite_novel

        # Continue counter
        progress_bar(
            counter,
            total,
            status='Enhance metabolites')
        counter += 1

    return metabolites_novel

"""Enhances information about a metabolite
arguments:
    metabolite_original (dict): information about a metabolite
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (dict): information about a metabolite
"""
def enhance_metabolite(
        metabolite_original=None,
        summary_hmdb=None):

    # Copy information.
    metabolite_novel = copy.deepcopy(metabolite_original)

    # Enhance metabolite's references.
    references_novel = enhance_metabolite_references(
        name=metabolite_novel['name'],
        references_original=metabolite_novel['references'],
        summary_hmdb=summary_hmdb)
    metabolite_novel['references'] = references_novel

    # Use name from HMDB.
    if len(references_novel['hmdb']) > 0:

        identifier_hmdb = references_novel['hmdb'][0]
        name = summary_hmdb[identifier_hmdb]['name']
        metabolite_novel['name'] = name

    return metabolite_novel

"""Enhances information about a metabolite by including references from HMDB
arguments:
    name (str): name of metabolite
    references_original (dict): references about a metabolite
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (dict): references about a metabolite
"""
def enhance_metabolite_references(
        name=None,
        references_original=None,
        summary_hmdb=None):

    # Copy information.
    references_novel = copy.deepcopy(references_original)

    # Enhance references to HMDB.
    references_hmdb_original = references_novel['hmdb']
    references_hmdb_novel = match_hmdb_entries_by_identifiers_names(
        identifiers=references_hmdb_original,
        names=[name],
        summary_hmdb=summary_hmdb)

    # Extract references from entries in HMDB
    hmdb_references = collect_hmdb_entries_references(
        keys=references_hmdb_novel,
        summary_hmdb=summary_hmdb)

    # Combine supplemental references to original references
    references_novel['hmdb'] = collect_unique_elements(
        hmdb_references['hmdb'])
    references_novel['pubchem'] = collect_unique_elements(
        references_original['pubchem'] + hmdb_references['pubchem'])
    references_novel['chebi'] = collect_unique_elements(
        references_original['chebi'] + hmdb_references['chebi'])
    references_novel['kegg'] = collect_unique_elements(
        references_original['kegg'] + hmdb_references['kegg'])

    return references_novel

"""Extracts references from entries in HMDB
arguments:
    keys (list<str>): keys of entries in HMDB
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (dict<list<str>>): references
"""
def collect_hmdb_entries_references(
        keys=None,
        summary_hmdb=None):

    references_hmdb = []
    references_pubchem = []
    references_chebi = []
    references_kegg = []

    for key in keys:

        # Only include primary accession identifiers in HMDB collection.
        record = summary_hmdb[key]
        references_hmdb.append(record['identifier'])

        # Only include valid identifiers in the collection.
        pubchem = record['reference_pubchem']
        chebi = record['reference_chebi']
        kegg = record['reference_kegg']

        if (pubchem is not None) and (len(pubchem) > 0):

            references_pubchem.append(pubchem)

        if (chebi is not None) and (len(chebi) > 0):

            references_chebi.append(chebi)

        if (kegg is not None) and (len(kegg) > 0):

            references_kegg.append(kegg)

    return {
        'hmdb': references_hmdb,
        'pubchem': references_pubchem,
        'chebi': references_chebi,
        'kegg': references_kegg}

"""Includes information about reactions' behavior
arguments:
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict>): information about reactions
"""
def include_reactions_behaviors(
        reactions_original=None):

    reactions_novel = {}
    counter = 1
    total = len(reactions_original)

    for reaction_original in reactions_original.values():

        reaction_novel = include_reaction_behavior(
            reaction_original=reaction_original)
        reactions_novel[reaction_novel['identifier']] = reaction_novel

        progress_bar(
            counter,
            total,
            status='Include reaction behavior')
        counter += 1

    return reactions_novel

"""Includes information about a reaction's behavior
arguments:
    reaction_original (dict): information about a reaction
returns:
    (dict): information about a reaction
"""
def include_reaction_behavior(
        reaction_original=None):

    # Determine whether reaction involves chemical conversion between reactants
    # and products
    conversion = determine_reaction_conversion(
        reaction=reaction_original)

    # Determine whether reaction involves participation of metabolites in
    # multiple compartments
    dispersal = determine_reaction_dispersal(
        reaction=reaction_original)

    # Determine reaction's transports
    transports = collect_reaction_transports(
        reaction=reaction_original)

    # Determine whether reaction involves transport of metabolites between
    # compartments
    transport = len(transports) > 0

    # Compile information
    reaction_novel = copy.deepcopy(reaction_original)
    reaction_novel['conversion'] = conversion
    reaction_novel['dispersal'] = dispersal
    reaction_novel['transports'] = transports
    reaction_novel['transport'] = transport

    return reaction_novel

"""Determines whether a reaction involves chemical conversion
arguments:
    reaction (dict): information about a reaction
returns:
    (bool): whether the reaction involves chemical conversion
"""
def determine_reaction_conversion(
        reaction=None):

    reactant_metabolites = collect_reaction_participants_value(
        key='metabolite',
        criteria={'roles': ['reactant']},
        participants=reaction['participants'])
    product_metabolites = collect_reaction_participants_value(
        key='metabolite',
        criteria={'roles': ['product']},
        participants=reaction['participants'])

    return not compare_lists_by_mutual_inclusion(
        list_one=reactant_metabolites,
        list_two=product_metabolites)

"""Determines whether a reaction involves metabolites in multiple compartments
arguments:
    reaction (dict): information about a reaction
returns:
    (bool): whether the reaction involves participation of metabolites in
        multiple compartments
"""
def determine_reaction_dispersal(
        reaction=None):

    compartments = collect_value_from_records(
        key='compartment',
        records=reaction['participants'])

    return len(compartments) > 0

"""Collects information about a reaction's transports.
This procedure applies an overly restrictive definition of transport that
requires chemically-identical metabolites in two separate compartments.
Some transports involve chemical conversion of substrates as part of
transport.
arguments:
    reaction (dict): information about a reaction
returns:
    (list<dict>): information about a reaction's transports
"""
def collect_reaction_transports(
        reaction=None):

    metabolites_reactant = collect_reaction_participants_value(
        key='metabolite',
        criteria={'roles': ['reactant']},
        participants=reaction['participants'])

    metabolites_product = collect_reaction_participants_value(
        key='metabolite',
        criteria={'roles': ['product']},
        participants=reaction['participants'])

    # Collect metabolites that participate as both reactants and products
    metabolites = filter_common_elements(
        list_one=metabolites_product,
        list_two=metabolites_reactant)
    transports = []

    for metabolite in metabolites:

        # Determine metabolite's compartments as reactant and product
        compartments_reactant = collect_reaction_participants_value(
            key='compartment',
            criteria={
                'metabolites': [metabolite],
                'roles': ['reactant']},
            participants=reaction['participants'])
        compartments_product = collect_reaction_participants_value(
            key='compartment',
            criteria={
                'metabolites': [metabolite],
                'roles': ['product']},
            participants=reaction['participants'])

        # Determine whether there is a difference between the metabolite's
        # compartments as reactant and product
        transport = not compare_lists_by_mutual_inclusion(
            list_one=compartments_reactant,
            list_two=compartments_product)

        if transport:

            compartments = compartments_reactant + compartments_product
            compartments_unique = collect_unique_elements(
                elements_original=compartments)
            record = {
                'metabolite': metabolite,
                'compartments': compartments_unique}
            transports.append(record)

    return transports

"""Includes information about reactions' transport processes
arguments:
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict>): information about reactions
"""
def include_reactions_transport_processes(
        reactions_original=None):

    # Collect information about compartments in which each metabolite
    # participates in each process
    processes_dispersal = collect_processes_metabolites_compartments(
        reactions=reactions_original)

    # Filter information for prospective transports of metabolites between
    # compartments in processes
    processes_transports = filter_processes_transports(
        processes_dispersal=processes_dispersal)

    reactions_novel = {}
    counter = 1
    total = len(reactions_original)

    for reaction_original in reactions_original.values():

        reaction_novel = include_reaction_transport_processes(
            reaction_original=reaction_original,
            processes_transports=processes_transports)
        reactions_novel[reaction_novel["identifier"]] = reaction_novel

        progress_bar(
            counter,
            total,
            status='Including reactions and transport processes')
        counter += 1

    return reactions_novel

"""Collects information about processes' metabolites and compartments
arguments:
    reactions (dict<dict>): information about reactions
returns:
    (dict<dict<list<str>>>): information about compartments in which each
        metabolite participates in each process
"""
def collect_processes_metabolites_compartments(
        reactions=None):

    collection = {}

    for reaction in reactions.values():

        processes = reaction['processes']

        for process in processes:

            if process not in collection:

                collection[process] = {}

            metabolites = collect_reaction_participants_value(
                key='metabolite',
                criteria={},
                participants=reaction['participants'])

            for metabolite in metabolites:

                if metabolite not in collection[process]:

                    collection[process][metabolite] = []

                compartments = collect_reaction_participants_value(
                    key='compartment',
                    criteria={'metabolites': [metabolite]},
                    participants=reaction['participants'])

                for compartment in compartments:

                    if compartment not in collection[process][metabolite]:

                        collection[process][metabolite].append(compartment)

    return collection

"""Collects information about transports in processes
arguments:
    processes_dispersal (dict<dict<list<str>>>): information about
        compartments in which each metabolite participates in each process
returns:
    (dict<dict<list<str>>>): information about transports in processes
"""
def filter_processes_transports(
        processes_dispersal=None):

    collection = {}

    for process in processes_dispersal.keys():

        collection[process] = {}

        for metabolite in processes_dispersal[process].keys():

            compartments = processes_dispersal[process][metabolite]

            if len(compartments) > 1:

                collection[process][metabolite] = copy.copy(compartments)

    return collection

"""Includes information about a reaction's transport processes
arguments:
    reaction_original (dict): information about a reaction
    processes_transports (dict<dict<list<str>>>): information about
        transports in processes
returns:
    (dict): information about a reaction
"""
def include_reaction_transport_processes(
        reaction_original=None,
        processes_transports=None):

    # Determine processes in which reaction participates by transport
    processes_transport = collect_reaction_transport_processes(
        reaction=reaction_original,
        processes_transports=processes_transports)
    processes_original = reaction_original['processes']
    processes_total = processes_original + processes_transport
    processes_unique = collect_unique_elements(
        elements_original=processes_total)

    # Compile information
    reaction_novel = copy.deepcopy(reaction_original)
    reaction_novel['processes'] = processes_unique

    return reaction_novel

"""Collects processes in which a reaction participates by transport
arguments:
    reaction (dict): information about a reaction
    processes_transports (dict<dict<list<str>>>): information about
        transports in processes
returns:
    (list<str>): identifiers of processes
"""
def collect_reaction_transport_processes(
        reaction=None,
        processes_transports=None):

    transports_reaction = reaction['transports']
    processes_transport = []

    for process in processes_transports.keys():

        # Determine whether reaction's transports match any of the process'
        # metabolites and compartments
        metabolites_process = processes_transports[process]

        for transport_reaction in transports_reaction:

            metabolite_reaction = transport_reaction['metabolite']
            compartments_reaction = transport_reaction['compartments']

            if metabolite_reaction in metabolites_process.keys():

                compartments_process = metabolites_process[metabolite_reaction]

                # Determine whether multiple compartments match between the
                # reaction and the process
                compartments = filter_common_elements(
                    list_one=compartments_reaction,
                    list_two=compartments_process)

                if len(compartments) > 1:

                    # Reaction participates in the process by transport
                    processes_transport.append(process)

    return processes_transport

"""Includes information about reactions' replications.
Replicate reactions involve participation of reactants and products that
are identical metabolites but not necessarily identical compartments.
Consideration of replication is necessary to avoid redundancy when
compartments are irrelevant.
arguments:
    reactions_original (dict<dict>): information about reactions
returns:
    (dict<dict>): information about reactions
"""
def include_reactions_replications(
        reactions_original=None):

    # Collect information about replicate reactions
    reactions_replicates = collect_reactions_replicates(
        reactions=reactions_original)
    reactions_novel = {}
    counter = 1
    total = len(reactions_original)

    for reaction_original in reactions_original.values():

        reaction_novel = include_reaction_replication(
            reaction_original=reaction_original,
            reactions_replicates=reactions_replicates)
        reactions_novel[reaction_novel['identifier']] = reaction_novel

        progress_bar(
            counter,
            total,
            status='Collect reactions replicates')
        counter += 1

    return reactions_novel

"""Collects information about reactions' replications
arguments:
    reactions (dict<dict>): information about reactions
returns:
    (list<dict>): information about reactions' replications
"""
def collect_reactions_replicates(
        reactions=None):

    reactions_replicates = []

    for reaction in reactions.values():

        identifier = reaction['identifier']

        # Collect identifiers of metabolites that participate as reactants and
        # products in the reaction
        reactants = collect_reaction_participants_value(
            key='metabolite',
            criteria={'roles': ['reactant']},
            participants=reaction['participants'])
        products = collect_reaction_participants_value(
            key='metabolite',
            criteria={'roles': ['product']},
            participants=reaction['participants'])

        # Determine whether collection includes a record for an identical
        # combination of reactants and products
        index = find_index_reactions_replicates_reactants_products(
            reactants=reactants,
            products=products,
            reactions_replicates=reactions_replicates)

        if index == -1:

            # Record does not exist
            # Create novel record
            record = {
                'reactions': [reaction['identifier']],
                'reactants': reactants,
                'products': products}
            reactions_replicates.append(record)

        else:

            # Record exists
            # Include reaction in record
            reactions_replicates[index]['reactions'].append(identifier)

    return reactions_replicates

"""Finds index of a record for replicate reactions by reactants and products
arguments:
    reactants (list<str>): identifiers of metabolites that participate in a
        reaction as reactants
    products (list<str>): identifiers of metabolites that participate in a
        reaction as products
    reactions_replicates (list<dict>): information about reactions'
        replications
returns:
    (int): index of record if it exists or -1 if it does not exist
"""
def find_index_reactions_replicates_reactants_products(
        reactants=None,
        products=None,
        reactions_replicates=None):

    def match (
            record=None):

        match_reactants = compare_lists_by_mutual_inclusion(
            list_one=reactants,
            list_two=record['reactants'])
        match_products = compare_lists_by_mutual_inclusion(
            list_one=products,
            list_two=record['products'])

        return match_reactants and match_products

    return find_index(match, reactions_replicates)

"""Includes information about a reaction's replication
arguments:
    reaction_original (dict): information about a reaction
    reactions_replicates (list<dict>): information about reactions'
        replications
returns:
    (dict): information about a reaction
"""
def include_reaction_replication(
        reaction_original=None,
        reactions_replicates=None):

    # Determine replicate reactions
    index = find_index_reactions_replicates_identifier(
        identifier=reaction_original['identifier'],
        reactions_replicates=reactions_replicates)

    if index == -1:

        replicates = []

    else:

        replicates = reactions_replicates[index]['reactions']

    # Compile information
    reaction_novel = copy.deepcopy(reaction_original)
    reaction_novel['replicates'] = replicates
    reaction_novel['replication'] = len(replicates) > 1

    return reaction_novel

"""Finds index of a record for replicate reactions by identifier
arguments:
    identifier (str): identifier of a reaction
    reactions_replicates (list<dict>): information about reactions'
        replications
returns:
    (int): index of record if it exists or -1 if it does not exist
"""
def find_index_reactions_replicates_identifier(
        identifier=None,
        reactions_replicates=None):

    def match(
            record=None):

        return identifier in record['reactions']

    return find_index(match, reactions_replicates)

"""Filters reactions by relevance to contextual metabolic network.
arguments:
    reactions (dict<dict>): information about reactions
returns:
    (list<dict>): information about reactions
"""
def filter_reactions(
        reactions_original=None):

    reactions_filter = []
    counter = 1
    total = len(reactions_original)

    for reaction in reactions_original.values():

        match, record = filter_reaction(reaction=reaction)

        if match:

            reactions_filter.append(record)

        progress_bar(counter, total, status='Filter reactions')
        counter += 1

    return reactions_filter

"""Filters a reaction by relevance to contextual metabolic network.
arguments:
    reaction (dict): information about a reaction
returns:
    (tuple<bool, dict>): whether reaction passes filters, and report
"""
def filter_reaction(
        reaction=None):

    # Name.
    # Biomass, protein assembly, and protein degradation are irrelevant.
    name = (
        (reaction['name'] == 'Generic human biomass reaction') \
        or (reaction['name'] == 'Protein degradation') \
        or ('Protein assembly' in reaction['name']))

    # Compartment.
    # Boundary is irrelevant.
    compartments = collect_value_from_records(
        key='compartment',
        records=reaction['participants'])
    compartment_boundary = ('BOUNDARY' in compartments)

    # Extracellular region is irrelevant.
    compartment_exterior = ('MNXC2' in compartments)

    # Metabolite.
    # Biomass is irrelevant.
    metabolites = collect_value_from_records(
        key='metabolite',
        records=reaction['participants'])
    metabolite = ('BIOMASS' in metabolites)

    # Reference.
    # MetaNetX reaction MNXR01 is for a meaningless proton exchange.
    reference = 'MNXR01' in reaction['references']['metanetx']

    # Determine whether reaction passes filters.
    filter = (
        name \
        or compartment_boundary \
        or compartment_exterior \
        or metabolite \
        or reference)

    # Prepare report.
    record = {
        'identifier_original': reaction['identifier'],
        'identifier_novel': 'null',
        'name_original': reaction['name'],
        'name_novel': reaction['name'],
        'custom': False}

    return (filter, record)

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
    path_compartments = os.path.join(directory, 'compartments.pickle')
    path_processes = os.path.join(directory, 'processes.pickle')
    path_reactions = os.path.join(directory, 'reactions.pickle')
    path_metabolites = os.path.join(directory, 'metabolites.pickle')
    path_metabolites_report = os.path.join(directory, 'metabolites.tsv')
    path_reactions_report = os.path.join(directory, 'reactions.tsv')
    path_reactions_filter = os.path.join(directory, 'reactions_filter.tsv')

    # Write information to file.
    with open(path_compartments, 'wb') as file_product:
        pickle.dump(information['compartments'], file_product)
    with open(path_processes, 'wb') as file_product:
        pickle.dump(information['processes'], file_product)
    with open(path_reactions, 'wb') as file_product:
        pickle.dump(information['reactions'], file_product)
    with open(path_metabolites, 'wb') as file_product:
        pickle.dump(information['metabolites'], file_product)

    write_file_table(
        information=information['metabolites_report'],
        path_file=path_metabolites_report,
        names=information['metabolites_report'][0].keys(),
        delimiter='\t')
    write_file_table(
        information=information['reactions_report'],
        path_file=path_reactions_report,
        names=information['reactions_report'][0].keys(),
        delimiter='\t')
    write_file_table(
        information=information['reactions_filter'],
        path_file=path_reactions_filter,
        names=information['reactions_filter'][0].keys(),
        delimiter='\t')

"""Function to execute module's main behavior.
The purpose of this procedure is to curate information about metabolic
entities and sets.
arguments:
    directory (str): path to directory for source and product files
"""
def execute_procedure(
        args_dict):

    # Read source information from file.
    source = read_source(
        args_dict['source'])

    # Enhance metabolites' references.
    print('Step 0/8: Starting enhancement. This may take a while...')
    print('Step 1/8: Enhancing metabolites...')
    metabolites = enhance_metabolites(
        metabolites_original=source["metabolites"],
        summary_hmdb=source["summary_hmdb"])

    # Include information about reactions' behavior.
    print('Step 2/8: Including reactions\' behaviors...')
    reactions_behavior = include_reactions_behaviors(
        reactions_original=source["reactions"])

    # Include transport reactions in processes.
    print('Step 3/8: Including reactions transport processes...')
    reactions_process = include_reactions_transport_processes(
        reactions_original=reactions_behavior)

    # Include information about reactions' replicates.
    print('Step 4/8: Including reactions replications...')
    reactions_replication = include_reactions_replications(
        reactions_original=reactions_process)

    # Prepare reports of information for review.
    print('Step 5/8: Converting metabolites text...')
    convert_one = convert_metabolites_text
    metabolites_report = convert_one(
        metabolites=metabolites)

    print('Step 6/8: Converting reactions text...')
    convert_two = convert_reactions_text
    reactions_report = convert_two(
        reactions=reactions_replication)

    # Filter reactions.
    print('Step 7/8: Filtering reactions...')
    reactions_filter = filter_reactions(
        reactions_original=reactions_replication)

    # Compile information.
    information = {
        "compartments": source["compartments"],
        "processes": source["processes"],
        "metabolites": metabolites,
        "reactions": reactions_replication,
        "metabolites_report": metabolites_report,
        "reactions_report": reactions_report,
        "reactions_filter": reactions_filter}

    #Write product information to file.
    print('Step 8/8: Writing outputs...')
    write_product(
        args_dict['enhance'],
        information=information)

    # Report.
    report = prepare_curation_report(
        compartments=source["compartments"],
        processes=source["processes"],
        reactions=reactions_replication,
        metabolites=metabolites)
    print(report)
