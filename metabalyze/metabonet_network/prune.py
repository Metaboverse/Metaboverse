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
import pickle
import copy

"""Import internal dependencies
"""
from metabalyze.metabonet_network.utils import read_file_table
from metabalyze.metabonet_network.utils import find_match
from metabalyze.metabonet_network.utils import collect_reaction_participants_value
from metabalyze.metabonet_network.utils import collect_records_targets_by_categories
from metabalyze.metabonet_network.utils import collect_unique_elements
from metabalyze.metabonet_network.utils import collect_value_from_records
from metabalyze.metabonet_network.utils import write_file_table

"""Set globals
"""
__customization__  =  os.path.dirname(os.path.realpath(__file__)) + '/customization/'

custom_compartments = 'filtration_compartments.tsv'
custom_processes = 'filtration_processes.tsv'
custom_reactions = 'simplification_reactions.tsv'
custom_metabolites = 'simplification_metabolites.tsv'

compartment_pickle = 'compartments.pickle'
process_pickle = 'processes.pickle'
reaction_pickle = 'reactions.pickle'
metabolite_pickle = 'metabolites.pickle'

metabolite_report = 'metabolites.tsv'
reaction_report = 'reactions.tsv'

"""Reads and organizes source information from file
arguments:
    directory (str): directory of source files
returns:
    (object): source information
"""
def read_source(
        directory):

    # Specify directories and files
    path_filtration_compartments = __customization__ + custom_compartments
    path_filtration_processes = __customization__ + custom_processes
    path_simplification_reactions = __customization__ + custom_reactions
    path_simplification_metabolites = __customization__ + custom_metabolites


    path_compartments = directory + compartment_pickle
    path_processes = directory + process_pickle
    path_reactions = directory + reaction_pickle
    path_metabolites = directory + metabolite_pickle

    # Read information from file
    with open(path_compartments, 'rb') as file_compartments:
        compartments = pickle.load(file_compartments)
    with open(path_processes, 'rb') as file_processes:
        processes = pickle.load(file_processes)
    with open(path_reactions, 'rb') as file_reactions:
        reactions = pickle.load(file_reactions)
    with open(path_metabolites, 'rb') as file_metabolites:
        metabolites = pickle.load(file_metabolites)

    filtration_compartments = read_file_table(
        path_file=path_filtration_compartments,
        names=None,
        delimiter='\t')
    filtration_processes = read_file_table(
        path_file=path_filtration_processes,
        names=None,
        delimiter='\t')
    simplification_reactions = read_file_table(
        path_file=path_simplification_reactions,
        names=None,
        delimiter='\t')
    simplification_metabolites = read_file_table(
        path_file=path_simplification_metabolites,
        names=None,
        delimiter='\t')

    return {
        'compartments': compartments,
        'processes': processes,
        'reactions': reactions,
        'metabolites': metabolites,
        'filtration_compartments': filtration_compartments,
        'filtration_processes': filtration_processes,
        'simplification_reactions': simplification_reactions,
        'simplification_metabolites': simplification_metabolites}

# Candidate reactions
"""Collects information about candidate reactions.
arguments:
    reactions (dict<dict>): information about reactions
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    filtration_processes (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific processes
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_reactions (list<dict<str>>): information about whether
        to simplify representations of specific reactions
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (dict<dict>): information about candidate reactions
"""
def collect_candidate_reactions(
        reactions,
        compartmentalization,
        filtration_compartments,
        filtration_processes,
        simplification,
        simplification_reactions,
        simplification_metabolites):

    reactions_candidacy = {}

    for reaction in reactions.values():

        candidacy, record = determine_reaction_candidacy(
            reaction_identifier=reaction['identifier'],
            reactions=reactions,
            reactions_candidacy=reactions_candidacy,
            compartmentalization=compartmentalization,
            filtration_compartments=filtration_compartments,
            filtration_processes=filtration_processes,
            simplification=simplification,
            simplification_reactions=simplification_reactions,
            simplification_metabolites=simplification_metabolites)

        if candidacy:
            reactions_candidacy[record['identifier']] = record

    return reactions_candidacy

"""
Determines whether a reaction is a candidate for representation in a
network.
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    reactions_candidacy (dict<dict>): information about candidate reactions
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    filtration_processes (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific processes
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_reactions (list<dict<str>>): information about whether
        to simplify representations of specific reactions
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (tuple<bool, dict>): whether reaction is a candidate, and information
        about the reaction candidate
"""
def determine_reaction_candidacy(
        reaction_identifier,
        reactions,
        reactions_candidacy,
        compartmentalization,
        filtration_compartments,
        filtration_processes,
        simplification,
        simplification_reactions,
        simplification_metabolites):

    # Relevance
    relevance = determine_reaction_relevance(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        filtration_processes=filtration_processes,
        simplification=simplification,
        simplification_reactions=simplification_reactions,
        simplification_metabolites=simplification_metabolites)

    # Redundancy
    redundancy = determine_reaction_redundancy(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        reactions_candidacy=reactions_candidacy,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        filtration_processes=filtration_processes,
        simplification=simplification,
        simplification_reactions=simplification_reactions,
        simplification_metabolites=simplification_metabolites)

    # Determine whether reaction is a candidate
    candidacy = relevance and not redundancy[0]

    # Prepare record for candidate reaction
    record = {
        'identifier': reaction_identifier,
        'reaction': reaction_identifier,
        'name': reactions[reaction_identifier]['name'],
        'replicates': redundancy[1],
        'reversibility': reactions[reaction_identifier]['reversibility']}

    return candidacy, record

"""Determines whether a reaction is relevant.
A reaction's relevance depends on the relevance of its behavior in the
context of the relevance of compartmentalization in general and the
relevance of specific compartments and processes.
A reaction's relevance also depends on the relevance of specific
metabolites that are its participants.
Factors in a reaction's relevance:
1. compartmentalization
2. filters by compartments
3. filters by processes
4. simplification of metabolites
5. simplification of reactions
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    filtration_processes (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific processes
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_reactions (list<dict<str>>): information about whether
        to simplify representations of specific reactions
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (bool): whether reaction is relevant
"""
def determine_reaction_relevance(
        reaction_identifier,
        reactions,
        compartmentalization,
        filtration_compartments,
        filtration_processes,
        simplification,
        simplification_reactions,
        simplification_metabolites):

    # Simplification
    simplification_match = determine_reaction_simplification(
        reaction_identifier=reaction_identifier,
        simplification=simplification,
        simplification_reactions=simplification_reactions)

    # Process
    process = determine_reaction_process_relevance(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        filtration_processes=filtration_processes)

    # Behavior
    behavior = determine_reaction_behavior_relevance(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        simplification=simplification,
        simplification_metabolites=simplification_metabolites)

    # Determine whether reaction is relevant
    relevance = (not simplification_match) and process and behavior

    return relevance

"""Determines whether a reaction has designation for simplification.
arguments:
    reaction_identifier (str): identifier of a reaction
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_reactions (list<dict<str>>): information about whether
        to simplify representations of specific reactions
returns:
    (bool): whether reaction has designation for simplification
"""
def determine_reaction_simplification(
        reaction_identifier,
        simplification,
        simplification_reactions):

    def match_reaction_simplification(
            record):

        return record['identifier'] == reaction_identifier

    # Determine whether to consider simplifications
    if simplification:
        match = find_match(
            match=match_reaction_simplification,
            sequence=simplification_reactions)

        if match is not None:
            omission = match['omission'] == 'True'

        else:
            omission = False

    else:
        omission = False

    return omission

"""Determines whether a reaction belongs to any relevant processes.
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    filtration_processes (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific processes
returns:
    (bool): whether reaction belongs to relevant process
"""
def determine_reaction_process_relevance(
        reaction_identifier,
        reactions,
        filtration_processes):

    reaction = reactions[reaction_identifier]
    reaction_processes = reaction["processes"]
    qualifiers = []

    for process in reaction_processes:

        relevance = determine_set_relevance(
            identifier=process,
            filters=filtration_processes)

        if relevance:
            qualifiers.append(process)

    return len(qualifiers) > 0

"""Determines whether a metabolic set, compartment or process, is relevant.
arguments:
    identifier (str): identifier of a set
    filters (list<dict<str>>): information about relevance of specific sets
returns:
    (bool): whether set is relevant
"""
def determine_set_relevance(
        identifier,
        filters):

    def match(
            record):

        return record['identifier'] == identifier

    record = find_match(
        match=match,
        sequence=filters)

    if record is not None:
        relevance = record['relevance'] == 'True'

        return relevance

    # If a set does not have a record in filters, then assume that it is relevant
    else:
        return True

"""Determines whether a reaction's behavior is relevant.
The relevance of a reaction's behavior depends on the relevance of
compartmentalization in general, and the relevance of specific compartments
and metabolites.
Factors in relevance of a reaction's behavior:
1. reaction's primary behavior, whether conversion or transport
2. compartmentalization
3. filters by compartments
4. simplification of metabolites
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (bool): whether reaction is a candidate
"""
def determine_reaction_behavior_relevance(
        reaction_identifier,
        reactions,
        compartmentalization,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    # Determine reaction's behavior.
    reaction = reactions[reaction_identifier]

    # Reaction involves chemical conversion.
    # Reaction's relevance requires a relevant reactant participant and
    # product participant that are chemically distinct.
    if reaction['conversion']:
        participation = determine_reaction_conversion_participation(
            reaction_identifier=reaction_identifier,
            reactions=reactions,
            filtration_compartments=filtration_compartments,
            simplification=simplification,
            simplification_metabolites=simplification_metabolites)

    # Reaction does not involve chemical conversion.
    # Reaction involves compartmental transport.
    # Reaction's relevance requires relevance of compartmentalization.
    elif reaction['transport']:

        # Reaction's relevance requires a relevant reactant participant and
        # product participant that are chemically identical and in separate
        # compartments
        if compartmentalization:
            participation = determine_reaction_transport_participation(
                reaction_identifier=reaction_identifier,
                reactions=reactions,
                filtration_compartments=filtration_compartments,
                simplification=simplification,
                simplification_metabolites=simplification_metabolites)

        else:
            participation = False

    return participation

"""Determines whether a reaction involves relevant participation for
conversion behavior.
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (bool): whether reaction involves relevant participation
"""
def determine_reaction_conversion_participation(
        reaction_identifier,
        reactions,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    # Determine relevant participants
    participants = determine_reaction_relevant_participants(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification=simplification,
        simplification_metabolites=simplification_metabolites)

    # Determine whether any reactant participants and product participants are relevant
    metabolites_reactant = collect_reaction_participants_value(
        key='metabolite',
        criteria={'roles': ['reactant']},
        participants=participants)
    metabolites_product = collect_reaction_participants_value(
        key='metabolite',
        criteria={'roles': ['product']},
        participants=participants)

    return (len(metabolites_reactant) > 0) and (len(metabolites_product) > 0)

"""Determines whether a reaction involves relevant participation for transport
behavior.
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (bool): whether reaction involves relevant participation
"""
def determine_reaction_transport_participation(
        reaction_identifier,
        reactions,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    # Determine relevant participants.
    participants = determine_reaction_relevant_participants(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification=simplification,
        simplification_metabolites=simplification_metabolites)

    # Determine whether any reactant participants and product participants match transport
    reaction = reactions[reaction_identifier]
    transports_original = reaction["transports"]
    transports_novel = []

    for transport in transports_original:

        metabolite = transport["metabolite"]
        compartments = transport["compartments"]
        compartments_reactant = collect_reaction_participants_value(
            key="compartment",
            criteria={
                "metabolites": [metabolite],
                "compartments": compartments,
                "roles": ["reactant"]},
            participants=participants)
        compartments_product = collect_reaction_participants_value(
            key="compartment",
            criteria={
                "metabolites": [metabolite],
                "compartments": compartments,
                "roles": ["product"]},
            participants=participants)

        if (len(compartments_reactant) > 0 and len(compartments_product) > 0):
            compartments_difference = (
                not utility.compare_lists_by_mutual_inclusion(
                    list_one=compartments_reactant,
                    list_two=compartments_product))

            if compartments_difference:
                transports_novel.append(transport)

    return len(transports_novel) > 0

"""Determines a reaction's relevant participants.
A participant's relevance depends on the relevance of its metabolite and
compartment.
Factors in relevance of a reaction's participant:
1. simplification of metabolites
2. filters by compartments
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (list<dict<str>>): reaction's relevant participants
"""
def determine_reaction_relevant_participants(
        reaction_identifier,
        reactions,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    reaction = reactions[reaction_identifier]
    participants_original = reaction["participants"]

    # Collect relevant participants
    participants_novel = []

    for participant in participants_original:

        # Determine relevance of participant's compartment
        compartment = participant["compartment"]
        compartment_relevance = determine_set_relevance(
            identifier=compartment,
            filters=filtration_compartments)

        # Determine relevance of participant's metabolite
        metabolite = participant["metabolite"]
        simplification_match = determine_metabolite_simplification(
            metabolite_identifier=metabolite,
            compartment_identifier=compartment,
            simplification=simplification,
            simplification_metabolites=simplification_metabolites)
        metabolite_relevance = (
            (not simplification_match["omission"]) and
            (not simplification_match["replication"]))

        # Determine whether participant is relevant
        if metabolite_relevance and compartment_relevance:
            participants_novel.append(participant)

    return participants_novel

"""Determines a metabolite's a designation for simplification.
A simplification must match both the metabolite and the compartment.
arguments:
    metabolite_identifier (str): identifier of a metabolite
    compartment_identifier (str): identifier of a compartment
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (dict<bool>): whether metabolite has designation for simplification by
        omission or replication
"""
def determine_metabolite_simplification(
        metabolite_identifier,
        compartment_identifier,
        simplification,
        simplification_metabolites):

    def match_metabolite(
            record):

        return record["metabolite"] == metabolite_identifier

    def match_compartment(
            record):

        return (
            (record["compartment"] == "all") or
            (record["compartment"] == compartment_identifier))

    # Determine whether to consider simplifications
    # There might be multiple records for a metabolite.
    if simplification:
        simplifications_metabolite = find_match(
            match=match_metabolite,
            sequence=simplification_metabolites)

        # Determine whether any of metabolite's simplifications match the compartment
        # There should be a single record for a metabolite in a specific compartment
        if simplifications_metabolite is not None:
            simplification_compartment = find_match(
                match=match_compartment,
                sequence=simplifications_metabolite)

            if simplification_compartment is not None:
                omission = simplification_compartment["omission"] == "True"
                replication = simplification_compartment["replication"] == "True"
                record = {
                    "omission": omission,
                    "replication": replication}

            else:
                record = {
                    "omission": False,
                    "replication": False}

        else:
            record = {
                "omission": False,
                "replication": False}

    else:
        record = {
            "omission": False,
            "replication": False}

    return record

"""Determines whether a reaction is redundant.
A reaction is redundant if its participants are identical to those of
another relevant reaction.
The relevance of compartmentalization determines whether compartments
influence this comparison of participants.
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    reactions_candidacy (dict<dict>): information about candidate reactions
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    filtration_processes (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific processes
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_reactions (list<dict<str>>): information about whether
        to simplify representations of specific reactions
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (tuple<bool, list<str>>): whether reaction is redundant, and relevant
        redundant replicate reactions
"""
def determine_reaction_redundancy(
        reaction_identifier,
        reactions,
        reactions_candidacy,
        compartmentalization,
        filtration_compartments,
        filtration_processes,
        simplification,
        simplification_reactions,
        simplification_metabolites):

    def match_relevance(
            replicate):

        identity = (replicate == reaction_identifier)
        relevance = determine_reaction_relevance(
            reaction_identifier=replicate,
            reactions=reactions,
            compartmentalization=compartmentalization,
            filtration_compartments=filtration_compartments,
            filtration_processes=filtration_processes,
            simplification=simplification,
            simplification_reactions=simplification_reactions,
            simplification_metabolites=simplification_metabolites)

        return (not identity) and relevance

    # Determine relevant replicates
    # Assume that all replicates correspond to records in reactions
    reaction = reactions[reaction_identifier]
    replicates = reaction["replicates"]
    replicates_relevant = list(filter(match_relevance, replicates))

    # Determine redundant replicates
    def match_redundancy(
            replicate):

        redundancy = determine_replicate_reactions_redundancy(
            reaction_one_identifier=reaction_identifier,
            reaction_two_identifier=replicate,
            reactions=reactions,
            compartmentalization=compartmentalization,
            filtration_compartments=filtration_compartments,
            simplification=simplification,
            simplification_metabolites=simplification_metabolites)

        return redundancy

    replicates_redundant = list(filter(match_redundancy, replicates_relevant))

    # Determine whether reaction is priority replicate.
    priority = determine_redundant_reaction_priority(
        reaction_identifier=reaction_identifier,
        reactions_replicates=replicates_redundant,
        reactions=reactions,
        reactions_candidacy=reactions_candidacy)
    redundancy = not priority

    return (
        redundancy,
        replicates_redundant)

"""Determines whether two replicate reactions are redundant.
arguments:
    reaction_one_identifier (str): identifier of a reaction
    reaction_two_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (bool): whether reactions are redundant
"""
def determine_replicate_reactions_redundancy(
        reaction_one_identifier,
        reaction_two_identifier,
        reactions,
        compartmentalization,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    reaction_one = reactions[reaction_one_identifier]
    reaction_two = reactions[reaction_two_identifier]

    # Determine whether reactions have redundant reversibility
    reversibility = reaction_one["reversibility"] == reaction_two["reversibility"]

    # Determine whether reactions have redundant participants
    # Determine relevant participants
    participation = determine_reactions_participation_redundancy(
        reaction_one_identifier=reaction_one_identifier,
        reaction_two_identifier=reaction_two_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        simplification=simplification,
        simplification_metabolites=simplification_metabolites)

    # Determine whether reactions are redundant
    return reversibility and participation

"""Determines whether two replicate reactions have redundant participation.
arguments:
    reaction_one_identifier (str): identifier of a reaction
    reaction_two_identifier (str): identifier of a reaction
    reactions (dict<dict>): information about reactions
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (bool): whether reactions' participants are redundant
"""
def determine_reactions_participation_redundancy(
        reaction_one_identifier,
        reaction_two_identifier,
        reactions,
        compartmentalization,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    # Determine relevant participants.
    participants_one = determine_reaction_relevant_participants(
        reaction_identifier=reaction_one_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification=simplification,
        simplification_metabolites=simplification_metabolites)
    participants_two = determine_reaction_relevant_participants(
        reaction_identifier=reaction_two_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification=simplification,
        simplification_metabolites=simplification_metabolites)

    # Determine whether participants are redundant.
    if compartmentalization:
        redundancy = determine_participants_attributes_mutual_redundancy(
            participants_one=participants_one,
            participants_two=participants_two,
            attributes=['metabolite', 'role', 'compartment'])

    else:
        redundancy = determine_participants_attributes_mutual_redundancy(
            participants_one=participants_one,
            participants_two=participants_two,
            attributes=['metabolite', 'role'])

    return redundancy

"""Determines whether participants of two reactions are redundant.
arguments:
    participants_one (str): participants of a reaction
    participants_two (str): participants of a reaction
    attributes (list<str>): names of attributes by which to compare
        participants
returns:
    (bool): whether reactions' participants are redundant
"""
def determine_participants_attributes_mutual_redundancy(
        participants_one,
        participants_two,
        attributes):

    comparison_one = determine_participants_attributes_redundancy(
        participants_one=participants_one,
        participants_two=participants_two,
        attributes=attributes)
    comparison_two = determine_participants_attributes_redundancy(
        participants_one=participants_two,
        participants_two=participants_one,
        attributes=attributes)

    return comparison_one and comparison_two

"""Determines whether participants of two reactions are redundant.
arguments:
    participants_one (str): participants of a reaction
    participants_two (str): participants of a reaction
    attributes (list<str>): names of attributes by which to compare
        participants
returns:
    (bool): whether reactions' participants are redundant
"""
def determine_participants_attributes_redundancy(
        participants_one,
        participants_two,
        attributes):

    iteration_one = []

    for participant_one in participants_one:

        iteration_two = []

        for participant_two in participants_two:

            iteration_three = []

            for attribute in attributes:

                match = participant_one[attribute] == participant_two[attribute]
                iteration_three.append(match)

            iteration_two.append(all(iteration_three))

        iteration_one.append(any(iteration_two))

    return all(iteration_one)

"""Determines whether a relevant, redundant reaction is a priority.
A reaction is a priority for candidacy if none of its other relevant,
redundant replicates are already candidates.
arguments:
    reaction_identifier (str): identifier of a reaction
    reactions_replicates (list<str>): identifiers of relevant, redundant
        reactions
    reactions (dict<dict>): information about reactions
    reactions_candidacy (dict<dict>): information about candidate reactions
returns:
    (bool): whether reaction is a priority
"""
def determine_redundant_reaction_priority(
        reaction_identifier,
        reactions_replicates,
        reactions,
        reactions_candidacy,):

    # Determine whether reaction has relevant, redundant replicates that are
    # already candidates.
    if len(reactions_replicates) > 0:
        candidates = []

        for replicate in reactions_replicates:

            candidate = replicate in reactions_candidacy
            candidates.append(candidate)

        priority = not any(candidates)

    else:
        priority = True

    return priority

# Candidate metabolites.
"""Collects information about candidate metabolites.
arguments:
    metabolites (dict<dict>): information about metabolites
    reactions (dict<dict>): information about reactions
    reactions_candidacy (dict<dict>): information about candidate reactions
    compartmentalization (bool): whether compartmentalization is relevant
    compartments (dict<dict>): information about compartments
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (dict<dict>): information about candidate reactions and candidate
        metabolites
"""
def collect_candidate_reactions_metabolites(
        metabolites,
        reactions,
        reactions_candidacy,
        compartmentalization,
        compartments,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    reactions_candidacy_metabolites = {}
    metabolites_candidacy = {}

    for reaction_candidacy in reactions_candidacy.values():

        candidates = collect_candidate_reaction_metabolites(
            reaction_candidacy_identifier=reaction_candidacy['identifier'],
            reactions=reactions,
            reactions_candidacy=reactions_candidacy,
            metabolites=metabolites,
            compartmentalization=compartmentalization,
            compartments=compartments,
            filtration_compartments=filtration_compartments,
            simplification=simplification,
            simplification_metabolites=simplification_metabolites)

        # Collect candidate reaction
        reactions_candidacy_metabolites[reaction_candidacy['identifier']] = candidates['reaction']

        # Determine whether candidate metabolites are novel
        for metabolite in candidates['metabolites']:

            if metabolite['identifier'] not in metabolites_candidacy.keys():
                metabolites_candidacy[metabolite['identifier']] = metabolite

    return {
        'reactions': reactions_candidacy_metabolites,
        'metabolites': metabolites_candidacy}

"""Collects information about a candidate reaction's candidate metabolites.
arguments:
    reaction_candidacy_identifier (str): identifier of a candidate reaction
    metabolites (dict<dict>): information about metabolites
    reactions (dict<dict>): information about reactions
    reactions_candidacy (dict<dict>): information about candidate reactions
    compartmentalization (bool): whether compartmentalization is relevant
    compartments (dict<dict>): information about compartments
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (dict): information about a candidate reaction and its candidate
        metabolites
"""
def collect_candidate_reaction_metabolites(
        reaction_candidacy_identifier,
        metabolites,
        reactions,
        reactions_candidacy,
        compartmentalization,
        compartments,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    # Access information
    reaction_candidacy = copy.deepcopy(reactions_candidacy[reaction_candidacy_identifier])
    reaction = reactions[reaction_candidacy['reaction']]
    participants = reaction['participants']

    # Determine candidate participants
    participants_candidacy = determine_reaction_candidate_participants(
        reaction_candidacy_identifier=reaction_candidacy_identifier,
        participants_original=participants,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        simplification=simplification,
        simplification_metabolites=simplification_metabolites)

    # Include information about candidate participants in information about
    # candidate reaction
    reaction_candidacy['participants'] = participants_candidacy

    # Determine candidate metabolites
    reaction_metabolites_candidacy = determine_reaction_candidate_metabolites(
        reaction_candidacy=reaction_candidacy,
        participants=participants_candidacy,
        metabolites=metabolites,
        compartmentalization=compartmentalization,
        compartments=compartments)

    # Include references in candidate reaction to candidate metabolites
    metabolites_candidacy_identifiers = collect_value_from_records(
        key='identifier',
        records=reaction_metabolites_candidacy)
    reaction_candidacy['metabolites_candidacy'] = collect_unique_elements(metabolites_candidacy_identifiers)
    reaction_candidacy['metabolites_candidacy_count'] = len(reaction_candidacy['metabolites_candidacy'])

    return {
        'reaction': reaction_candidacy,
        'metabolites': reaction_metabolites_candidacy}

"""Determines a reaction's candidate participants.
Candidate participants include metabolites with designations for
simplification by replication.
arguments:
    reaction_candidacy_identifier (str): identifier of a candidate reaction
    participants_original (list<dict<str>>): information about metabolites
        and compartments that participate in a reaction
    compartmentalization (bool): whether compartmentalization is relevant
    filtration_compartments (list<dict<str>>): information about whether to
        remove metabolites and reactions relevant to specific compartments
    simplification (bool): whether to simplify representations of specific
        entities in network
    simplification_metabolites (list<dict<str>>): information about whether
        to simplify representations of specific metabolites
returns:
    (list<dict<str>>): reaction's candidate participants
"""
def determine_reaction_candidate_participants(
        reaction_candidacy_identifier,
        participants_original,
        compartmentalization,
        filtration_compartments,
        simplification,
        simplification_metabolites):

    # Collect candidate participants
    participants_novel = []

    for participant in participants_original:

        # Determine relevance of participant's compartment
        compartment = participant['compartment']
        compartment_relevance = determine_set_relevance(
            identifier=compartment,
            filters=filtration_compartments)

        # Determine candidacy of participant's metabolite
        metabolite = participant['metabolite']
        simplification_match = determine_metabolite_simplification(
            metabolite_identifier=metabolite,
            compartment_identifier=compartment,
            simplification=simplification,
            simplification_metabolites=simplification_metabolites)
        participant['replication'] = simplification_match['replication']

        # Designate identifier for candidate metabolites for participant
        metabolite_candidacy = determine_candidate_metabolite_identifier(
            compartmentalization=compartmentalization,
            replication=simplification_match['replication'],
            metabolite_identifier=metabolite,
            compartment_identifier=compartment,
            reaction_candidacy_identifier=reaction_candidacy_identifier)
        participant['metabolite_candidacy'] = metabolite_candidacy

        # Determine whether participant is a candidate
        # Metabolites with designation for simplification by replication are
        # still candidates
        if compartment_relevance and not simplification_match['omission']:
            participants_novel.append(participant)

    return participants_novel

"""Determines a reaction's candidate metabolites.
arguments:
    reaction_candidacy (dict): information about a candidate reaction
    participants (list<dict<str>>): information about metabolites and
        compartments that participate in a reaction
    metabolites (dict<dict>): information about metabolites
    compartmentalization (bool): whether compartmentalization is relevant
    compartments (dict<dict>): information about compartments
returns:
    (list<dict<str>>): reaction's candidate metabolites
"""
def determine_reaction_candidate_metabolites(
        reaction_candidacy,
        participants,
        metabolites,
        compartmentalization,
        compartments):

    metabolites_candidacy = []

    for participant in participants:

        # Access information
        metabolite_identifier = participant['metabolite']
        compartment_identifier = participant['compartment']
        metabolite = metabolites[metabolite_identifier]
        compartment = compartments[compartment_identifier]
        metabolite_name = metabolite['name']
        compartment_name = compartment['name']
        replication=participant['replication']

        # Candidate metabolite identifier
        identifier = determine_candidate_metabolite_identifier(
            compartmentalization=compartmentalization,
            replication=replication,
            metabolite_identifier=metabolite_identifier,
            compartment_identifier=compartment_identifier,
            reaction_candidacy_identifier=reaction_candidacy['identifier'])

        # Candidate metabolite name
        name = determine_candidate_metabolite_name(
            compartmentalization=compartmentalization,
            replication=replication,
            metabolite_name=metabolite_name,
            compartment_name=compartment_name,
            reaction_candidacy_name=reaction_candidacy['name'])

        # Compartment
        if compartmentalization:
            compartment_reference = compartment_identifier

        else:
            compartment_reference = 'null'

        # Compile information
        record = {
            'identifier': identifier,
            'name': name,
            'metabolite': metabolite_identifier,
            'compartment': compartment_reference,
            'replication': replication}
        metabolites_candidacy.append(record)

    return metabolites_candidacy

"""Determines a candidate metabolite's identifier.
arguments:
    compartmentalization (bool): whether compartmentalization is relevant
    replication (bool): whether candidate metabolite has simplification by
        replication
    metabolite_identifier (str): identifier of a metabolite
    compartment_identifier (str): identifier of a compartment
    reaction_candidacy_identifier (str): identifier of a candidate reaction
returns:
    (str): identifier of a candidate metabolite
"""
def determine_candidate_metabolite_identifier(
        compartmentalization,
        replication,
        metabolite_identifier,
        compartment_identifier,
        reaction_candidacy_identifier):

    if compartmentalization:
        identifier_one = metabolite_identifier + '_' + compartment_identifier

    else:
        identifier_one = metabolite_identifier

    if replication:
        identifier_two = identifier_one + '_' + reaction_candidacy_identifier

    else:
        identifier_two = identifier_one

    return identifier_two

"""Determines a candidate metabolite's identifier.
arguments:
    compartmentalization (bool): whether compartmentalization is relevant
    metabolite_identifier (str): identifier of a metabolite
    compartment_identifier (str): identifier of a compartment
    reaction_candidacy_name (str): name of a candidate reaction
returns:
    (str): name of a candidate metabolite
"""
def determine_candidate_metabolite_name(
        compartmentalization,
        replication,
        metabolite_name,
        compartment_name,
        reaction_candidacy_name):

    if compartmentalization:
        name_one = metabolite_name + ' in ' + compartment_name

    else:
        name_one = metabolite_name

    if replication:
        name_two = name_one + ' for ' + reaction_candidacy_name

    else:
        name_two = name_one

    return name_two

"""Collects information about candidate metabolites.
arguments:
    metabolites_candidacy (dict<dict>): information about candidate
        metabolites
    reactions_candidacy (dict<dict>): information about candidate reactions
returns:
    (dict<dict>): information about candidate metabolites
"""
def collect_candidate_metabolites_reactions(
        metabolites_candidacy,
        reactions_candidacy):

    # Collect references to candidate reactions in which each candidate
    # metabolite participates
    metabolites_reactions = collect_records_targets_by_categories(
        target='identifier',
        category='metabolites_candidacy',
        records=list(reactions_candidacy.values()))
    metabolites_candidacy_reactions = {}

    for metabolite, reactions in metabolites_reactions.items():

        metabolite_candidacy = copy.deepcopy(metabolites_candidacy[metabolite])

        # Collect unique candidate reaction
        reactions_candidacy_unique = collect_unique_elements(reactions)
        metabolite_candidacy['reactions_candidacy'] = reactions_candidacy_unique
        metabolite_candidacy['reactions_candidacy_count'] = len(reactions_candidacy_unique)
        metabolites_candidacy_reactions[metabolite] = metabolite_candidacy

    return metabolites_candidacy_reactions

"""Converts information about reactions to text format.
arguments:
    reactions (dict<dict>): information about reactions
returns:
    (list<dict>): information about reactions
"""
def convert_reactions_text(
        reactions):

    records = []

    for reaction in reactions.values():

        # Compile information
        record = {
            'identifier': reaction['identifier'],
            'reaction': reaction['reaction'],
            'name': reaction['name'],
            'replicates': ';'.join(reaction['replicates']),
            'metabolites_candidacy': ';'.join(reaction['metabolites_candidacy']),
            'metabolites_candidacy_count': (reaction['metabolites_candidacy_count'])}
        records.append(record)

    records.sort(
        key=lambda record: record['metabolites_candidacy_count'],
        reverse=True)

    return records

"""Converts information about metabolites to text format.
arguments:
    metabolites (dict<dict>): information about metabolites
returns:
    (list<dict>): information about metabolites
"""
def convert_metabolites_text(
        metabolites):

    records = []

    for metabolite in metabolites.values():

        record = {
            'identifier': metabolite['identifier'],
            'name': metabolite['name'],
            'metabolite': metabolite['metabolite'],
            'compartment': metabolite['compartment'],
            'reactions_candidacy': ';'.join(metabolite['reactions_candidacy']),
            'reactions_candidacy_count': metabolite['reactions_candidacy_count']}
        records.append(record)

    records.sort(
        key=lambda record: record['reactions_candidacy_count'],
        reverse=True)

    return records

"""Writes product information to file
arguments:
    directory (str): directory for product files
    information (object): information to write to file
"""
def write_product(
        directory,
        information):

    # Specify directories and files
    path_reactions = directory + reaction_pickle
    path_metabolites = directory + metabolite_pickle
    path_metabolites_report = directory + metabolite_report
    path_reactions_report = directory + reaction_report

    # Write information to file
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

"""Function to execute module's main behavior.
The purpose of this procedure is to evaluate the candidacy of metabolites
and reactions for representation in a network.
arguments:
    compartmentalization (bool): whether compartmentalization is relevant
    simplification (bool): whether to simplify representations of specific
        entities in network
    directory (str): path to directory for source and product files
raises:
returns:
"""
def __main__(
        args_dict,
        compartmentalization=False,
        simplification=False):

    # Read source information from file
    print('Step 0/4: Reading in source data...')
    source = read_source(
        directory=args_dict['model'])

    # Collect candidate reactions
    print('Step 1/4: Collecting reactions...')
    reactions_candidacy = collect_candidate_reactions(
        reactions=source['reactions'],
        compartmentalization=compartmentalization,
        filtration_compartments=source['filtration_compartments'],
        filtration_processes=source['filtration_processes'],
        simplification=simplification,
        simplification_reactions=source['simplification_reactions'],
        simplification_metabolites=source['simplification_metabolites'])

    # Collect candidate metabolites
    # Include references to candidate metabolites with information about
    # candidate reactions
    print('Step 2/4: Collecting reactions...')
    reactions_metabolites_candidacy = collect_candidate_reactions_metabolites(
        metabolites=source['metabolites'],
        reactions=source['reactions'],
        reactions_candidacy=reactions_candidacy,
        compartmentalization=compartmentalization,
        compartments=source['compartments'],
        filtration_compartments=source['filtration_compartments'],
        simplification=simplification,
        simplification_metabolites=source['simplification_metabolites'])
    reactions_report = convert_reactions_text(
        reactions=reactions_metabolites_candidacy['reactions'])

    # Include references to candidate reactions with information about
    # candidate metabolites
    print('Step 3/4: Collecting metabolites...')
    metabolites_candidacy = collect_candidate_metabolites_reactions(
        reactions_candidacy=reactions_metabolites_candidacy['reactions'],
        metabolites_candidacy=reactions_metabolites_candidacy['metabolites'])
    metabolites_report = convert_metabolites_text(
        metabolites=metabolites_candidacy)

    # Compile information
    information = {
        'metabolites': metabolites_candidacy,
        'reactions': reactions_metabolites_candidacy['reactions'],
        'metabolites_report': metabolites_report,
        'reactions_report': reactions_report}

    #Write product information to file
    print('\nStep 4/4: Writing output...')
    write_product(
        directory=args_dict['candidates'],
        information=information)
