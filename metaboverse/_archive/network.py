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
import networkx

"""Import internal dependencies
"""
from metabalyze.metabonet_network.utils import collect_value_from_records
from metabalyze.metabonet_network.utils import collect_unique_elements
from metabalyze.metabonet_network.utils import collect_values_from_records_in_reference
from metabalyze.metabonet_network.utils import filter_entries_identifiers

"""Set globals
"""
compartment_pickle = 'compartments.pickle'
process_pickle = 'processes.pickle'
reaction_pickle = 'reactions.pickle'
metabolite_pickle = 'metabolites.pickle'

candidate_reaction_pickle = 'reactions.pickle'
candidate_metabolite_pickle = 'metabolites.pickle'

node_reactions_pickle = 'nodes_reactions.pickle'
node_metabolites_pickle = 'nodes_metabolites.pickle'
links_pickle = 'links.pickle'

network_pickle = 'network.pickle'

"""Reads and organizes source information from file
arguments:
    directory (str): directory of source files
returns:
    (object): source information
"""
def read_source(
        model_directory,
        candidates_directory):

    # Specify directories and files
    path_compartments = model_directory + compartment_pickle
    path_processes = model_directory + process_pickle
    path_reactions = model_directory + reaction_pickle
    path_metabolites = model_directory + metabolite_pickle

    path_reactions_candidacy = candidates_directory + candidate_reaction_pickle
    path_metabolites_candidacy = candidates_directory + candidate_metabolite_pickle

    # Read information from file
    with open(path_compartments, 'rb') as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, 'rb') as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, 'rb') as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, 'rb') as file_source:
        metabolites = pickle.load(file_source)

    with open(path_reactions_candidacy, 'rb') as file_source:
        reactions_candidacy = pickle.load(file_source)
    with open(path_metabolites_candidacy, 'rb') as file_source:
        metabolites_candidacy = pickle.load(file_source)

    return {
        'compartments': compartments,
        'processes': processes,
        'reactions': reactions,
        'metabolites': metabolites,
        'reactions_candidacy': reactions_candidacy,
        'metabolites_candidacy': metabolites_candidacy}

"""Collects information about reactions' nodes.
arguments:
    reactions_candidacy (dict<dict>): information about candidate reactions
    reactions (dict<dict>): information about reactions
    metabolites (dict<dict>): information about metabolites
    compartments (dict<dict>): information about compartments
    processes (dict<dict>): information about processes
returns:
    (dict<dict>): information about reactions' nodes
"""
def collect_reactions_nodes(
        reactions_candidacy,
        reactions,
        metabolites,
        compartments,
        processes):

    reactions_nodes = {}

    for reaction_candidacy in reactions_candidacy.values():

        reaction_node = define_reaction_node(
            reaction_candidacy_identifier=reaction_candidacy['identifier'],
            reactions_candidacy=reactions_candidacy,
            reactions=reactions,
            metabolites=metabolites,
            compartments=compartments,
            processes=processes)
        reactions_nodes[reaction_node['identifier']] = reaction_node

    return reactions_nodes

"""Defines information about a reaction's node.
arguments:
    reaction_candidacy_identifier (str): identifier of a candidate reaction
    reactions_candidacy (dict<dict>): information about candidate reactions
    reactions (dict<dict>): information about reactions
    metabolites (dict<dict>): information about metabolites
    compartments (dict<dict>): information about compartments
    processes (dict<dict>): information about processes
returns:
    (dict): information about a reaction's node
"""
def define_reaction_node(
        reaction_candidacy_identifier,
        reactions_candidacy,
        reactions,
        metabolites,
        compartments,
        processes):

    # Access information
    reaction_candidacy = reactions_candidacy[reaction_candidacy_identifier]
    reaction = reactions[reaction_candidacy['reaction']]

    # Compartments
    compartments_reaction = collect_value_from_records(
        key='compartment',
        records=reaction['participants'])
    compartments_unique = collect_unique_elements(
        elements_original=compartments_reaction)
    compartments_names = collect_values_from_records_in_reference(
        key='name',
        identifiers=compartments_unique,
        reference=compartments)

    # Processes
    processes_names = collect_values_from_records_in_reference(
        key='name',
        identifiers=reaction['processes'],
        reference=processes)

    # Metabolites
    metabolites_reaction = collect_value_from_records(
        key='metabolite',
        records=reaction['participants'])
    metabolites_unique = collect_unique_elements(
        elements_original=metabolites_reaction)
    metabolites_names = collect_values_from_records_in_reference(
        key='name',
        identifiers=metabolites_unique,
        reference=metabolites)

    # Compile information
    reaction_node = {
        'identifier': reaction_candidacy['identifier'],
        'type': 'reaction',
        'entity': reaction['identifier'],
        'name': reaction_candidacy['name'],
        'reversibility': reaction_candidacy['reversibility'],
        'metabolites': ';'.join(metabolites_names),
        'compartments': ';'.join(compartments_names),
        'processes': ';'.join(processes_names),
        'replicates': reaction_candidacy['replicates']}

    return reaction_node

"""Collects information about metabolites' nodes.
arguments:
    metabolites_candidacy (dict<dict>): information about candidate
        metabolites
    reactions (dict<dict>): information about reactions
    metabolites (dict<dict>): information about metabolites
    compartments (dict<dict>): information about compartments
    processes (dict<dict>): information about processes
returns:
    (dict<dict>): information about metabolites' nodes
"""
def collect_metabolites_nodes(
        metabolites_candidacy,
        reactions,
        metabolites,
        compartments,
        processes):

    metabolites_nodes = {}

    for metabolite_candidacy in metabolites_candidacy.values():

        metabolite_node = define_metabolite_node(
            metabolite_candidacy_identifier=metabolite_candidacy['identifier'],
            metabolites_candidacy=metabolites_candidacy,
            reactions=reactions,
            metabolites=metabolites,
            compartments=compartments,
            processes=processes)
        metabolites_nodes[metabolite_node['identifier']] = metabolite_node

    return metabolites_nodes

"""Defines information about a metabolite's node.
arguments:
    metabolite_candidacy_identifier (str): identifier of a candidate
        metabolite
    metabolites_candidacy (dict<dict>): information about candidate
        metabolites
    reactions (dict<dict>): information about reactions
    metabolites (dict<dict>): information about metabolites
    compartments (dict<dict>): information about compartments
    processes (dict<dict>): information about processes
returns:
    (dict): information about a metabolite's node
"""
def define_metabolite_node(
        metabolite_candidacy_identifier,
        metabolites_candidacy,
        reactions,
        metabolites,
        compartments,
        processes):

    # Access information
    metabolite_candidacy = metabolites_candidacy[metabolite_candidacy_identifier]
    metabolite = metabolites[metabolite_candidacy['metabolite']]
    compartment = metabolite_candidacy['compartment']

    if compartment != 'null':
        compartment_name = compartments[compartment]

    else:
        compartment_name = 'null'

    # Compile information
    metabolite_node = {
        'identifier': metabolite_candidacy['identifier'],
        'type': 'metabolite',
        'entity': metabolite['identifier'],
        'name': metabolite_candidacy['name'],
        'compartment': compartment_name,
        'formula': metabolite['formula'],
        'mass': metabolite['mass'],
        'charge': metabolite['charge'],
        'reference_hmdb': metabolite['references']['hmdb'],
        'reference_pubchem': metabolite['references']['pubchem'],
        'replication': metabolite_candidacy['replication']}

    return metabolite_node

"""Collects information about reactions' nodes.
arguments:
    reactions_candidacy (dict<dict>): information about candidate reactions
    metabolites_candidacy (dict<dict>): information about candidate
        metabolites
returns:
    (dict<dict>): information about links between reactions and metabolites
"""
def collect_links(
        reactions_candidacy,
        metabolites_candidacy):

    links = {}

    for reaction_candidacy in reactions_candidacy.values():

        links_reaction = define_reaction_links(
            reaction_candidacy_identifier=reaction_candidacy['identifier'],
            reactions_candidacy=reactions_candidacy,
            metabolites_candidacy=metabolites_candidacy)

        # Include links in collection
        for link in links_reaction:

            candidacy = determine_link_candidacy(
                link=link,
                links=links,
                reactions_candidacy=reactions_candidacy,
                metabolites_candidacy=metabolites_candidacy)

            if candidacy:
                links[link['identifier']] = link

    return links

"""Defines information about a reaction's links.
arguments:
    reaction_candidacy_identifier (str): identifier of a candidate reaction
    reactions_candidacy (dict<dict>): information about candidate reactions
    metabolites_candidacy (dict<dict>): information about candidate
        metabolites
returns:
    (list<dict>): information about reaction's links
"""
def define_reaction_links(
        reaction_candidacy_identifier,
        reactions_candidacy,
        metabolites_candidacy):

    # Access information
    reaction_candidacy = reactions_candidacy[reaction_candidacy_identifier]
    participants = reaction_candidacy['participants']

    # Create links for reaction's participants
    links_reaction = []

    for participant in participants:

        links_participant = define_reaction_participant_links(
            participant=participant,
            reversibility=reaction_candidacy['reversibility'],
            reaction_candidacy_identifier=reaction_candidacy_identifier)
        links_reaction.extend(links_participant)

    return links_reaction

"""Defines information about a reaction's participant's links.
arguments:
    participant (dict<str>): information about a metabolite and a
        compartment that participate in a reaction
    reversibility (bool): whether a reaction is reversible
    reaction_candidacy_identifier (str): identifier of a candidate reaction
returns:
    (list<dict>): information about reaction's participant's links
"""
def define_reaction_participant_links(
        participant,
        reversibility,
        reaction_candidacy_identifier):

    # Access information
    metabolite_candidacy_identifier = participant['metabolite_candidacy']
    replication = participant['replication']
    role = participant['role']

    # Determine which links to define
    if reversibility:
        # Reaction is reversible.
        # Define links for participant as both reactant and product.
        link_forward = define_link(
            source=metabolite_candidacy_identifier,
            target=reaction_candidacy_identifier,
            role=role,
            replication=replication)
        link_reverse = define_link(
            source=reaction_candidacy_identifier,
            target=metabolite_candidacy_identifier,
            role=role,
            replication=replication)
        links_participant = [link_forward, link_reverse]

    else:
        # Reaction is irreversible.
        # Define link for participant according to its role as reactant or
        # product.
        if role == 'reactant':
            link = define_link(
                source=metabolite_candidacy_identifier,
                target=reaction_candidacy_identifier,
                role=role,
                replication=replication)

        elif role == 'product':
            link = define_link(
                source=reaction_candidacy_identifier,
                target=metabolite_candidacy_identifier,
                role=role,
                replication=replication)

        links_participant = [link]

    return links_participant

"""Defines information about a link.
arguments:
    source (str): identifier of a node that is link's source
    target (str): identifier of a node that is link's target
    role (str): role of link's participant, either reactant or product
    replication (bool): whether candidate metabolite has simplification by
        replication
returns:
    (dict<str>): information about a link
"""
def define_link(
        source,
        target,
        role,
        replication):

    identifier = source + '_-_' + target
    link = {
        'identifier': identifier,
        'source': source,
        'target': target,
        'role': role,
        'replication': replication}

    return link

"""Determines whether link is a candidate.
arguments:
    link (dict<str>): information about a link
    links (dict<dict>): information about links between reactions and
        metabolites
    reactions_candidacy (dict<dict>): information about candidate reactions
    metabolites_candidacy (dict<dict>): information about candidate
        metabolites
returns:
    (bool): whether link is a candidate
"""
def determine_link_candidacy(
        link,
        links,
        reactions_candidacy,
        metabolites_candidacy):

    # Determine whether link refers to valid nodes
    source_reactions = link['source'] in reactions_candidacy.keys()
    source_metabolites = link['source'] in metabolites_candidacy.keys()
    source = source_reactions or source_metabolites
    target_reactions = link['target'] in reactions_candidacy.keys()
    target_metabolites = link['target'] in metabolites_candidacy.keys()
    target = target_reactions or target_metabolites

    # Determine whether link is novel in collection
    novelty = link['identifier'] not in links.keys()

    # Determine whether link is a candidate
    candidacy = source and target and novelty

    return candidacy

"""Determines whether to select the main component of a network.
arguments:
    component (bool): whether to select network's main component
    nodes_reactions (dict<dict>): information about reactions' nodes
    nodes_metabolites (dict<dict>): information about metabolites' nodes
    links (dict<dict>): information about links between nodes for reactions
        and metabolites
returns:
    (dict<dict<dict>>): information about network's nodes and links
"""
def select_network_component_elements(
        nodes_reactions,
        nodes_metabolites,
        links):

    # Convert information format for export to NetworkX
    networkx = conversion.convert_networkx(
        nodes_reactions=nodes_reactions,
        nodes_metabolites=nodes_metabolites,
        links=links)

    # Instantiate network in NetworkX
    network = analysis.instantiate_networkx(
        nodes=networkx['nodes'],
        links=networkx['links'])

    # Select network's main component
    components = networkx.algorithms.components.connected_components(network.to_undirected())
    components_sort = sorted(
        components,
        key=len,
        reverse=True)
    component_main = components_sort[0]

    # Induce subnetwork for component
    subnetwork = networkx.DiGraph.subgraph(
        network,
        list(component_main))

    # Extract identifiers of nodes and links from subnetwork
    nodes_identifiers = []

    for node, data in subnetwork.nodes.items():

        nodes_identifiers.append(data['identifier'])

    links_identifiers = []

    for link, data in subnetwork.edges.items():

        links_identifiers.append(data['identifier'])

    # Filter nodes and links by identifiers
    nodes_reactions_component = filter_entries_identifiers(
        identifiers=nodes_identifiers,
        entries_original=nodes_reactions)
    nodes_metabolites_component = filter_entries_identifiers(
        identifiers=nodes_identifiers,
        entries_original=nodes_metabolites)
    links_component = filter_entries_identifiers(
        identifiers=links_identifiers,
        entries_original=links)

    return {
        'nodes_reactions': nodes_reactions_component,
        'nodes_metabolites': nodes_metabolites_component,
        'links': links_component}

"""Writes product information to file
arguments:
    directory (str): directory for product files
    information (object): information to write to file
"""
def write_intermediates(
        directory,
        information):

    # Specify directories and files
    path_nodes_reactions = directory + node_reactions_pickle
    path_nodes_metabolites = directory + node_metabolites_pickle
    path_links = directory + links_pickle

    # Write information to file
    with open(path_nodes_reactions, 'wb') as file_product:
        pickle.dump(information['nodes_reactions'], file_product)
    with open(path_nodes_metabolites, 'wb') as file_product:
        pickle.dump(information['nodes_metabolites'], file_product)
    with open(path_links, 'wb') as file_product:
        pickle.dump(information['links'], file_product)

    return {
        'nodes_reactions': information['nodes_reactions'],
        'nodes_metabolites': information['nodes_metabolites'],
        'links': information['links']}

"""Converts information about network's nodes and links to format for
NetworkX.
Network is bipartite.
Store information about separate groups of nodes for reactions and
metabolites.
arguments:
    nodes_reactions (dict<dict>): information about reactions' nodes
    nodes_metabolites (dict<dict>): information about metabolites' nodes
    links (dict<dict>): information about links between nodes for reactions
        and metabolites
returns:
    (dict<list<tuple>>): information about network's nodes and links
"""
def convert_networkx(
        nodes_reactions,
        nodes_metabolites,
        links):

    nodes_networkx = convert_nodes_networkx(
        nodes_reactions=nodes_reactions,
        nodes_metabolites=nodes_metabolites)
    nodes_reactions_identifiers = collect_value_from_records(
        key='identifier',
        records=nodes_reactions.values())
    nodes_metabolites_identifiers = collect_value_from_records(
        key='identifier',
        records=nodes_metabolites.values())
    links_networkx = convert_links_networkx(links=links)

    return {
        'nodes': nodes_networkx,
        'nodes_reactions_identifiers': nodes_reactions_identifiers,
        'nodes_reactions': nodes_reactions,
        'nodes_metabolites_identifiers': nodes_metabolites_identifiers,
        'nodes_metabolites': nodes_metabolites,
        'links': links_networkx}

"""Converts information about network's nodes to format for NetworkX.
arguments:
    nodes_reactions (dict<dict>): information about network's nodes for
        reactions
    nodes_metabolites (dict<dict>): information about network's nodes for
        metabolites
returns:
    (list<tuple<str,dict>>): information about network's nodes
"""
def convert_nodes_networkx(
        nodes_reactions,
        nodes_metabolites):

    nodes_networkx = []

    for node_reaction in nodes_reactions.values():

        node_networkx = (node_reaction['identifier'], node_reaction)
        nodes_networkx.append(node_networkx)

    for node_metabolite in nodes_metabolites.values():

        node_networkx = (node_metabolite['identifier'], node_metabolite)
        nodes_networkx.append(node_networkx)

    return nodes_networkx

"""Converts information about network's links to format for NetworkX.
arguments:
    links (dict<dict>): information about links between nodes for reactions
        and metabolites
returns:
    (list<tuple<str,str,dict>>): information about network's links
"""
def convert_links_networkx(
        links):

    links_networkx = []

    for link in links.values():

        link_networkx = (link['source'], link['target'], link)
        links_networkx.append(link_networkx)

    return links_networkx

"""Writes product information to file
arguments:
    directory (str): directory for product files
    information (object): information to write to file
"""
def write_product(
        directory,
        network):

    # Specify directories and files
    path_networkx = directory + network_pickle

    # Write information to file.
    with open(path_networkx, 'wb') as file_product:
        pickle.dump(network, file_product)

"""The purpose of this procedure is to define network's elements and convert
information about network's elements to versatile formats, specifically for
compatibility with NetworkX.
arguments:
    directory (str): path to directory for source and product files
    component (bool): whether to select network's main component
"""
def __main__(
        args_dict):

    # Read source information from file
    print('Step 0/7: Reading in source data...')
    source = read_source(
        model_directory=args_dict['model'],
        candidates_directory=args_dict['candidates'])

    # Define network's nodes for reactions
    print('Step 1/7: Collecting reaction nodes...')
    nodes_reactions = collect_reactions_nodes(
        reactions_candidacy=source["reactions_candidacy"],
        reactions=source["reactions"],
        metabolites=source["metabolites"],
        compartments=source["compartments"],
        processes=source["processes"])

    # Define network's nodes for metabolites
    print('Step 2/7: Collecting metabolite nodes...')
    nodes_metabolites = collect_metabolites_nodes(
        metabolites_candidacy=source["metabolites_candidacy"],
        reactions=source["reactions"],
        metabolites=source["metabolites"],
        compartments=source["compartments"],
        processes=source["processes"])

    # Define network's links
    print('Step 3/7: Collecting node link information...')
    links = collect_links(
        reactions_candidacy=source["reactions_candidacy"],
        metabolites_candidacy=source["metabolites_candidacy"])

    # Determine whether to select network's main component
    print('Step 4/7: Selecting network components...')
    if 'component' in args_dict \
    and args_dict['component'] == True:
        component_elements = select_network_component_elements(
            nodes_reactions=nodes_reactions,
            nodes_metabolites=nodes_metabolites,
            links=links)
    else:
        print('--> Skipping as not specified...')
        component_elements = {
            'nodes_reactions': nodes_reactions,
            'nodes_metabolites': nodes_metabolites,
            'links': links}

    # Compile information
    print('Step 5/7: Saving intermediate network data...')
    intermediate_information = {
        "nodes_reactions": component_elements["nodes_reactions"],
        "nodes_metabolites": component_elements["nodes_metabolites"],
        "links": component_elements["links"]}

    #Write product information to file.
    output = write_intermediates(
        directory=args_dict['components'],
        information=intermediate_information)

    # Convert information format for export to NetworkX
    print('Step 6/7: Converting network for compatibility with NetworkX...')
    networkx = convert_networkx(
        nodes_reactions=output["nodes_reactions"],
        nodes_metabolites=output["nodes_metabolites"],
        links=output["links"])

    #Write product information to file.
    print('Step 7/7: Outputting network data...')
    write_product(
        directory=args_dict['network'],
        network=networkx)
