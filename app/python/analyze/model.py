"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
    Copyright (C) 2019  Jordan A. Berg
    jordan <dot> berg <at> biochem <dot> utah <dot> edu

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.
    You should have received a copy of the GNU General Public License along with
    this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""Import dependencies
"""
import os
from datetime import date
import pandas as pd
import numpy as np
from math import sqrt
from ast import literal_eval
import json
import pickle
import networkx as nx
from networkx.readwrite import json_graph
import matplotlib
import matplotlib.pyplot as plt
cmap = matplotlib.cm.get_cmap('seismic')
pmap = matplotlib.cm.get_cmap('Reds')

"""Import internal dependencies
"""
from analyze.collapse import collapse_nodes
from analyze.collapse import generate_updated_dictionary
from analyze.utils import convert_rgba
from utils import progress_feed

def test():

    args_dict = {'output': '/Users/jordan/Desktop/'}
    output_file = "/Users/jordan/Desktop/test.json"
    species_id = "HSA"
    network_url = "/Users/jordan/Desktop/HSA_metaboverse_db.pickle"
    with open(network_url, 'rb') as network_file:
        network = pickle.load(network_file)

"""Graph utils
"""
def name_graph(
        output_file,
        species_id):
    """Name graph
    """

    if output_file[-5:].lower() == '.json':
        graph_name = output_file
    else:
        graph_name = output_file + species_id + '_global_reactions.json'

    return graph_name

"""Graph building
"""
def build_graph(
        network,
        pathway_database,
        species_reference,
        name_reference,
        protein_reference,
        chebi_mapper,
        uniprot_reference,
        complexes,
        species_id,
        gene_reference,
        compartment_reference,
        component_database):
    """Build graph
    - Add nodes and edges
    - Map names to objects in the graph for display
    - Calculate degree for each node
    """

    # Initialize graph object
    G = nx.DiGraph()
    key_hash = set()
    remove_keys = []
    for reactome_id in network.keys():

        G, network, key_hash, remove_keys = process_reactions(
            graph=G,
            reactome_id=reactome_id,
            network=network,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            chebi_mapper=chebi_mapper,
            uniprot_reference=uniprot_reference,
            complex_reference=complexes,
            species_id=species_id,
            gene_reference=gene_reference,
            compartment_reference=compartment_reference,
            component_database=component_database,
            key_hash=key_hash,
            remove_keys=remove_keys)

    # Clean up duplicate reactions by ID
    for k in remove_keys:
        del network[k]

    for k, v in pathway_database.items():
        reactions = pathway_database[k]['reactions']
        updated_reactions = []
        for r in reactions:
            if r not in remove_keys:
                updated_reactions.append(r)
        pathway_database[k]['reactions'] = updated_reactions

    return G, network, pathway_database

def process_reactions(
        graph,
        reactome_id,
        network,
        species_reference,
        name_reference,
        protein_reference,
        chebi_mapper,
        uniprot_reference,
        complex_reference,
        species_id,
        gene_reference,
        compartment_reference,
        component_database,
        key_hash,
        remove_keys):
    """
    """
    new_components = []

    # Get reaction name and check if already addded
    reaction_name = network[reactome_id]['name']
    sort_name = ''.join(sorted(reaction_name.lower().replace(" ","")))

    if sort_name in key_hash:
        remove_keys.append(reactome_id)
    else:
        key_hash.add(sort_name)

        reaction_id = network[reactome_id]['id']
        reaction_rev = network[reactome_id]['reversible']
        reaction_notes = network[reactome_id]['notes']

        reactants = network[reactome_id]['reactants']
        products = network[reactome_id]['products']
        modifiers = network[reactome_id]['modifiers'] # ordered list
        compartment_id = network[reactome_id]['compartment']
        compartment_name = compartment_reference[compartment_id]

        prot_ref = {}
        for k, v in uniprot_reference.items():
            prot_ref[k] = v
            prot_ref[v] = k
        uniprot_reference = prot_ref

        flipped_ensembl = {}
        for k, v in gene_reference.items():
            flipped_ensembl[k] = v
            flipped_ensembl[v] = k

        # Add reaction node
        graph.add_node(reaction_id)
        graph.nodes()[reaction_id]['id'] = reactome_id
        graph.nodes()[reaction_id]['map_id'] = 'none'
        graph.nodes()[reaction_id]['name'] = reaction_name
        graph.nodes()[reaction_id]['reversible'] = reaction_rev
        graph.nodes()[reaction_id]['notes'] = reaction_notes
        graph.nodes()[reaction_id]['type'] = 'reaction'
        graph.nodes()[reaction_id]['sub_type'] = 'reaction'
        graph.nodes()[reaction_id]['compartment'] = compartment_id
        graph.nodes()[reaction_id]['compartment_display'] = compartment_name

        # Add vanilla element nodes and their edges
        for reactant in reactants:

            if component_database[reactant]['is'] != '':
                map_id = component_database[reactant]['is']
            else:
                map_id = 'none'

            graph = add_node_edge(
                graph=graph,
                id=reactant,
                map_id=map_id,
                name=component_database[reactant]['name'],
                compartment=component_database[reactant]['compartment'],
                reaction_membership=reaction_id,
                type='reactant',
                sub_type=component_database[reactant]['type'],
                reversible=reaction_rev,
                complex_reference=complex_reference,
                species_reference=species_reference,
                name_reference=name_reference,
                protein_reference=protein_reference,
                compartment_reference=compartment_reference)

            if len(component_database[reactant]['hasPart']) > 0:
                graph.nodes()[reactant]['complex'] = 'true'
                graph, additional_components = check_complexes(
                    species_id=species_id,
                    graph=graph,
                    id=reactant,
                    complex_reference=complex_reference,
                    species_reference=species_reference,
                    name_reference=name_reference,
                    protein_reference=protein_reference,
                    chebi_mapper=chebi_mapper,
                    uniprot_reference=uniprot_reference,
                    gene_reference=gene_reference,
                    component_database=component_database,
                    compartment_reference=compartment_reference)
                for x in additional_components:
                    new_components.append(x)
            else:
                graph.nodes()[reactant]['complex'] = 'false'

            if component_database[reactant]['type'] == 'protein_component':
                try:
                    uniprot_id = component_database[reactant]['is']
                    gene = flipped_ensembl[uniprot_reference[uniprot_id]]
                    new_components.append(gene)

                    graph = add_node_edge(
                        graph=graph,
                        id=gene,
                        map_id=gene,
                        name=gene_reference[gene],
                        compartment='none',
                        reaction_membership=reactant,
                        type='gene_component',
                        sub_type='gene',
                        reversible='false',
                        complex_reference=complex_reference,
                        species_reference=species_reference,
                        name_reference=name_reference,
                        protein_reference=protein_reference,
                        compartment_reference=compartment_reference)
                except:
                    pass

        for product in products:

            if component_database[product]['is'] != '':
                map_id = component_database[product]['is']
            else:
                map_id = 'none'

            graph = add_node_edge(
                graph=graph,
                id=product,
                map_id=map_id,
                name=component_database[product]['name'],
                compartment=component_database[product]['compartment'],
                reaction_membership=reaction_id,
                type='product',
                sub_type=component_database[product]['type'],
                reversible=reaction_rev,
                complex_reference=complex_reference,
                species_reference=species_reference,
                name_reference=name_reference,
                protein_reference=protein_reference,
                compartment_reference=compartment_reference)

            if len(component_database[product]['hasPart']) > 0:
                graph.nodes()[product]['complex'] = 'true'
                graph, additional_components = check_complexes(
                    species_id=species_id,
                    graph=graph,
                    id=product,
                    complex_reference=complex_reference,
                    species_reference=species_reference,
                    name_reference=name_reference,
                    protein_reference=protein_reference,
                    chebi_mapper=chebi_mapper,
                    uniprot_reference=uniprot_reference,
                    gene_reference=gene_reference,
                    component_database=component_database,
                    compartment_reference=compartment_reference)
                for x in additional_components:
                    new_components.append(x)
            else:
                graph.nodes()[product]['complex'] = 'false'

                if component_database[product]['type'] == 'protein_component':
                    try:
                        uniprot_id = component_database[product]['is']
                        gene = flipped_ensembl[uniprot_reference[uniprot_id]]
                        new_components.append(gene)

                        graph = add_node_edge(
                            graph=graph,
                            id=gene,
                            map_id=gene,
                            name=gene_reference[gene],
                            compartment='none',
                            reaction_membership=product,
                            type='gene_component',
                            sub_type='gene',
                            reversible='false',
                            complex_reference=complex_reference,
                            species_reference=species_reference,
                            name_reference=name_reference,
                            protein_reference=protein_reference,
                            compartment_reference=compartment_reference)
                    except:
                        pass

        for modifier in modifiers:

            # Extract modifier type
            # Formatted as a list of lists with first index of each sub list being
            # the species ID and the second index being the modifier type
            # ex: [['species_0', 'inhibitor'], ['species_1', 'catalyst']]
            # Labeling the edge should allow for differentiation between the same
            # modifier node acting as a catalyst or inhibitor
            id = modifier[0]
            type = modifier[1]

            if component_database[id]['is'] != '':
                map_id = component_database[id]['is']
            else:
                map_id = 'none'

            graph = add_node_edge(
                graph=graph,
                id=id,
                map_id=map_id,
                name=component_database[id]['name'],
                reaction_membership=reaction_id,
                compartment=component_database[id]['compartment'],
                type=type,
                sub_type=component_database[id]['type'],
                reversible='false',
                complex_reference=complex_reference,
                species_reference=species_reference,
                name_reference=name_reference,
                protein_reference=protein_reference,
                compartment_reference=compartment_reference)

            if len(component_database[id]['hasPart']) > 0:
                graph.nodes()[id]['complex'] = 'true'
                graph, additional_components = check_complexes(
                    species_id=species_id,
                    graph=graph,
                    id=id,
                    complex_reference=complex_reference,
                    species_reference=species_reference,
                    name_reference=name_reference,
                    protein_reference=protein_reference,
                    chebi_mapper=chebi_mapper,
                    uniprot_reference=uniprot_reference,
                    gene_reference=gene_reference,
                    component_database=component_database,
                    compartment_reference=compartment_reference)
                for x in additional_components:
                    new_components.append(x)
            else:
                graph.nodes()[id]['complex'] = 'false'

                if component_database[id]['type'] == 'protein_component':
                    try:
                        uniprot_id = component_database[id]['is']
                        gene = flipped_ensembl[uniprot_reference[uniprot_id]]
                        new_components.append(gene)
                        graph = add_node_edge(
                            graph=graph,
                            id=gene,
                            map_id=gene,
                            name=gene_reference[gene],
                            compartment='none',
                            reaction_membership=id,
                            type='gene_component',
                            sub_type='gene',
                            reversible='false',
                            complex_reference=complex_reference,
                            species_reference=species_reference,
                            name_reference=name_reference,
                            protein_reference=protein_reference,
                            compartment_reference=compartment_reference)
                    except:
                        pass

        network[reactome_id]['additional_components'] = new_components

    return graph, network, key_hash, remove_keys

def add_node_edge(
        graph,
        id, # node id
        map_id, # used for mapping user data
        name,
        compartment,
        reaction_membership,
        type,
        sub_type,
        reversible,
        complex_reference,
        species_reference,
        name_reference,
        protein_reference,
        compartment_reference):
    """Add node and edge information to graph
    """

    graph.add_node(id)
    graph.nodes()[id]['id'] = id
    graph.nodes()[id]['map_id'] = map_id
    graph.nodes()[id]['name'] = name
    graph.nodes()[id]['type'] = type
    graph.nodes()[id]['sub_type'] = sub_type
    graph.nodes()[id]['inferred'] = 'false'
    try:
        graph.nodes()[id]['compartment'] = compartment
        graph.nodes()[id]['compartment_name'] = compartment_reference[compartment]
    except:
        graph.nodes()[id]['compartment'] = 'none'
        graph.nodes()[id]['compartment_name'] = 'none'

    if type == 'reactant':
        graph.add_edges_from([
            (id, reaction_membership)])
        graph.edges()[(id, reaction_membership)]['type'] = type
        graph.edges()[(id, reaction_membership)]['sub_type'] = sub_type

        if reversible == 'true':
            graph.add_edges_from([
                (reaction_membership, id)])
            graph.edges()[(reaction_membership, id)]['type'] = type
            graph.edges()[(reaction_membership, id)]['sub_type'] = sub_type

    elif type == 'product':
        graph.add_edges_from([
            (reaction_membership, id)])
        graph.edges()[(reaction_membership, id)]['type'] = type
        graph.edges()[(reaction_membership, id)]['sub_type'] = sub_type

        if reversible == 'true':
            graph.add_edges_from([
                (id, reaction_membership)])
            graph.edges()[(id, reaction_membership)]['type'] = type
            graph.edges()[(id, reaction_membership)]['sub_type'] = sub_type

    else:
        graph.add_edges_from([
            (id, reaction_membership)])
        graph.edges()[(id, reaction_membership)]['type'] = type
        graph.edges()[(id, reaction_membership)]['sub_type'] = sub_type

    return graph

def check_complexes(
        species_id,
        graph,
        id,
        complex_reference,
        species_reference,
        name_reference,
        protein_reference,
        chebi_mapper,
        uniprot_reference,
        gene_reference,
        component_database,
        compartment_reference):
    """Check if species being added is in complex dictionary
    - If record exists, add nodes and edges for the new relationship.
    - If record contains a UniProt ID, cross reference with Ensembl database
    - If complex, label true; else label as false
    """

    add_components = []

    for x in component_database[id]['hasPart']:

        if x in uniprot_reference:
            name = uniprot_reference[x]
            type = 'complex_component'
            sub_type = 'protein_component'
            add_components.append(x)
            display_name = uniprot_reference[x]

            graph = add_node_edge(
                graph=graph,
                id=x,
                map_id=x,
                name=display_name,
                compartment='none',
                reaction_membership=id,
                type='complex_component',
                sub_type=sub_type,
                reversible='false',
                complex_reference=complex_reference,
                species_reference=species_reference,
                name_reference=name_reference,
                protein_reference=protein_reference,
                compartment_reference=compartment_reference)

            try:
                if x in protein_reference.keys():
                    gene = protein_reference[x]
                    name = gene_reference[gene]
                else:
                    gene = name = display_name
                add_components.append(gene)

                graph = add_node_edge(
                    graph=graph,
                    id=gene,
                    map_id=gene,
                    name=name,
                    compartment='none',
                    reaction_membership=x,
                    type='gene_component',
                    sub_type='gene',
                    reversible='false',
                    complex_reference=complex_reference,
                    species_reference=species_reference,
                    name_reference=name_reference,
                    protein_reference=protein_reference,
                    compartment_reference=compartment_reference)
            except:
                print('---------------')
                print(x)
                print(name)
                print(display_name)
                print('=================')

        else:
            if 'chebi' in x.lower():
                map_id = name = x
                sub_type = 'metabolite_component'

            elif x.lower in chebi_mapper.keys():
                name = x
                map_id = chebi_mapper[x]
                sub_type = 'metabolite_component'

            elif 'mi' in x.lower():
                map_id = name = x
                sub_type = 'mirna_component'

            else:
                map_id = name = x
                sub_type = 'other'

            if x.lower in chebi_mapper.keys():
                component_id = map_id
                display_name = x
            else:
                try:
                    component_id = name_reference[name]
                    display_name = species_reference[component_id]
                except:
                    display_name = component_id = x

            add_components.append(component_id)

            graph = add_node_edge(
                graph=graph,
                id=component_id,
                map_id=map_id,
                name=display_name,
                compartment='none',
                reaction_membership=id,
                type='complex_component',
                sub_type=sub_type,
                reversible='false',
                complex_reference=complex_reference,
                species_reference=species_reference,
                name_reference=name_reference,
                protein_reference=protein_reference,
                compartment_reference=compartment_reference)

    return graph, add_components

def uniprot_ensembl_reference(
        uniprot_reference,
        ensembl_reference):
    """Build cross-referencing dictionary to convert uniprot ID to
    corresponding Ensembl ID
    """

    new_dict = {}
    for k, v in uniprot_reference.items():
        try:
            new_dict[k] = ensembl_reference[v]
        except:
            pass

    return new_dict

def parse_attributes(
        dataframe,
        name_reference,
        chebi_mapper,
        metabolite_mapper,
        ignore_enantiomers=True):

    dataframe_dict = {}
    non_mappers = []
    for index, row in dataframe.iterrows():
        dataframe_dict[index] = {
            'values': list(row),
            'mapping_synonyms': [],
            'display_synonyms': [],
            'chebi_ids': set()
        }

        if index in name_reference.keys():
            dataframe_dict[index]['mapping_synonyms'] = [index, name_reference[index]]
            dataframe_dict[index]['display_synonyms'] = [index, name_reference[index]]

            dataframe_dict[name_reference[index]] = {
                'values': list(row),
                'mapping_synonyms': [],
                'display_synonyms': [],
                'chebi_ids': set()
            }
            dataframe_dict[name_reference[index]]['mapping_synonyms'] = [index, name_reference[index]]
            dataframe_dict[name_reference[index]]['display_synonyms'] = [index, name_reference[index]]

        else:
            _i = ''.join(c.lower() for c in str(index) if c.isalnum())

            try:
                __i = metabolite_mapper['mapping_dictionary'][_i]
                dataframe_dict[index]['mapping_synonyms'] = \
                    metabolite_mapper['hmdb_dictionary'][__i]
                dataframe_dict[index]['display_synonyms'] = \
                    metabolite_mapper['display_dictionary'][__i]
            except:

                try:
                    # remove last letter
                    __i = metabolite_mapper['mapping_dictionary'][_i[:-1]]
                    dataframe_dict[index]['mapping_synonyms'] = \
                        metabolite_mapper['hmdb_dictionary'][__i]
                    dataframe_dict[index]['display_synonyms'] = \
                        metabolite_mapper['display_dictionary'][__i]
                except:
                    if _i[0] == 'd' or _i[0] == 'l' and ignore_enantiomers == True:
                        try:
                            __i = metabolite_mapper['mapping_dictionary'][_i[1:]]
                            dataframe_dict[index]['mapping_synonyms'] = \
                                metabolite_mapper['hmdb_dictionary'][__i]
                            dataframe_dict[index]['display_synonyms'] = \
                                metabolite_mapper['display_dictionary'][__i]
                        except:
                            non_mappers.append(index)

                    else:
                        non_mappers.append(index)

    dataframe_mapper = {}
    chebi_synonyms = {}
    for k, v in dataframe_dict.items():

        if k in chebi_mapper:
            dataframe_dict[k]['chebi_ids'].add(chebi_mapper[k])

        for s in dataframe_dict[k]['mapping_synonyms']:
            if s in chebi_mapper:
                dataframe_dict[k]['chebi_ids'].add(chebi_mapper[s])
            else:
                dataframe_dict[k]['chebi_ids'].add(s)
        for id in list(dataframe_dict[k]['chebi_ids']):
            dataframe_mapper[id] = dataframe_dict[k]['values']
            chebi_synonyms[id] = dataframe_dict[k]['display_synonyms']

    return dataframe_mapper, chebi_synonyms, non_mappers

def map_attributes(
        graph,
        data,
        stats,
        name_reference,
        degree_dictionary,
        chebi_mapper,
        metabolite_mapper):
    """Data overlay
    - Map repo id to species_id
    - If a node is a complex, take average of neighbors that are not
    To do:
    - Currently, many metabolites that should map are not found in name
    database
    """

    n = len(data.columns.tolist())
    reaction_color = (0.75, 0.75, 0.75, 1)
    missing_color = (1, 1, 1, 1)

    # Re-index data and stats
    data_renamed = data.copy()
    #data_renamed = data.rename(index=name_reference)
    data_renamed = data_renamed.loc[data_renamed.index.dropna()]
    data_max = abs(data_renamed).max().max()

    stats_renamed = stats.copy()
    #stats_renamed = stats.rename(index=name_reference)
    stats_renamed = stats_renamed.loc[stats_renamed.index.dropna()]
    stats_logged = -1 * np.log10(stats_renamed + 1e-100)

    stats_max = abs(stats_logged).max().max()

    data_dict, chebi_synonyms, non_mappers1 = parse_attributes(
        dataframe=data_renamed,
        name_reference=name_reference,
        chebi_mapper=chebi_mapper,
        metabolite_mapper=metabolite_mapper,
        ignore_enantiomers=True)

    stats_dict, chebi_synonyms, non_mappers2 = parse_attributes(
        dataframe=stats_renamed,
        name_reference=name_reference,
        chebi_mapper=chebi_mapper,
        metabolite_mapper=metabolite_mapper,
        ignore_enantiomers=True)

    non_mappers = [set(non_mappers1 + non_mappers2)]

    data_renamed = None
    stats_renamed = None
    chebi_mapper = None
    metabolite_mapper = None

    for current_id in list(graph.nodes()):

        x = current_id
        map_id = graph.nodes()[current_id]['map_id']
        backup_mapper = graph.nodes()[current_id]['name']

        # Add degree
        graph.nodes()[x]['degree'] = degree_dictionary[current_id]
        if map_id != 'none' and map_id in chebi_synonyms:
            graph.nodes()[x]['synonyms'] = chebi_synonyms[map_id]
        else:
            graph.nodes()[x]['synonyms'] = []

        if graph.nodes()[x]['type'] == 'reaction':
            colors = [reaction_color for x in range(n)]

            graph.nodes()[x]['values'] = [None for x in range(n)]
            graph.nodes()[x]['values_rgba'] = colors
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=colors)

            graph.nodes()[x]['stats'] = [None for x in range(n)]
            graph.nodes()[x]['stats_rgba'] = colors
            graph.nodes()[x]['stats_js'] = convert_rgba(
                rgba_tuples=colors)

        elif map_id in set(data_dict.keys()) \
        and map_id in set(stats_dict.keys()) \
        and map_id != 'none':
            graph.nodes()[x]['values'] = data_dict[map_id]
            graph.nodes()[x]['values_rgba'] = extract_value(
                value_array=data_dict[map_id],
                max_value=data_max)
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['values_rgba'])

            graph.nodes()[x]['stats'] = stats_dict[map_id]
            graph.nodes()[x]['stats_rgba'] = extract_value(
                value_array=stats_dict[map_id],
                max_value=stats_max,
                type="stats")
            graph.nodes()[x]['stats_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['stats_rgba'])

        elif backup_mapper in set(data_dict.keys()) \
        and backup_mapper in set(stats_dict.keys()) \
        and backup_mapper != 'none':
            graph.nodes()[x]['values'] = data_dict[backup_mapper]
            graph.nodes()[x]['values_rgba'] = extract_value(
                value_array=data_dict[backup_mapper],
                max_value=data_max)
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['values_rgba'])

            graph.nodes()[x]['stats'] = stats_dict[backup_mapper]
            graph.nodes()[x]['stats_rgba'] = extract_value(
                value_array=stats_dict[backup_mapper],
                max_value=stats_max,
                type="stats")
            graph.nodes()[x]['stats_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['stats_rgba'])

        else:
            colors = [missing_color for x in range(n)]

            graph.nodes()[x]['values'] = [None for x in range(n)]
            graph.nodes()[x]['values_rgba'] = colors
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=colors)

            graph.nodes()[x]['stats'] = [None for x in range(n)]
            graph.nodes()[x]['stats_rgba'] = colors
            graph.nodes()[x]['stats_js'] = convert_rgba(
                rgba_tuples=colors)

    return graph, data_max, stats_max, non_mappers

def extract_value(
        value_array,
        max_value,
        type="value"):
    """Extract expression value
    """

    rgba = []

    if type == "value":

        for x in value_array:

            position = (x + max_value) / (2 * max_value)
            rgba_tuple = cmap(position)
            rgba.append(rgba_tuple)

    else:

        for x in value_array:

            x = -1 * np.log10(x + 1e-100)
            position = x / max_value
            rgba_tuple = pmap(position)
            rgba.append(rgba_tuple)

    return rgba

def output_graph(
        graph,
        output_name,
        pathway_dictionary,
        collapsed_pathway_dictionary,
        super_pathways,
        reaction_dictionary,
        collapsed_reaction_dictionary,
        motif_reaction_dictionary,
        mod_collapsed_pathways,
        degree_dictionary,
        max_value,
        max_stat,
        categories,
        labels,
        blacklist,
        metadata,
        unmapped):
    """Output graph and necessary metadata
    """

    data = json_graph.node_link_data(graph)
    data['pathway_dictionary'] = pathway_dictionary
    data['collapsed_pathway_dictionary'] = collapsed_pathway_dictionary
    data['super_pathways'] = super_pathways
    data['reaction_dictionary'] = reaction_dictionary
    data['collapsed_reaction_dictionary'] = collapsed_reaction_dictionary
    data['motif_reaction_dictionary'] = motif_reaction_dictionary
    data['mod_collapsed_pathways'] = mod_collapsed_pathways
    data['degree_dictionary'] = degree_dictionary
    data['max_value'] = max_value
    data['max_stat'] = max_stat
    data['categories'] = categories
    data['labels'] = labels
    data['blacklist'] = blacklist
    data['metadata'] = metadata
    data['unmapped'] = unmapped

    with open(output_name, 'w') as f:
        json.dump(data, f, indent=4) # Parse out as array for javascript

def compile_pathway_degree(
        pathways,
        scale_factor=200):
    """Compile database of large pathways
    """

    super_pathways = {}

    for key in list(pathways.keys()):

        if len(pathways[key]["reactions"]) > scale_factor:
            super_pathways[key] = pathways[key]

    return super_pathways

def compile_node_degrees(
        graph):
    """Retrieve degree metrics for each node
    """

    d = {}

    deg_dict = graph.degree
    for k,v in deg_dict:

        d[k] = v

    return d

def remove_nulls(values):

    if all(None in v for v in values) == True:
        return []

    else:
        for v in values:
            if None in v:
                values.remove(v)

        return values

def infer_protein_values(values, length):

    protein_vals = []
    for i in range(length):

        pos = []

        for j in range(len(values)):

            if values[j][i] != None:
                pos.append(values[j][i])

        protein_vals.append(sum(pos) / len(pos))

    return protein_vals

def infer_protein_stats(stats, length):

    protein_stats = []
    for i in range(length):

        pos = []

        for j in range(len(stats)):

            if stats[j][i] != None:
                pos.append(stats[j][i])

        protein_stats.append(max(pos))

    return protein_stats

def broadcast_values(
        graph,
        categories,
        max_value,
        max_stat,
        broadcast_genes=True):
    """
    """

    u = graph.to_undirected()
    length = len(categories)

    if broadcast_genes == True:
        for x in graph.nodes():

            if None not in graph.nodes()[x]['values'] \
            and None not in graph.nodes()[x]['stats']:
                pass

            elif graph.nodes()[x]['type'] == 'reaction':
                pass

            else:

                # 1. sub_type == 'protein_component' && type == 'complex_component'
                #       find edges
                #       find nodes type == 'gene_component'
                #       for x in values:
                #           take min, max, avg of genes to broadcast
                #           mark as inferred

                if graph.nodes()[x]['sub_type'] == 'protein_component' \
                and graph.nodes()[x]['type'] == 'complex_component':

                    gene_values = []
                    gene_stats = []
                    for neighbor in u[x]:

                        if graph.nodes()[neighbor]['sub_type'] == 'gene':
                            gene_values.append(graph.nodes()[neighbor]['values'])
                            gene_stats.append(graph.nodes()[neighbor]['stats'])

                    # Remove None and avg
                    gene_values = remove_nulls(gene_values)
                    gene_stats = remove_nulls(gene_stats)

                    # Mark as inferred if this is possible
                    if gene_values != []:
                        inferred_values = infer_protein_values(gene_values, length)

                        graph.nodes()[x]['inferred'] = 'true'
                        graph.nodes()[x]['values'] = inferred_values
                        graph.nodes()[x]['values_rgba'] = extract_value(
                            value_array=inferred_values,
                            max_value=max_value)
                        graph.nodes()[x]['values_js'] = convert_rgba(
                            rgba_tuples=graph.nodes()[x]['values_rgba'])

                    if gene_stats != []:
                        inferred_stats = infer_protein_stats(gene_stats, length)

                        graph.nodes()[x]['inferred'] = 'true'
                        graph.nodes()[x]['stats'] = inferred_stats
                        graph.nodes()[x]['stats_rgba'] = extract_value(
                            value_array=inferred_stats,
                            max_value=max_stat,
                            type="stats")
                        graph.nodes()[x]['stats_js'] = convert_rgba(
                            rgba_tuples=graph.nodes()[x]['stats_rgba'])

    for x in graph.nodes():

        if None not in graph.nodes()[x]['values'] \
        and None not in graph.nodes()[x]['stats']:
            pass

        elif graph.nodes()[x]['type'] == 'reaction':
            pass

        else:

            # 2. graph.nodes()[id]['complex'] == 'true'
            #       find edges
            #       find nodes type == 'complex_component'
            #       for x in values:
            #           take min, max, avg of genes to broadcast
            #           mark as inferred

            if 'complex' in graph.nodes()[x] \
            and graph.nodes()[x]['complex'] == 'true':

                gene_values = []
                gene_stats = []
                for neighbor in u[x]:

                    if graph.nodes()[neighbor]['type'] == 'complex_component':
                        gene_values.append(graph.nodes()[neighbor]['values'])
                        gene_stats.append(graph.nodes()[neighbor]['stats'])

                # Remove None and avg
                gene_values = remove_nulls(gene_values)
                gene_stats = remove_nulls(gene_stats)

                # Mark as inferred if this is possible
                if gene_values != []:
                    inferred_values = infer_protein_values(gene_values, length)

                    graph.nodes()[x]['inferred'] = 'true'
                    graph.nodes()[x]['values'] = inferred_values
                    graph.nodes()[x]['values_rgba'] = extract_value(
                        value_array=inferred_values,
                        max_value=max_value)
                    graph.nodes()[x]['values_js'] = convert_rgba(
                        rgba_tuples=graph.nodes()[x]['values_rgba'])

                if gene_stats != []:
                    inferred_stats = infer_protein_stats(gene_stats, length)

                    graph.nodes()[x]['inferred'] = 'true'
                    graph.nodes()[x]['stats'] = inferred_stats
                    graph.nodes()[x]['stats_rgba'] = extract_value(
                        value_array=inferred_stats,
                        max_value=max_stat,
                        type="stats")
                    graph.nodes()[x]['stats_js'] = convert_rgba(
                        rgba_tuples=graph.nodes()[x]['stats_rgba'])

    return graph

def make_motif_reaction_dictionary(
        network,
        updated_reactions,
        updated_pathway_dictionary):

    simp_dict = {}
    for k,v in network['reaction_database'].items():
        simp_dict[k] = []

    for k,v in network['pathway_database'].items():

        for x in network['pathway_database'][k]['reactions']:
            try:
                simp_dict[x].append(network['pathway_database'][k]['id'])
            except:
                simp_dict[x] = []
                simp_dict[x].append(network['pathway_database'][k]['id'])

    motif_reaction_dictionary = {}
    for k,v in updated_reactions.items():
        motif_reaction_dictionary[k] = []

    for k,v in updated_pathway_dictionary.items():

        for x in updated_pathway_dictionary[k]['reactions']:
            try:
                motif_reaction_dictionary[x].append(updated_pathway_dictionary[k]['id'])
            except:
                motif_reaction_dictionary[x] = []
                motif_reaction_dictionary[x].append(updated_pathway_dictionary[k]['id'])

    for k,v in motif_reaction_dictionary.items():
        if len(v) == 0:
            comps = k.split('_reaction_')
            comps = [x.replace('reaction_','') for x in comps]

            for y in comps:
                _rxn = 'reaction_' + y
                _paths = simp_dict[_rxn]
                for z in _paths:
                    motif_reaction_dictionary[k].append(z)

    return motif_reaction_dictionary

def make_metabolite_synonym_dictionary(
        network,
        output_dir,
        ref_url='https://sourceforge.net/projects/metaboverse/files/utils/metabolite_mapping.pickle.zip/download'):

    output_dir = '/Users/jordan/Desktop/'

    print("Downloading metabolite mapper...")
    os.system('curl -L ' + ref_url + ' -o ' + output_dir + 'metabolite_mapping.pickle.zip')
    print("Unzipping metabolite mapper...")
    os.system('unzip ' + output_dir + 'metabolite_mapping.pickle.zip -d ' + output_dir)
    print("Parsing HMDB metabolite records...")
    with open(output_dir + 'metabolite_mapping.pickle', 'rb') as ref_file:
        metabolite_mapper = pickle.load(ref_file)
    os.remove(output_dir + 'metabolite_mapping.pickle.zip')
    os.remove(output_dir + 'metabolite_mapping.pickle')

    return metabolite_mapper

def __main__(
        args_dict,
        network,
        data,
        stats,
        species_id,
        output_file,
        unmapped,
        flag_data=False):
    """Generate graph object for visualization
    """

    print('Preparing metadata...')
    # Generate output file name
    graph_name = name_graph(
        output_file=output_file,
        species_id=species_id)

    # Prepare uniprot to ensembl name mapper
    reverse_genes = {v:k for k,v in network['ensembl_synonyms'].items()}
    protein_dictionary = uniprot_ensembl_reference(
        uniprot_reference=network['uniprot_synonyms'],
        ensembl_reference=reverse_genes)
    progress_feed(args_dict, "model", 1)

    chebi_mapper = {}
    for k, v in network['chebi_synonyms'].items():
        _k = ''.join(c.lower() for c in str(k) if c.isalnum())
        chebi_mapper[_k] = v
        chebi_mapper[k] = v
    for k, v in network['uniprot_metabolites'].items():
        _k = ''.join(c.lower() for c in str(k) if c.isalnum())
        chebi_mapper[_k] = v
        chebi_mapper[k] = v

    # Generate graph
    # Name mapping
    print('Building network...')
    G, network['reaction_database'], network['pathway_database'] = build_graph(
        network=network['reaction_database'],
        pathway_database=network['pathway_database'],
        species_reference=network['species_database'],
        name_reference=network['name_database'],
        protein_reference=protein_dictionary,
        chebi_mapper=chebi_mapper,
        uniprot_reference=network['uniprot_synonyms'],
        complexes=network['complex_dictionary'],
        species_id=species_id,
        gene_reference=network['ensembl_synonyms'],
        compartment_reference=network['compartment_dictionary'],
        component_database=network['components_database'])
        #additional_reactions=args_dict['additional_reactions'])
    progress_feed(args_dict, "model", 9)

    # For gene and protein components, add section to reaction database
    #for additional_components and list
    # Pull those in with everything else in JS

    # Overlay data and stats, calculate heatmap values for p-value
    # and expression value
    print('Mapping user data...')
    degree_dictionary = compile_node_degrees(
        graph=G)

    name_reference = {}
    for k, v in network['ensembl_synonyms'].items():
        name_reference[v] = k
        name_reference[k] = k
    for k, v in network['uniprot_synonyms'].items():
        name_reference[v] = k
        name_reference[k] = k

    # add any mapping IDs
    # Add synonyms
    # Change name to user provided if available
    metabolite_mapper = make_metabolite_synonym_dictionary(
        network=network,
        output_dir=args_dict['output'])

    G, max_value, max_stat, non_mappers = map_attributes(
        graph=G,
        data=data,
        stats=stats,
        name_reference=name_reference,
        degree_dictionary=degree_dictionary,
        chebi_mapper=chebi_mapper,
        metabolite_mapper=metabolite_mapper)
    progress_feed(args_dict, "graph", 5)

    if flag_data == True:
        max_value = 5
        max_stat = 1

    print('Broadcasting values where available...')
    categories = data.columns.tolist()

    if 'broadcast_genes' in args_dict and args_dict['broadcast_genes'] == True:
        broadcast_genes = True
    else:
        broadcast_genes = False

    G = broadcast_values(
        graph=G,
        categories=categories,
        max_value=max_value,
        max_stat=max_stat,
        broadcast_genes=broadcast_genes)
    progress_feed(args_dict, "graph", 10)

    print('Compiling collapsed reaction reference...')
    # Collapse reactions
    G, updated_reactions, changed_reactions, removed_reaction = collapse_nodes(
        graph=G,
        reaction_dictionary=network['reaction_database'],
        samples=len(categories),
        collapse_with_modifiers=args_dict['collapse_with_modifiers'])
    updated_pathway_dictionary = generate_updated_dictionary(
        original_database=network['pathway_database'],
        update_dictionary=changed_reactions,
        removed_reaction=removed_reaction)
    progress_feed(args_dict, "graph", 8)

    # Generate list of super pathways (those with more than 200 reactions)
    print('Compiling super pathways...')
    scale_factor = int(len(network['reaction_database'].keys()) * 0.0157)
    super_pathways = compile_pathway_degree(
        pathways=network['pathway_database'],
        scale_factor=scale_factor)

    motif_reaction_dictionary = make_motif_reaction_dictionary(
        network=network,
        updated_reactions=updated_reactions,
        updated_pathway_dictionary=updated_pathway_dictionary)

    mod_collapsed_pathways = {}
    for k, v in updated_pathway_dictionary.items():
        mod_collapsed_pathways[v['id']] = v

    # Final name mapping for chebi
    reverse_chebi = {}
    for k, v in network['chebi_synonyms'].items():
        reverse_chebi[v] = k
    for node in G.nodes():
        if 'chebi' in G.nodes()[node]['name'].lower():
            if G.nodes()[node]['name'] in reverse_chebi.keys():
                G.nodes()[node]['name'] = reverse_chebi[G.nodes()[node]['name']]

    # Export graph, pathway membership, pathway degree, other refs
    #unmapped['metabolites_unmapped'] = []
    #for x in non_mappers:
    #    if x not in data.index.tolist():
    #        unmapped['metabolites_unmapped'].append(x)
    #print(unmapped)

    print('Exporting graph...')
    args_dict["max_value"] = max_value
    args_dict["max_stat"] = max_stat
    args_dict["database_date"] = date.today().strftime('%Y-%m-%d')
    args_dict["curation_date"] = network["curation_date"]
    args_dict["reactome_version"] = network["reactome_version"]
    output_graph(
        graph=G,
        output_name=graph_name,
        pathway_dictionary=network['pathway_database'],
        collapsed_pathway_dictionary=updated_pathway_dictionary,
        super_pathways=super_pathways,
        reaction_dictionary=network['reaction_database'],
        collapsed_reaction_dictionary=updated_reactions,
        motif_reaction_dictionary=motif_reaction_dictionary,
        mod_collapsed_pathways=mod_collapsed_pathways,
        degree_dictionary=degree_dictionary,
        max_value=max_value,
        max_stat=max_stat,
        categories=categories,
        labels=args_dict['labels'],
        blacklist=args_dict['blacklist'],
        metadata=args_dict,
        unmapped=unmapped)
    print('Graphing complete.')
    progress_feed(args_dict, "graph", 2)

    return graph_name
