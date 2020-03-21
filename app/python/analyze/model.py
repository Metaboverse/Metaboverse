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
from python.analyze.collapse import collapse_nodes
from python.analyze.collapse import generate_updated_dictionary
from python.analyze.utils import convert_rgba
from python.utils import progress_feed

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
        species_reference,
        name_reference,
        protein_reference,
        uniprot_reference,
        complexes,
        species_id,
        gene_reference,
        compartment_reference):
    """Build graph
    - Add nodes and edges
    - Map names to objects in the graph for display
    - Calculate degree for each node
    """

    # Initialize graph object
    G = nx.DiGraph()

    for reactome_id in network.keys():

        G, network = process_reactions(
            graph=G,
            reactome_id=reactome_id,
            network=network,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            uniprot_reference=uniprot_reference,
            complex_reference=complexes,
            species_id=species_id,
            gene_reference=gene_reference,
            compartment_reference=compartment_reference)

    return G, network

def process_reactions(
        graph,
        reactome_id,
        network,
        species_reference,
        name_reference,
        protein_reference,
        uniprot_reference,
        complex_reference,
        species_id,
        gene_reference,
        compartment_reference):
    """
    """
    new_components = []

    # Get reaction name
    reaction_id = network[reactome_id]['id']
    reaction_name = network[reactome_id]['name']
    reaction_rev = network[reactome_id]['reversible']
    reaction_notes = network[reactome_id]['notes']

    reactants = network[reactome_id]['reactants']
    products = network[reactome_id]['products']
    modifiers = network[reactome_id]['modifiers'] # ordered list
    compartment_id = network[reactome_id]['compartment']
    compartment_name = compartment_reference[compartment_id]

    # Add reaction node
    graph.add_node(reaction_id)
    graph.nodes()[reaction_id]['id'] = reactome_id
    graph.nodes()[reaction_id]['name'] = reaction_name
    graph.nodes()[reaction_id]['reversible'] = reaction_rev
    graph.nodes()[reaction_id]['notes'] = reaction_notes
    graph.nodes()[reaction_id]['type'] = 'reaction'
    graph.nodes()[reaction_id]['sub_type'] = 'reaction'
    graph.nodes()[reaction_id]['compartment'] = compartment_id
    graph.nodes()[reaction_id]['compartment_display'] = compartment_name

    # Add vanilla element nodes and their edges
    for reactant in reactants:

        graph = add_node_edge(
            graph=graph,
            id=reactant,
            name=species_reference[reactant],
            reaction_membership=reaction_id,
            type='reactant',
            sub_type='reactant',
            reversible=reaction_rev,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference)

        graph, additional_components = check_complexes(
            graph=graph,
            id=reactant,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            uniprot_reference=uniprot_reference,
            gene_reference=gene_reference)
        for x in additional_components:
            new_components.append(x)

    for product in products:

        graph = add_node_edge(
            graph=graph,
            id=product,
            name=species_reference[product],
            reaction_membership=reaction_id,
            type='product',
            sub_type='product',
            reversible=reaction_rev,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference)

        graph, additional_components = check_complexes(
            graph=graph,
            id=product,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            uniprot_reference=uniprot_reference,
            gene_reference=gene_reference)
        for x in additional_components:
            new_components.append(x)

    for modifier in modifiers:

        # Extract modifier type
        # Formatted as a list of lists with first index of each sub list being
        # the species ID and the second index being the modifier type
        # ex: [['species_0', 'inhibitor'], ['species_1', 'catalyst']]
        # Labeling the edge should allow for differentiation between the same
        # modifier node acting as a catalyst or inhibitor
        id = modifier[0]
        type = modifier[1]

        graph = add_node_edge(
            graph=graph,
            id=id,
            name=species_reference[id],
            reaction_membership=reaction_id,
            type=type,
            sub_type='modifier',
            reversible='false',
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference)

        graph, additional_components = check_complexes(
            graph=graph,
            id=id,
            complex_reference=complex_reference,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            uniprot_reference=uniprot_reference,
            gene_reference=gene_reference)
        for x in additional_components:
            new_components.append(x)

    network[reactome_id]['additional_components'] = new_components

    return graph, network

def add_node_edge(
        graph,
        id, # node id
        name, # display name
        reaction_membership,
        type,
        sub_type,
        reversible,
        complex_reference,
        species_reference,
        name_reference,
        protein_reference):
    """Add node and edge information to graph
    """

    graph.add_node(id)
    graph.nodes()[id]['id'] = id
    graph.nodes()[id]['name'] = name
    graph.nodes()[id]['type'] = type
    graph.nodes()[id]['sub_type'] = sub_type
    graph.nodes()[id]['inferred'] = 'false'
    graph.nodes()[id]['compartment'] = 'none'

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
        graph,
        id,
        complex_reference,
        species_reference,
        name_reference,
        protein_reference,
        uniprot_reference,
        gene_reference):
    """Check if species being added is in complex dictionary
    - If record exists, add nodes and edges for the new relationship.
    - If record contains a UniProt ID, cross reference with Ensembl database
    - If complex, label true; else label as false
    """

    prot_ref = {}
    for k, v in uniprot_reference.items():
        prot_ref[k] = v
        prot_ref[v] = k
    uniprot_reference = prot_ref

    add_components = []

    if id in complex_reference.keys():
        graph.nodes()[id]['complex'] = 'true'

        participants = complex_reference[id]['participants']
        for p in participants.keys():

            for x in participants[p]:

                if p.lower() == 'chebi':
                    name = 'CHEBI:' + x
                    sub_type = 'metabolite_component'
                    component_id = name_reference[name]
                    add_components.append(component_id)
                    display_name = species_reference[component_id]

                else:

                    if p.lower() == 'uniprot':
                        name = uniprot_reference[x]
                        sub_type = 'protein_component'
                        component_id = x
                        add_components.append(x)
                        display_name = name
                    elif p.lower() == 'mirna':
                        name = x
                        sub_type = 'mirna_component'
                        component_id = name_reference[name]
                        add_components.append(component_id)
                        display_name = species_reference[component_id]
                    elif p.lower() == 'ensembl':
                        name = x
                        sub_type = 'gene_component'
                        component_id = x
                        add_components.append(component_id)
                        display_name = gene_reference[name]
                    else:
                        name = x
                        sub_type = 'other'
                        try:
                            component_id = name_reference[name]
                            display_name = species_reference[component_id]
                        except:
                            display_name = component_id = x
                        add_components.append(component_id)

                try:

                    graph = add_node_edge(
                        graph=graph,
                        id=component_id,
                        name=display_name,
                        reaction_membership=id,
                        type='complex_component',
                        sub_type=sub_type,
                        reversible='false',
                        complex_reference=complex_reference,
                        species_reference=species_reference,
                        name_reference=name_reference,
                        protein_reference=protein_reference)

                    if p.lower() == 'uniprot':

                        # Get protein's corresponding gene ID
                        gene = protein_reference[component_id]
                        add_components.append(gene)

                        graph = add_node_edge(
                            graph=graph,
                            id=gene,
                            name=gene_reference[gene],
                            reaction_membership=component_id,
                            type='gene_component',
                            sub_type='gene',
                            reversible='false',
                            complex_reference=complex_reference,
                            species_reference=species_reference,
                            name_reference=name_reference,
                            protein_reference=protein_reference)

                except:
                    pass
                    #print('Could not retrieve components for', name)

    else:
        graph.nodes()[id]['complex'] = 'false'

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

def map_attributes(
        graph,
        data,
        stats,
        name_reference,
        degree_dictionary):
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
    data_renamed = data.rename(index=network['name_database'])
    data_renamed = data_renamed.loc[data_renamed.index.dropna()]
    data_max = abs(data_renamed).max().max()

    stats_renamed = stats.copy()
    stats_renamed = stats.rename(index=network['name_database'])
    stats_renamed = stats_renamed.loc[stats_renamed.index.dropna()]
    stats_logged = -1 * np.log10(stats_renamed + 1e-100)

    stats_max = abs(stats_logged).max().max()

    data_dict = {}
    for index, row in data_renamed.iterrows():
        data_dict[index] = list(row)

    stats_dict = {}
    for index, row in stats_renamed.iterrows():
        stats_dict[index] = list(row)

    for x in list(graph.nodes()):

        current_id = graph.nodes()[x]['id']

        # Add degree
        graph.nodes()[x]['degree'] = degree_dictionary[current_id]

        if current_id in set(data_dict.keys()):

            graph.nodes()[x]['values'] = data_dict[current_id]
            graph.nodes()[x]['values_rgba'] = extract_value(
                value_array=data_dict[current_id],
                max_value=data_max)
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['values_rgba'])

            graph.nodes()[x]['stats'] = stats_dict[current_id]
            graph.nodes()[x]['stats_rgba'] = extract_value(
                value_array=stats_dict[current_id],
                max_value=stats_max,
                type="stats")
            graph.nodes()[x]['stats_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['stats_rgba'])

        else:

            if graph.nodes()[x]['type'] == 'reaction':
                colors = [reaction_color for x in range(n)]

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

    return graph, data_max, stats_max

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
        max_value,
        max_stat,
        categories):
    """Output graph and necessary metadata
    """

    data = json_graph.node_link_data(graph)
    data['pathway_dictionary'] = pathway_dictionary
    data['collapsed_pathway_dictionary'] = collapsed_pathway_dictionary
    data['super_pathways'] = super_pathways
    data['reaction_dictionary'] = reaction_dictionary
    data['collapsed_reaction_dictionary'] = collapsed_reaction_dictionary
    data['max_value'] = max_value
    data['max_stat'] = max_stat
    data['categories'] = categories

    with open(output_name, 'w') as f:
        json.dump(data, f, indent=4) # Parse out as array for javascript

def compile_pathway_degree(
        pathways):
    """Compile database of large pathways
    """

    super_pathways = {}

    for key in list(pathways.keys()):

        if len(pathways[key]["reactions"]) > 200:
            super_pathways[key] = pathways[key]

    return super_pathways

def compile_node_degrees(
        graph):
    """Retrieve degree metrics for each node
    """

    d = {}

    deg_dict = G.degree
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

def broadcast_values(
        graph,
        categories,
        max_value,
        max_stat):
    """
    """

    u = graph.to_undirected()
    length = len(categories)

    for x in graph.nodes():

        if None not in graph.nodes()[x]['values'] \
        and None not in graph.nodes()[x]['stats']:
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
                    inferred_stats = infer_protein_values(gene_stats, length)

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
                    inferred_stats = infer_protein_values(gene_stats, length)

                    graph.nodes()[x]['inferred'] = 'true'
                    graph.nodes()[x]['stats'] = inferred_stats
                    graph.nodes()[x]['stats_rgba'] = extract_value(
                        value_array=inferred_stats,
                        max_value=max_stat,
                        type="stats")
                    graph.nodes()[x]['stats_js'] = convert_rgba(
                        rgba_tuples=graph.nodes()[x]['stats_rgba'])


    return graph

def __main__(
        args_dict,
        network,
        data,
        stats,
        species_id,
        output_file):
    """Generate graph object for visualization
    """

    #############################
    def read_network(
            network_url):
        """Read in network from previous curation module
        - was provided as a URL to the file and saved to args_dict['network']  in "curate" sub-module
        """
        import pickle
        with open(network_url, 'rb') as network_file:
            network = pickle.load(network_file)

        return network

    network = read_network(
        network_url='/Users/jordan/Desktop/HSA_metaboverse_db.pickle')

    data = pd.read_csv(
        '/Users/jordan/Desktop/metaboverse/app/python/analyze/test/cat_data.txt',
        sep='\t',
        index_col=0)

    stats = pd.read_csv(
        '/Users/jordan/Desktop/metaboverse/app/python/analyze/test/cat_stats.txt',
        sep='\t',
        index_col=0)

    output_file = '/Users/jordan/Desktop/HSA_global_reactions.json'

    species_id = 'HSA'

    """Start of pancancer-necessary code
    To run, place all graph-making code within the for loop below
    """
    """
    cancers = [
    ['breast', 'brca'],
    ['colon', 'coad'],
    ['colon','read'],
    ['bladder','blca'],
    ['cervix','cesc'],
    ['kidney','kirc'],
    ['kidney','kirp'],
    ['liver','lihc'],
    ['lung','luad'],
    ['lung','lusc'],
    ['prostate','prad'],
    ['stomach','stad'],
    ['thyroid','thca'],
    ['uterus','ucec']]

    for x in [x[1] for x in cancers]:
        output_file = '/Users/jordan/Desktop/HSA_global_reactions_' + x + '.json'
        print(output_file)
        data = stats = pd.read_csv(
            '/Users/jordan/Desktop/cancer_network_topology/tables/' + x + '_z_table.txt',
            sep='\t',
            index_col=0)
    """


    #############################

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
    progress_feed(args_dict, "model", 5)

    # Generate graph
    # Name mapping
    print('Building network...')
    G, network['reaction_database'] = build_graph(
        network=network['reaction_database'],
        species_reference=network['species_database'],
        name_reference=network['name_database'],
        protein_reference=protein_dictionary,
        uniprot_reference=network['uniprot_synonyms'],
        complexes=network['complex_dictionary'],
        species_id=species_id,
        gene_reference=network['ensembl_synonyms'],
        compartment_reference=network['compartment_dictionary'])
    progress_feed(args_dict, "model", 10)

    # For gene and protein components, add section to reaction database
    #for additional_components and list
    # Pull those in with everything else in JS

    # Overlay data and stats, calculate heatmap values for p-value
    # and expression value
    print('Mapping user data...')
    degree_dictionary = compile_node_degrees(
        graph=G)
    G, max_value, max_stat = map_attributes(
        graph=G,
        data=data,
        stats=stats,
        name_reference=network['name_database'],
        degree_dictionary=degree_dictionary)
    progress_feed(args_dict, "graph", 5)

    print('Broadcasting values where available...')
    categories = data.columns.tolist()
    G = broadcast_values(
        graph=G,
        categories=categories,
        max_value=max_value,
        max_stat=max_stat)
    progress_feed(args_dict, "graph", 10)

    # Generate list of super pathways (those with more than 200 reactions)
    print('Compiling super pathways...')
    super_pathways = compile_pathway_degree(
        pathways=network['pathway_database'])


    print('Compiling collapsed reaction reference...')
    # Collapse reactions
    G, updated_reactions, changed_reactions = collapse_nodes(
        graph=G,
        reaction_dictionary=network['reaction_database'],
        samples=len(categories))
    updated_pathway_dictionary = generate_updated_dictionary(
        original_database=network['pathway_database'],
        update_dictionary=changed_reactions)
    progress_feed(args_dict, "graph", 8)

    ###
    # For pancancer
    """
    metabolism = network['pathway_database']['R-HSA-1430728']
    network['pathway_database'] = {}
    network['pathway_database']['R-HSA-1430728'] = metabolism

    with open('/Users/jordan/Desktop/cancer_network_topology/network/cancer_reactions.json') as f:
        file = json.load(f)
        for k,v in file.items():

            network['pathway_database'][k] = {}
            network['pathway_database'][k]['id'] = k
            network['pathway_database'][k]['name'] = k
            network['pathway_database'][k]['reactions'] = v
    """
    ###

    # Export graph, pathway membership, pathway degree, other refs
    print('Exporting graph...')
    output_graph(
        graph=G,
        output_name=graph_name,
        pathway_dictionary=network['pathway_database'],
        collapsed_pathway_dictionary=updated_pathway_dictionary,
        super_pathways=super_pathways,
        reaction_dictionary=network['reaction_database'],
        collapsed_reaction_dictionary=updated_reactions,
        max_value=max_value,
        max_stat=max_stat,
        categories=categories)
    print('Graphing complete.')
    progress_feed(args_dict, "graph", 2)
