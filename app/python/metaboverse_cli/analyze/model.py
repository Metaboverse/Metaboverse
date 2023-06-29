"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) 2022 Metaboverse

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
from __future__ import print_function
from networkx.readwrite import json_graph
from collections import Counter
from datetime import date
from scipy.stats import gmean
import networkx as nx
import pandas as pd
import numpy as np
import zipfile
import pickle
import math
import json
import re
import os

"""Import internal dependencies
"""
try:
    from analyze.collapse import collapse_nodes
    from analyze.collapse import generate_updated_dictionary
    from analyze.mpl_colormaps import get_mpl_colormap
    from analyze.utils import convert_rgba, remove_defective_reactions
    from utils import progress_feed, track_progress, get_metaboverse_cli_version
except:
    import importlib.util
    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "analyze", "collapse.py"
                     ))
    spec = importlib.util.spec_from_file_location("", module_path)
    collapse = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(collapse)
    collapse_nodes = collapse.collapse_nodes
    generate_updated_dictionary = collapse.generate_updated_dictionary

    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "analyze", "mpl_colormaps.py"
                     ))
    spec = importlib.util.spec_from_file_location("", module_path)
    mpl_colormaps = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mpl_colormaps)
    get_mpl_colormap = mpl_colormaps.get_mpl_colormap

    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "analyze", "utils.py"
                     ))
    spec = importlib.util.spec_from_file_location("", module_path)
    analyze_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyze_utils)
    convert_rgba = analyze_utils.convert_rgba
    remove_defective_reactions = analyze_utils.remove_defective_reactions

    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "utils.py"))
    spec = importlib.util.spec_from_file_location("", module_path)
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    track_progress = utils.track_progress
    get_metaboverse_cli_version = utils.get_metaboverse_cli_version


CMAP = get_mpl_colormap('seismic')
REACTION_COLOR = (0.75, 0.75, 0.75, 1)
MISSING_COLOR = (1, 1, 1, 1)


def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None


def name_graph(
        output_file,
        species_id,
        template=False):
    """Name graph
    """

    if template == True:
        graph_name = species_id + '_template.mvrs'
    elif output_file[-5:].lower() == '.mvrs':
        graph_name = output_file
    elif output_file[-5:].lower() == '.eldb':
        graph_name = output_file
    else:
        graph_name = species_id + '_global_reactions.mvrs'

    return graph_name


def build_graph(
        args_dict,
        network,
        pathway_database,
        species_reference,
        name_reference,
        protein_reference,
        chebi_dictionary,
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

    counter = 0
    reaction_number = len(list(network.keys()))
    for reactome_id in network.keys():
        counter = track_progress(args_dict, counter, reaction_number, 5)
        G, network, key_hash, remove_keys = process_reactions(
            graph=G,
            reactome_id=reactome_id,
            network=network,
            species_reference=species_reference,
            name_reference=name_reference,
            protein_reference=protein_reference,
            chebi_dictionary=chebi_dictionary,
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
        chebi_dictionary,
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
    sort_name = ''.join(sorted(reaction_name.lower().replace(" ", "")))

    if sort_name in key_hash:
        remove_keys.append(reactome_id)
    else:
        key_hash.add(sort_name)

        reaction_id = network[reactome_id]['id']
        reaction_rev = network[reactome_id]['reversible']
        reaction_notes = network[reactome_id]['notes']

        reactants = network[reactome_id]['reactants']
        products = network[reactome_id]['products']
        modifiers = network[reactome_id]['modifiers']  # ordered list
        compartment_id = network[reactome_id]['compartment']
        try:
            compartment_name = compartment_reference[compartment_id]
        except:
            # for non-Reactome models where reactions do not have a compartment annotation
            compartment_name = ''

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
            if reactant in component_database \
            and component_database[reactant]['is'] != '':
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
                    chebi_dictionary=chebi_dictionary,
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
                    chebi_dictionary=chebi_dictionary,
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
                    chebi_dictionary=chebi_dictionary,
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
        id,  # node id
        map_id,  # used for mapping user data
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
        graph.nodes()[
            id]['compartment_display'] = compartment_reference[compartment]
    except:
        graph.nodes()[id]['compartment'] = 'none'
        graph.nodes()[id]['compartment_display'] = 'none'

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
        chebi_dictionary,
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

            elif x.lower in chebi_dictionary.keys():
                name = x
                map_id = chebi_dictionary[x]
                sub_type = 'metabolite_component'

            elif 'mi' in x.lower():
                map_id = name = x
                sub_type = 'mirna_component'

            else:
                map_id = name = x
                sub_type = 'other'

            if x.lower in chebi_dictionary.keys():
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


def reindex_data(
        data,
        stats):
    """Fix duplicate labels, etc
    """

    data_renamed = data.copy()
    data_renamed = data_renamed.loc[data_renamed.dropna(
        axis=0).index.drop_duplicates(keep=False)]
    d_cols = data_renamed.columns
    data_renamed[d_cols] = data_renamed[d_cols].apply(
        pd.to_numeric, errors='coerce')

    stats_renamed = stats.copy()
    stats_renamed = stats_renamed.loc[stats_renamed.dropna(
        axis=0).index.drop_duplicates(keep=False)]
    s_cols = stats_renamed.columns
    
    # Only force stats to numeric if not list type (list type indicates confidence interval arrays)
    if type(stats_renamed.iloc[0,0]) != list:
        stats_renamed[s_cols] = stats_renamed[s_cols].apply(
            pd.to_numeric, errors='coerce')

    if len(data_renamed.index.tolist()) != len(data.index.tolist()) \
            or len(stats_renamed.index.tolist()) != len(stats.index.tolist()):
        print('Warning: Duplicate row names were found. All duplicate data were \nremoved from downstream processing. Choose one row per name to \nmaintain the duplicated entity in downstream processing.')
        merge = list(set(data.index.tolist() + stats.index.tolist()))
        merge_after = list(
            set(data_renamed.index.tolist() + stats_renamed.index.tolist()))
        print("\n Unmapped entities:")
        for x in [y for y in merge if y not in merge_after]:
            print("\t-> " + str(x))

    return data_renamed, stats_renamed


def gather_synonyms(
        map_id,
        init_syns,
        metabolite_mapper,
        uniprot_mapper,
        ignore_enantiomers):

    if map_id in uniprot_mapper:
        init_syns.append(uniprot_mapper[map_id])
        init_syns.append(uniprot_mapper[map_id].lower())
        _u = ''.join(
            c.lower() for c in str(uniprot_mapper[map_id]) if c.isalnum())
        init_syns.append(_u)

    parsed_syns = set()
    for s in init_syns:
        parsed_syns.add(s)
        parsed_syns.add(s.lower())
        _s = ''.join(c.lower() for c in str(s) if c.isalnum())
        parsed_syns.add(_s)

        right_splits = [
            ', ',
            ' monocation',
            ' anion',
            ' monoanion',
            ', potassium',
            ', aluminum']
        left_splits = [
            'hydrogen ']
        for r in right_splits:
            if r in s.lower():
                parsed_syns.add(s.lower().split(r)[0])
        for l in left_splits:
            if l in s.lower():
                parsed_syns.add(s.lower().split(l)[1])

        if ignore_enantiomers == True:
            # if _s[0] == 'l' or _s[0] == 'd':
            #    parsed_syns.add(_s[1:])
            if s[0:2] == 'L-' or s[0:2] == 'D-':
                parsed_syns.add(s[2:])
            if s.lower()[0:2] == 'l-' or s.lower()[0:2] == 'd-':
                parsed_syns.add(s.lower()[2:])

    # Step 1: If ignore_enantiomers, trim off from beginning in data table and dictionaries, but search both with and without
    mapper_id = None
    check_keys = []
    search_keys = []
    log_keys = []
    for p in list(parsed_syns):
        if len(p) <= 1:  # ignore electrons, etc in mapping
            pass
        elif p in metabolite_mapper['mapping_dictionary']:
            mapper_id = metabolite_mapper['mapping_dictionary'][p]
        elif map_id in uniprot_mapper:
            if uniprot_mapper[map_id] in metabolite_mapper['mapping_dictionary']:
                mapper_id = metabolite_mapper['mapping_dictionary'][uniprot_mapper[map_id]]
        else:
            pass
        if mapper_id != None:
            log_keys.append(mapper_id)
    if len(log_keys) > 0:
        mapper_id = [word for word, word_count in Counter(
            log_keys).most_common(1)][0]
    else:
        mapper_id = None

    parsed_syns_list = list(parsed_syns)
    if mapper_id != None:
        parsed_syns_list.append(mapper_id)
    if mapper_id in metabolite_mapper['hmdb_dictionary']:
        for m in metabolite_mapper['hmdb_dictionary'][mapper_id]:
            parsed_syns_list.append(m)

    return mapper_id, parsed_syns_list


def prepare_mapping_data(graph, data, stats):
    """
    """

    n = len(data.columns.tolist())

    # Re-index data and stats
    data_renamed, stats_renamed = reindex_data(data, stats)
    data_max = abs(data_renamed).max().max()
    
    # Allow for lists of stats
    if type(stats_renamed.iloc[0,0]) == list:
        stats_logged = -1 
        stats_max = -1 
    else:
        stats_logged = -1 * np.log10(stats_renamed + 1e-100)
        stats_max = abs(stats_logged).max().max()

    # Make data dict with values and whether or not used
    temp_idx = [
        ''.join(
            c.lower() for c in str(i) if c.isalnum())
        for i in data_renamed.index.tolist()]
    temp_idx_set = set(temp_idx)

    # Get cross-species CHEBI synonyms
    chebi_mapping = {}
    for current_id in list(graph.nodes()):
        map_id = graph.nodes()[current_id]['map_id']
        name = graph.nodes()[current_id]['name']
        if 'chebi' in map_id.lower():
            if name in chebi_mapping:
                chebi_mapping[name].add(map_id)
            else:
                chebi_mapping[name] = set()
                chebi_mapping[name].add(map_id)

    return data_renamed, stats_renamed, data_max, stats_max, n, \
        temp_idx, temp_idx_set, chebi_mapping


def map_attributes(
        args_dict,
        graph,
        data,
        stats,
        name_reference,
        degree_dictionary,
        chebi_dictionary,
        chebi_synonyms,
        uniprot_mapper,
        metabolite_mapper,
        ignore_enantiomers=True):
    """Data overlay
    - Map repo id to species_id
    - If a node is a complex, take average of neighbors that are not
    To do:
    - Currently, many metabolites that should map are not found in name
    database
    """

    print("Mapping data onto network...")
    print("Input data dimensions: " + str(data.shape))
    data_renamed, stats_renamed, data_max, stats_max, n, \
    temp_idx, temp_idx_set, chebi_mapping = prepare_mapping_data(
        graph=graph,
        data=data,
        stats=stats)
    print("Pre-processed data dimensions: " + str(data_renamed.shape))
    
    mapped_nodes = []

    counter = 0
    node_number = len(list(graph.nodes()))

    # Map values to nodes
    for current_id in list(graph.nodes()):
        counter = track_progress(args_dict, counter, node_number, 5)

        x = current_id
        map_id = graph.nodes()[current_id]['map_id']
        backup_mapper = graph.nodes()[current_id]['name']

        # Add degree
        graph.nodes()[x]['degree'] = degree_dictionary[current_id]
        graph.nodes()[x]['synonyms'] = []

        # Add synonyms to node
        if map_id != 'none' \
        and map_id in name_reference \
        and (graph.nodes()[x]['type'] == 'protein_component'
        or graph.nodes()[x]['type'] == 'gene_component'):
            graph.nodes()[x]['synonyms'].append(map_id)

        elif map_id != 'none' \
        and backup_mapper in name_reference \
        and (graph.nodes()[x]['type'] == 'protein_component'
                 or graph.nodes()[x]['type'] == 'gene_component'):
            graph.nodes()[x]['synonyms'].append(backup_mapper)

        elif map_id != 'none' \
        and map_id in chebi_synonyms:
            graph.nodes()[x]['type'] = 'metabolite_component'
            graph.nodes()[x]['synonyms'].append(map_id)

            # Step 0: Get all CHEBI
            init_syns = chebi_synonyms[map_id]
            for s in init_syns:
                graph.nodes()[x]['synonyms'].append(s)

            # Step 2: Get initial CHEBI synonyms
            _mapper, _synonyms = gather_synonyms(
                map_id,
                init_syns,
                metabolite_mapper,
                uniprot_mapper,
                ignore_enantiomers
            )
            graph.nodes()[x]['hmdb_mapper'] = _mapper

            if len(_synonyms) > 0:
                for _s in _synonyms:
                    graph.nodes()[x]['synonyms'].append(_s)

        else:
            graph.nodes()[x]['synonyms'] = [x]
            
        if graph.nodes()[x]['sub_type'] == 'reaction':
            graph.nodes()[x]['type'] = 'reaction'
            colors = [REACTION_COLOR for x in range(n)]
            graph.nodes()[x]['values'] = [None for x in range(n)]
            graph.nodes()[x]['values_rgba'] = colors
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=colors)
            graph.nodes()[x]['stats'] = [None for x in range(n)]

        elif map_id in set(data_renamed.index.tolist()) \
        and map_id in set(stats_renamed.index.tolist()) \
        and map_id != 'none' \
        and len(map_id) > 1 \
        and graph.nodes()[x]['type'] != 'metabolite_component':
            graph.nodes()[x]['user_label'] = map_id
            graph.nodes()[x]['values'] = data_renamed.loc[map_id].tolist()
            graph.nodes()[x]['values_rgba'] = extract_value(
                value_array=data_renamed.loc[map_id].tolist(),
                max_value=data_max)
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['values_rgba'])
            graph.nodes()[x]['stats'] = stats_renamed.loc[map_id].tolist()
            mapped_nodes.append(map_id)

        elif backup_mapper in set(data_renamed.index.tolist()) \
        and backup_mapper in set(stats_renamed.index.tolist()) \
        and backup_mapper != 'none' \
        and len(backup_mapper) > 1 \
        and graph.nodes()[x]['type'] != 'metabolite_component':
            graph.nodes()[x]['user_label'] = backup_mapper
            graph.nodes()[
                x]['values'] = data_renamed.loc[backup_mapper].tolist()
            graph.nodes()[x]['values_rgba'] = extract_value(
                value_array=data_renamed.loc[backup_mapper].tolist(),
                max_value=data_max)
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=graph.nodes()[x]['values_rgba'])
            graph.nodes()[
                x]['stats'] = stats_renamed.loc[backup_mapper].tolist()
            mapped_nodes.append(backup_mapper)

        # elif map_id in chebi_synonyms
        elif ('chebi' in map_id.lower() or map_id in chebi_synonyms) \
        and graph.nodes()[x]['type'] == 'metabolite_component':
            _idx = None
            all_synonyms = graph.nodes()[x]['synonyms']
            if graph.nodes()[x]['name'] in chebi_mapping:
                all_synonyms.extend(chebi_mapping[graph.nodes()[x]['name']])

            for a in all_synonyms:
                _a = ''.join(c.lower() for c in a if c.isalnum())
                if _a in temp_idx_set:
                    _idx = data_renamed.index.tolist()[
                        temp_idx.index(_a)]
                elif ignore_enantiomers == True:
                    if 'D-' + str(a) in data_renamed.index.tolist():
                        _idx = 'D-' + str(a)
                    if 'd-' + str(a) in data_renamed.index.tolist():
                        _idx = 'd-' + str(a)
                    if 'L-' + str(a) in data_renamed.index.tolist():
                        _idx = 'L-' + str(a)
                    if 'l-' + str(a) in data_renamed.index.tolist():
                        _idx = 'l-' + str(a)
                    if 'N-' + str(a) in data_renamed.index.tolist():
                        _idx = 'N-' + str(a)
                    if 'n-' + str(a) in data_renamed.index.tolist():
                        _idx = 'n-' + str(a)
                else:
                    pass

            if graph.nodes()[x]['hmdb_mapper'] != None:
                graph.nodes()[x]['synonyms'] = metabolite_mapper['display_dictionary'][graph.nodes()[
                    x]['hmdb_mapper']]

            if _idx != None \
            and len(_idx) > 1:
                graph.nodes()[x]['user_label'] = _idx
                graph.nodes()[x]['values'] = data_renamed.loc[_idx].tolist()
                graph.nodes()[x]['values_rgba'] = extract_value(
                    value_array=data_renamed.loc[_idx].tolist(),
                    max_value=data_max)
                graph.nodes()[x]['values_js'] = convert_rgba(
                    rgba_tuples=graph.nodes()[x]['values_rgba'])
                graph.nodes()[x]['stats'] = stats_renamed.loc[_idx].tolist()
                mapped_nodes.append(_idx)
            elif map_id in set(data_renamed.index.tolist()) \
            and map_id in set(stats_renamed.index.tolist()) \
            and map_id != 'none' \
            and len(map_id) > 1:
                graph.nodes()[x]['user_label'] = map_id
                graph.nodes()[x]['values'] = data_renamed.loc[map_id].tolist()
                graph.nodes()[x]['values_rgba'] = extract_value(
                    value_array=data_renamed.loc[map_id].tolist(),
                    max_value=data_max)
                graph.nodes()[x]['values_js'] = convert_rgba(
                    rgba_tuples=graph.nodes()[x]['values_rgba'])
                graph.nodes()[x]['stats'] = stats_renamed.loc[map_id].tolist()
                mapped_nodes.append(map_id)
            else:
                colors = [MISSING_COLOR for x in range(n)]
                graph.nodes()[x]['values'] = [None for x in range(n)]
                graph.nodes()[x]['values_rgba'] = colors
                graph.nodes()[x]['values_js'] = convert_rgba(
                    rgba_tuples=colors)
                graph.nodes()[x]['stats'] = [None for x in range(n)]

        else:
            colors = [MISSING_COLOR for x in range(n)]
            graph.nodes()[x]['values'] = [None for x in range(n)]
            graph.nodes()[x]['values_rgba'] = colors
            graph.nodes()[x]['values_js'] = convert_rgba(
                rgba_tuples=colors)
            graph.nodes()[x]['stats'] = [None for x in range(n)]

    non_mappers = [x for x in data.index.tolist() if x not in mapped_nodes]

    return graph, data_max, stats_max, non_mappers


def extract_value(
        value_array,
        max_value,
        type="value"):
    """Extract expression value
    """

    def get_key_value(d, key):
        if key in d:
            return d[key]
        else:
            return d[max([x for x in d.keys() if x < key])]

    rgba = []
    for x in value_array:
        position = (x + max_value) / (2 * max_value)
        rgba_tuple = get_key_value(CMAP, round(position, 3))
        rgba.append(tuple(rgba_tuple))

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
        blocklist,
        species_blocklist,
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
    data['blocklist'] = blocklist
    data['species_blocklist'] = species_blocklist
    data['metadata'] = metadata
    data['unmapped'] = unmapped

    with open(output_name, 'w') as f:
        json.dump(data, f, indent=4)  # Parse out as array for javascript


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
    for k, v in deg_dict:

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

        protein_vals.append(median(pos))

    return protein_vals


def infer_protein_stats(stats, length, stat_type="float"):

    # Do not infer confidence intervals
    if stat_type == "array":
        protein_stats = [None for x in range(length)]
    else:
        protein_stats = []
        for i in range(length):

            pos = []

            for j in range(len(stats)):

                if stats[j][i] != None:
                    pos.append(stats[j][i])

            if len(pos) == 1:
                this_stat = pos[0]
            else:
                this_stat = (math.e * gmean(pos))

            if this_stat > 1.0:
                this_stat = 1.0
            protein_stats.append(this_stat)

    return protein_stats


def broadcast_values(
        args_dict,
        graph,
        categories,
        max_value,
        max_stat,
        broadcast_genes=True,
        broadcast_metabolites=True,
        stat_type="float"):
    """
    """

    u = graph.to_undirected()
    length = len(categories)

    if broadcast_genes == True:
        counter = 0
        node_number = len(list(graph.nodes()))

        for x in graph.nodes():
            counter = track_progress(args_dict, counter, node_number, 5)

            try:
                graph.nodes()[x]['values']
            except:
                print(graph.nodes()[x])
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

                if graph.nodes()[x]['sub_type'] == 'protein_component':

                    gene_values = []
                    gene_stats = []
                    for neighbor in u[x]:

                        if graph.nodes()[neighbor]['sub_type'] == 'gene':
                            gene_values.append(
                                graph.nodes()[neighbor]['values'])
                            gene_stats.append(graph.nodes()[neighbor]['stats'])

                    # Remove None and avg
                    gene_values = remove_nulls(gene_values)
                    gene_stats = remove_nulls(gene_stats)

                    # Mark as inferred if this is possible
                    if gene_values != []:
                        inferred_values = infer_protein_values(
                            gene_values, length)

                        graph.nodes()[x]['inferred'] = 'true'
                        graph.nodes()[x]['values'] = inferred_values
                        graph.nodes()[x]['values_rgba'] = extract_value(
                            value_array=inferred_values,
                            max_value=max_value)
                        graph.nodes()[x]['values_js'] = convert_rgba(
                            rgba_tuples=graph.nodes()[x]['values_rgba'])

                    if gene_stats != []:
                        inferred_stats = infer_protein_stats(
                            gene_stats, length, stat_type=stat_type)

                        graph.nodes()[x]['inferred'] = 'true'
                        graph.nodes()[x]['stats'] = inferred_stats
    else:
        progress_feed(args_dict, "graph", 5)


    counter = 0
    node_number = len(list(graph.nodes()))

    for x in graph.nodes():
        counter = track_progress(args_dict, counter, node_number, 5)

        try:
            graph.nodes()[x]['values']
        except:
            print(graph.nodes()[x])
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
                    elif broadcast_metabolites == True \
                            and graph.nodes()[neighbor]['type'] == 'metabolite_component':
                        gene_values.append(graph.nodes()[neighbor]['values'])
                        gene_stats.append(graph.nodes()[neighbor]['stats'])
                    else:
                        pass

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
                    inferred_stats = infer_protein_stats(
                        gene_stats, length, stat_type=stat_type)

                    graph.nodes()[x]['inferred'] = 'true'
                    graph.nodes()[x]['stats'] = inferred_stats

    return graph


def make_motif_reaction_dictionary(
        network,
        updated_reactions,
        updated_pathway_dictionary):

    simp_dict = {}
    for k, v in network['reaction_database'].items():
        simp_dict[k] = []

    for k, v in network['pathway_database'].items():

        for x in network['pathway_database'][k]['reactions']:
            try:
                simp_dict[x].append(network['pathway_database'][k]['id'])
            except:
                simp_dict[x] = []
                simp_dict[x].append(network['pathway_database'][k]['id'])

    motif_reaction_dictionary = {}
    for k, v in updated_reactions.items():
        motif_reaction_dictionary[k] = []

    for k, v in updated_pathway_dictionary.items():

        for x in updated_pathway_dictionary[k]['reactions']:
            try:
                motif_reaction_dictionary[x].append(
                    updated_pathway_dictionary[k]['id'])
            except:
                motif_reaction_dictionary[x] = []
                motif_reaction_dictionary[x].append(
                    updated_pathway_dictionary[k]['id'])

    for k, v in motif_reaction_dictionary.items():
        if len(v) == 0:
            comps = k.split('_reaction_')
            comps = [x.replace('reaction_', '') for x in comps]

            unity_paths = []
            for y in comps:
                _rxn = 'reaction_' + y
                _paths = simp_dict[_rxn]
                unity_paths.append(_paths)

            # For collapsed reactions, only include pathways that all
            # component reactions are a member of
            result = set(unity_paths[0])
            for u in unity_paths[1:]:
                result.intersection_update(set(u))
            motif_reaction_dictionary[k].append(list(result))

    return motif_reaction_dictionary


def load_metabolite_synonym_dictionary(
        dir=os.path.join(os.path.dirname(__file__), 'data'),
        file='metabolite_mapping.pickle'):

    print("Reading metabolite mapper...")
    with zipfile.ZipFile(dir + os.path.sep + file + '.zip', 'r') as zip_ref:
        metabolite_mapper = pickle.load(
            zip_ref.open(file)
        )

    return metabolite_mapper


def build_name_reference(
        ensembl,
        uniprot):

    name_reference = {}
    for k, v in ensembl.items():
        if 'phospho-' in v and '-phospho-' not in v:
            v = v.replace('phospho-', '')
        if '(' in v and ')' in v:
            v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
        if '  ' in v:
            v = v.replace('  ', ' ')
        v = v.strip()
        name_reference[v] = k
        name_reference[k] = k
    for k, v in uniprot.items():
        if 'phospho-' in v and '-phospho-' not in v:
            v = v.replace('phospho-', '')
        if '(' in v and ')' in v:
            v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
        if '  ' in v:
            v = v.replace('  ', ' ')
        v = v.strip()
        name_reference[v] = k
        name_reference[k] = k

    return name_reference


def build_chebi_reference(
        chebi,
        uniprot):

    chebi_dictionary = {}
    for k, v in chebi.items():
        _k = ''.join(c.lower() for c in str(k) if c.isalnum())
        chebi_dictionary[_k] = v
        chebi_dictionary[k] = v
    for k, v in uniprot.items():
        _k = ''.join(c.lower() for c in str(k) if c.isalnum())
        chebi_dictionary[_k] = v
        chebi_dictionary[k] = v

    return chebi_dictionary


def load_references(
        args_dict,
        ensembl,
        uniprot,
        chebi,
        uniprot_metabolites):
    """Load and prepare reference databases
    """

    # Prepare uniprot to ensembl name mapper
    reverse_genes = {v: k for k, v in ensembl.items()}
    protein_dictionary = uniprot_ensembl_reference(
        uniprot_reference=uniprot,
        ensembl_reference=reverse_genes)
    progress_feed(args_dict, "graph", 1)

    chebi_dictionary = build_chebi_reference(
        chebi=chebi,
        uniprot=uniprot_metabolites)

    name_reference = build_name_reference(
        ensembl=ensembl,
        uniprot=uniprot)

    uniprot_mapper = {}
    for k, v in uniprot_metabolites.items():
        uniprot_mapper[v] = k

    return reverse_genes, protein_dictionary, chebi_dictionary, \
        name_reference, uniprot_mapper


def __template__(
        args_dict,
        network,
        species_id,
        output_file):
    """
    """

    print('Preparing metadata...')
    # Generate output file name
    graph_name = name_graph(
        output_file=output_file,
        species_id=species_id,
        template=True)

    print('Preparing references...')
    reverse_genes, protein_dictionary, chebi_dictionary, \
        name_reference, uniprot_mapper = load_references(
            args_dict=args_dict,
            ensembl=network['ensembl_synonyms'],
            uniprot=network['uniprot_synonyms'],
            chebi=network['chebi_mapper'],
            uniprot_metabolites=network['uniprot_metabolites'])
    metabolite_mapper = load_metabolite_synonym_dictionary()

    # Generate graph and name mapping
    print('Building network...')
    G, network['reaction_database'], network['pathway_database'] = build_graph(
        args_dict=args_dict,
        network=network['reaction_database'],
        pathway_database=network['pathway_database'],
        species_reference=network['species_database'],
        name_reference=network['name_database'],
        protein_reference=protein_dictionary,
        chebi_dictionary=chebi_dictionary,
        uniprot_reference=network['uniprot_synonyms'],
        complexes=network['complex_dictionary'],
        species_id=species_id,
        gene_reference=network['ensembl_synonyms'],
        compartment_reference=network['compartment_dictionary'],
        component_database=network['components_database'])
    # additional_reactions=args_dict['additional_reactions'])
    progress_feed(args_dict, "graph", 1)

    # Generate list of super pathways (those with more than 200 reactions)
    print('Compiling super pathways...')
    scale_factor = int(len(network['reaction_database'].keys()) * 0.0157)
    super_pathways = compile_pathway_degree(
        pathways=network['pathway_database'],
        scale_factor=scale_factor)

    print('Compiling network degree database...')
    degree_dictionary = compile_node_degrees(
        graph=G)

    print('Exporting template...')
    args_dict["max_value"] = 0
    args_dict["max_stat"] = 1
    args_dict["model_version"] = get_metaboverse_cli_version()
    args_dict["model_date"] = date.today().strftime('%Y-%m-%d')
    args_dict["curation_date"] = network["curation_date"]
    args_dict["curation_version"] = network["metaboverse-curate_version"]
    args_dict['template_url'] = os.path.join(args_dict['output'], graph_name)
    args_dict['template_version'] = get_metaboverse_cli_version()
    args_dict['template_date'] = date.today().strftime('%Y-%m-%d')

    output_graph(
        graph=G,
        output_name=os.path.join(args_dict['output'], graph_name),
        pathway_dictionary=network['pathway_database'],
        collapsed_pathway_dictionary=network['pathway_database'],
        super_pathways=super_pathways,
        reaction_dictionary=network['reaction_database'],
        collapsed_reaction_dictionary=network['reaction_database'],
        motif_reaction_dictionary=network['reaction_database'],
        mod_collapsed_pathways={},
        degree_dictionary=degree_dictionary,
        max_value=0,
        max_stat=1,
        categories=[],
        labels=args_dict['labels'],
        blocklist=args_dict['blocklist'],
        species_blocklist=[],
        metadata=args_dict,
        unmapped=[])
    print('Graphing complete.')

    return G, args_dict, network, name_reference, degree_dictionary, \
        super_pathways, chebi_dictionary, uniprot_mapper, metabolite_mapper


def __model__(
        graph,
        args_dict,
        network,
        data,
        stats,
        species_id,
        output_file,
        neighbors_dictionary,
        name_reference,
        degree_dictionary,
        chebi_dictionary,
        uniprot_mapper,
        metabolite_mapper,
        super_pathways,
        unmapped,
        flag_data=False):
    """Generate graph object for visualization
    """

    # Generate output file name
    graph_name = name_graph(
        output_file=output_file,
        species_id=species_id)

    print('Post-processing graph metadata...')
    # Remove disease reactions that are "defective" from reaction collapse
    no_defective_reactions = remove_defective_reactions(
        network=network)

    print('Mapping user data...')
    G, max_value, max_stat, non_mappers = map_attributes(
        args_dict=args_dict,
        graph=graph,
        data=data,
        stats=stats,
        name_reference=name_reference,
        degree_dictionary=degree_dictionary,
        chebi_dictionary=chebi_dictionary,
        chebi_synonyms=network['chebi_synonyms'],
        uniprot_mapper=uniprot_mapper,
        metabolite_mapper=metabolite_mapper)
    
    print('Outputting unmapped metabolomics values (if any exist)...')
    if args_dict['metabolomics'].lower() != 'none':
        m_data = pd.read_csv(
            args_dict['metabolomics'],
            sep='\t',
            index_col=0)

        m_non_mapper = m_data[m_data.index.isin(non_mappers)]
        print("\t- Outputting " + str(len(m_non_mapper.index.tolist())) + " unmapped metabolites for reference.")
        if len(m_non_mapper.index.tolist()) > 0:
            m_non_mapper.to_csv(
                args_dict['metabolomics'][:-4] + '_unmapped.txt',
                sep='\t')

    print('Broadcasting values where available...')
    if 'broadcast_genes' in args_dict \
            and args_dict['broadcast_genes'] == True:
        broadcast_genes = True
    else:
        broadcast_genes = False

    if 'broadcast_metabolites' in args_dict \
            and args_dict['broadcast_metabolites'] == True:
        broadcast_metabolites = True
    else:
        broadcast_metabolites = False

    if flag_data == True:
        max_value = 5
        max_stat = 1

    categories = data.columns.tolist()
    if type(stats.iloc[0,0]) == list:
        args_dict["stat_type"] = 'array'
    else:
        args_dict["stat_type"] = 'float'
    G = broadcast_values(
        args_dict=args_dict,
        graph=G,
        categories=categories,
        max_value=max_value,
        max_stat=max_stat,
        broadcast_genes=broadcast_genes,
        broadcast_metabolites=broadcast_metabolites, 
        stat_type=args_dict["stat_type"])
    progress_feed(args_dict, "graph", 5)

    print('Compiling collapsed reaction reference...')
    # Get hub threshold
    degrees = []
    for k in degree_dictionary.keys():
        if 'reaction' not in k:
            degrees.append(degree_dictionary[k])
    if len(degrees) > 0:
        degree_threshold = np.percentile(degrees, 98)
    else:
        degree_threshold = 0


    if 'blocklist' in args_dict \
            and isinstance(args_dict['blocklist'], str):
        named_blocklist = args_dict['blocklist'].replace(' ', '').split(',')
    elif 'blocklist' in args_dict \
            and isinstance(args_dict['blocklist'], list):
        named_blocklist = args_dict['blocklist']
    else:
        named_blocklist = []

    species_blocklist = []
    for b in named_blocklist:
        if b in network['name_database']:
            species_blocklist.append(network['name_database'][b])

    # Collapse reactions
    G, updated_reactions, changed_reactions, \
    removed_reaction = collapse_nodes(
        args_dict=args_dict,
        graph=G,
        reaction_dictionary=no_defective_reactions,
        neighbors_dictionary=neighbors_dictionary,
        degree_dictionary=degree_dictionary,
        samples=len(categories),
        collapse_with_modifiers=args_dict['collapse_with_modifiers'],
        blocklist=species_blocklist,
        degree_threshold=degree_threshold,
        collapse_threshold=args_dict['collapse_threshold'])
    updated_pathway_dictionary = generate_updated_dictionary(
        original_database=network['pathway_database'],
        update_dictionary=changed_reactions,
        removed_reaction=removed_reaction)
    progress_feed(args_dict, "graph", 2)

    motif_reaction_dictionary = make_motif_reaction_dictionary(
        network=network,
        updated_reactions=updated_reactions,
        updated_pathway_dictionary=updated_pathway_dictionary)

    mod_collapsed_pathways = {}
    for k, v in updated_pathway_dictionary.items():
        mod_collapsed_pathways[v['id']] = v

    # Final name mapping for chebi
    reverse_chebi = {}
    for k, v in network['chebi_mapper'].items():
        reverse_chebi[v] = k
    for node in G.nodes():
        if 'chebi' in G.nodes()[node]['name'].lower():
            if G.nodes()[node]['name'] in reverse_chebi.keys():
                G.nodes()[node]['name'] = reverse_chebi[G.nodes()[node]['name']]

    print('Exporting graph...')
    args_dict["max_value"] = max_value
    args_dict["max_stat"] = max_stat 
    args_dict["curation_url"] = os.path.join(
        args_dict['output'],
        args_dict['curation'])
    args_dict["curation_version"] = network["metaboverse-curate_version"]
    args_dict["curation_date"] = network["curation_date"]
    args_dict["database_version"] = network["database_version"]
    args_dict["model_version"] = get_metaboverse_cli_version()
    args_dict["model_date"] = date.today().strftime('%Y-%m-%d')
    args_dict['neighbors_url'] = neighbors_dictionary['nbdb-Metaboverse-url']
    args_dict['neighbors_version'] = neighbors_dictionary['nbdb-Metaboverse-version']
    args_dict['neighbors_date'] = neighbors_dictionary['nbdb-Metaboverse-date']

    output_graph(
        graph=G,
        output_name=os.path.join(args_dict['output'], graph_name),
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
        blocklist=named_blocklist,
        species_blocklist=species_blocklist,
        metadata=args_dict,
        unmapped=non_mappers)
    print('Graphing complete.')
    progress_feed(args_dict, "graph", 1)

    return os.path.join(args_dict['output'], graph_name)
