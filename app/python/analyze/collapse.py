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
from analyze.utils import convert_rgba

def generate_updated_dictionary(
        original_database,
        update_dictionary,
        removed_reaction):
    """Back-reference collapsed reactions from pathways
    """

    updated_pathway_dictionary = {}
    for k, v in original_database.items():

        key = reactome = v['reactome']
        id = v['id']
        name = v['name']
        reactions = v['reactions']
        new_reactions = []

        for k, v in update_dictionary.items():
            if type(k) is tuple:
                l = list(k)
                if set(l) <= set(reactions):
                    if v not in removed_reaction:
                        new_reactions.append(v)
                        reactions = list(set(reactions) - set(v))
            else:
                l = [k]
                if l in reactions or k in reactions:
                    if v not in removed_reaction:
                        new_reactions.append(v)
                        reactions = list(set(reactions) - set(v))

        updated_pathway_dictionary[reactome] = {
            'id': id,
            'reactome': reactome,
            'name': name,
            'reactions': new_reactions
        }

    return updated_pathway_dictionary

def find_values(
        graph,
        reaction_dictionary,
        neighbor,
        current_inputs,
        collapse_with_modifiers):
    """Are there any values for either side of the reaction?
    - There can only be one side with values for this situation to be valid in
    this context
    """

    eval_items = [] # What will be returned
    values = []

    # Check if values exist for the new side
    # Allow modifier values to count
    if set(current_inputs) == set(reaction_dictionary[neighbor]['reactants']):
        side_key = 'products'
        if collapse_with_modifiers == True or collapse_with_modifiers == "True":
            mod_key = 'catalyst'
        else:
            # Using this to only consider main components during
            # reaction collapse
            mod_key = 'impossible_key_to_match'
    else:
        side_key = 'reactants'
        if collapse_with_modifiers == True or collapse_with_modifiers == "True":
            mod_key = 'inhibitor'
        else:
            # Using this to only consider main components during
            # reaction collapse
            mod_key = 'impossible_key_to_match'

    # Get values for chosen side
    for y in (reaction_dictionary[neighbor][side_key] \
            + (
                [x[0]
                for x in reaction_dictionary[neighbor]['modifiers']
                    if x[1] == mod_key])):
        for z in graph.nodes()[y]['values']:
            values.append(z)

    # Return the only the names of reactants or products (no modifiers) if
    # there are at least one measured value for the reaction
    if len(values) > 0 \
            and any([False if v is None else True for v in values]):
        eval_items = reaction_dictionary[neighbor][side_key]
    else:
        eval_items = []

    return eval_items

def add_collapsed_components(
        graph,
        rxn,
        ref,
        samples):
    """Add new collapsed reaction node along with all
    needed edges to related components
    """

    reaction_color = (0.75, 0.75, 0.75, 1) # Grey
    colors = [reaction_color for x in range(samples)]

    # Add basic collapsed reaction info
    graph.add_node(rxn)
    graph.nodes()[rxn]['id'] = rxn
    graph.nodes()[rxn]['name'] = ref[rxn]['name']
    graph.nodes()[rxn]['reversible'] = 'false'
    graph.nodes()[rxn]['notes'] = ref[rxn]['notes']
    graph.nodes()[rxn]['type'] = 'collapsed'
    graph.nodes()[rxn]['sub_type'] = 'collapsed'

    # Add reaction node default values for color
    graph.nodes()[rxn]['values'] = [
        None for x in range(samples)]
    graph.nodes()[rxn]['values_rgba'] = colors
    graph.nodes()[rxn]['values_js'] = convert_rgba(
        rgba_tuples=colors)
    graph.nodes()[rxn]['stats'] = [
        None for x in range(samples)]
    graph.nodes()[rxn]['stats_rgba'] = colors
    graph.nodes()[rxn]['stats_js'] = convert_rgba(
        rgba_tuples=colors)

    for x in ref[rxn]['reactants'] + ref[rxn]['products']:
        graph.add_edges_from([
            (x, rxn)])
        graph.edges()[(x, rxn)]['type'] = 'collapsed'
        graph.edges()[(x, rxn)]['sub_type'] = 'collapsed'

    for x in ref[rxn]['modifiers']:
        graph.add_edges_from([
            (x[0], rxn)])
        graph.edges()[(x[0], rxn)]['type'] = 'collapsed_' + x[1]
        graph.edges()[(x[0], rxn)]['sub_type'] = 'collapsed_' + x[1]

    return graph

def collapse_nodes(
        graph,
        reaction_dictionary,
        samples,
        collapse_with_modifiers):
    """After values are broadcast, collapse network by creating new reaction
    dictionary
    Methods:
    1) If at least one value is on both side of the reaction, keep as is
        - inhibitors' values are counted in the reactants (more inhibitor,
        more build-up of input)
        - catalysts' values are counted in the products (more catalyst, more
        product)
    2) If one side of the reaction has a value, but the other doesn't, search
    for neighbors with matching reactants/products for whatever is missing. If
    the other side of any of those matching neighbors has a value, collapse.
    3) If both sides have no values, search neighbors for both sides. As in #2,
    if both have a distal value, collapse the three reactions together.
    """

    updated_reactions = {} # Collapsed reaction dictionary for plotting
    changed_reactions = {} # For mapping collapsed reactions post-processing
    # This makes sure that collapsed reactions are not added twice
    collapsed_options = []
    removed_reaction = set() # for reactions that are collapsed, make sure the
    # original reactions are removed from the final reaction dictionary

    for rxn in list(reaction_dictionary.keys()):

        # Parse out reaction metadata
        key = rxn
        compartment = reaction_dictionary[rxn]['compartment']
        id = reaction_dictionary[rxn]['id']
        name = reaction_dictionary[rxn]['name']
        reversible = reaction_dictionary[rxn]['reversible']
        notes = reaction_dictionary[rxn]['notes']
        reactants = reaction_dictionary[rxn]['reactants']
        products = reaction_dictionary[rxn]['products']
        modifiers = reaction_dictionary[rxn]['modifiers']
        additional_components = (
            reaction_dictionary[rxn]['additional_components'])

        # Get effective components, where inputs and outputs include modifiers
        # Get real components for reaction end matching
        if collapse_with_modifiers == "True" or collapse_with_modifiers == True:
            real_reactants = reaction_dictionary[rxn]['reactants']
            effective_reactants = (
                reaction_dictionary[rxn]['reactants'] \
                + (
                    [x[0]
                    for x in reaction_dictionary[rxn]['modifiers']
                        if x[1] == 'inhibitor']
                )
            )

            real_products = reaction_dictionary[rxn]['products']
            effective_products = (
                reaction_dictionary[rxn]['products'] \
                + (
                    [x[0]
                    for x in reaction_dictionary[rxn]['modifiers']
                        if x[1] == 'catalyst']
                )
            )
            real_modifiers = reaction_dictionary[rxn]['modifiers']
        else:
            real_reactants = effective_reactants = (
                reaction_dictionary[rxn]['reactants']
            )
            real_products = effective_products = (
                reaction_dictionary[rxn]['products']
            )
            real_modifiers = reaction_dictionary[rxn]['modifiers']

        # Collect values for effective inputs and outputs
        inputs = []
        for r in effective_reactants:
            for x in graph.nodes()[r]['values']:
                inputs.append(x)
        reactants_exist = any([False if x is None else True for x in inputs])

        outputs = []
        for p in effective_products:
            for y in graph.nodes()[p]['values']:
                outputs.append(y)
        products_exist = any([False if y is None else True for y in outputs])

        # If inputs and outputs both have at least one value, push to new dict
        # as is
        if reactants_exist and products_exist:
            changed_reactions[(key)] = key
            updated_reactions[key] = {
                'collapsed': 'false',
                'collapsed_reactions': [],
                'compartment': compartment,
                'id': id,
                'name': name,
                'reversible': reversible,
                'notes': notes,
                'reactants': reactants,
                'products': products,
                'modifiers': modifiers,
                'additional_components': additional_components
            }

        # Check other conditions where reactions may be able to be collapsed
        else:
            # Get the matching neighbors of this reaction for next decision tree
            input_neighbors = []
            output_neighbors = []

            # Check for reactions with matching sides
            for neighbor_key in reaction_dictionary.keys():

                # Parse potential bridge inputs and outputs
                neighbor_reactants = (
                    reaction_dictionary[neighbor_key]['reactants'])
                neighbor_products = (
                    reaction_dictionary[neighbor_key]['products'])
                neighbor_modifiers = (
                    reaction_dictionary[neighbor_key]['modifiers'])

                # If the current reaction's neighbors inputs or outputs match,
                # append that reaction to the current reaction's neighbors list

                ### Testing: Making sure the reaction inputs and outputs are not identical
                if (real_reactants == neighbor_reactants \
                        or real_reactants == neighbor_products) \
                        and real_reactants + real_products != \
                            neighbor_reactants + neighbor_products \
                        and real_modifiers != neighbor_modifiers \
                        and neighbor_key != key:
                    input_neighbors.append(neighbor_key)
                if (real_products == neighbor_reactants \
                        or real_products == neighbor_products) \
                        and real_reactants + real_products != \
                            neighbor_reactants + neighbor_products \
                        and real_modifiers != neighbor_modifiers \
                        and neighbor_key != key:
                    output_neighbors.append(neighbor_key)

            # Run one-sided bridging for reactions where inputs exist and
            # outputs have neighbors (could be with the neighbor's reactants
            # or products)
            if reactants_exist and len(output_neighbors) != 0:

                # For all output reaction neighbor options...
                for o in output_neighbors:
                    # Find side that doesn't match the current reaction outputs
                    if (set(real_products) \
                            == set(reaction_dictionary[o]['reactants'])):
                        side_key = 'products'
                        if collapse_with_modifiers == True or collapse_with_modifiers == "True":
                            mod_key = 'catalyst'
                        else:
                            mod_key = 'impossible_key_to_match'
                    else:
                        side_key = 'reactants'
                        if collapse_with_modifiers == True or collapse_with_modifiers == "True":
                            mod_key = 'inhibitor'
                        else:
                            mod_key = 'impossible_key_to_match'

                    # Check if values exist for the new side
                    # Allow modifier values to count
                    outputs = []
                    for oo in (reaction_dictionary[o][side_key] \
                            + (
                                [x[0]
                                for x in reaction_dictionary[o]['modifiers']
                                    if x[1] == mod_key]
                            )):
                        for y in graph.nodes()[oo]['values']:
                            outputs.append(y)

                    # If a neighbor has output values, create compressed
                    # reaction
                    if any([False if y is None else True for y in outputs]):

                        if set([rxn, o]) not in collapsed_options:
                            collapsed_options.append(set([rxn, o]))
                            add = id + '_' + reaction_dictionary[o]['id']
                            changed_reactions[(rxn, o)] = add
                            removed_reaction.add(id)
                            removed_reaction.add(reaction_dictionary[o]['id'])
                            updated_reactions[add] = {
                                'collapsed': 'true',
                                'collapsed_reactions': [rxn, o],
                                'compartment': compartment,
                                'id': id + '_' + reaction_dictionary[o]['id'],
                                'name': name \
                                    + ' // ' + reaction_dictionary[o]['name'],
                                'reversible': 'false',
                                'notes': 'Compressed reaction between ' \
                                    + name \
                                    + ' // ' + reaction_dictionary[o]['name'],
                                'reactants': reactants,
                                'products': reaction_dictionary[o]['products'],
                                'modifiers': modifiers \
                                    + reaction_dictionary[o]['modifiers'],
                                'additional_components': additional_components \
                                    + reaction_dictionary[o][
                                        'additional_components']
                            }
                            graph = add_collapsed_components(
                                graph=graph,
                                rxn=add,
                                ref=updated_reactions,
                                samples=samples)

                    # Create basic reaction
                    else:
                        changed_reactions[(key)] = key
                        updated_reactions[key] = {
                            'collapsed': 'false',
                            'collapsed_reactions': [],
                            'compartment': compartment,
                            'id': id,
                            'name': name,
                            'reversible': reversible,
                            'notes': notes,
                            'reactants': reactants,
                            'products': products,
                            'modifiers': modifiers,
                            'additional_components': additional_components
                        }

            # Run one-sided bridging for reactions where outputs exist and
            # inputs have neighbors
            elif products_exist and len(input_neighbors) != 0:

                # For all input reaction neighbor options...
                for i in input_neighbors:

                    # Find side that doesn't match the current reaction outputs
                    if (set(reactants) \
                            == set(reaction_dictionary[i]['reactants'])):
                        side_key = 'products'
                        if collapse_with_modifiers == True or collapse_with_modifiers == "True":
                            mod_key = 'catalyst'
                        else:
                            mod_key = 'impossible_key_to_match'
                    else:
                        side_key = 'reactants'
                        if collapse_with_modifiers == True or collapse_with_modifiers == "True":
                            mod_key = 'inhibitor'
                        else:
                            mod_key = 'impossible_key_to_match'

                    # Check if values exist for the new side
                    # Allow modifier values to count
                    inputs = []
                    for ii in (reaction_dictionary[i][side_key] \
                        + (
                            [x[0]
                            for x in reaction_dictionary[i]['modifiers']
                                if x[1] == mod_key]
                        )):
                        for z in graph.nodes()[ii]['values']:
                            inputs.append(z)

                    # If a neighbor has output values, create compressed
                    if any([False if z is None else True for z in inputs]):

                        if set([rxn, i]) not in collapsed_options:
                            collapsed_options.append(set([rxn, i]))
                            add = id + '_' + reaction_dictionary[i]['id']
                            changed_reactions[(rxn, i)] = add
                            removed_reaction.add(id)
                            removed_reaction.add(reaction_dictionary[i]['id'])
                            updated_reactions[add] = {
                                'collapsed': 'true',
                                'collapsed_reactions': [rxn, i],
                                'compartment': compartment,
                                'id': id + '_' \
                                    + reaction_dictionary[i]['id'],
                                'name': name \
                                    + ' // ' + reaction_dictionary[i]['name'],
                                'reversible': 'false',
                                'notes': 'Compressed reaction between ' \
                                    + name \
                                    + ' and ' + reaction_dictionary[i]['name'],
                                'reactants': \
                                    reaction_dictionary[i]['reactants'],
                                'products': products,
                                'modifiers': modifiers \
                                    + reaction_dictionary[i]['modifiers'],
                                'additional_components': additional_components \
                                    + reaction_dictionary[i]['additional_components']
                            }
                            graph = add_collapsed_components(
                                graph=graph,
                                rxn=add,
                                ref=updated_reactions,
                                samples=samples)

                    # Create basic reaction
                    else:
                        changed_reactions[(key)] = key
                        updated_reactions[key] = {
                            'collapsed': 'false',
                            'collapsed_reactions': [],
                            'compartment': compartment,
                            'id': id,
                            'name': name,
                            'reversible': reversible,
                            'notes': notes,
                            'reactants': reactants,
                            'products': products,
                            'modifiers': modifiers,
                            'additional_components': additional_components
                        }

            # Handle cases where both sides are missing values by searching to
            # see both neighbors have values that can be used to fill in
            else:

                # If both neighbors can connect, do so
                if len(input_neighbors) != 0 and len(output_neighbors) != 0:

                    # Cycle through all possible combinations of input and
                    # output neighbors and get their components
                    for i in input_neighbors:
                        for j in output_neighbors:

                            eval_i = find_values(
                                graph=graph,
                                reaction_dictionary=reaction_dictionary,
                                neighbor=i,
                                current_inputs=reactants,
                                collapse_with_modifiers=collapse_with_modifiers)
                            eval_j = find_values(
                                graph=graph,
                                reaction_dictionary=reaction_dictionary,
                                neighbor=j,
                                current_inputs=products,
                                collapse_with_modifiers=collapse_with_modifiers)

                            # If both neighbor reactions have values and the
                            # the collapsed reaction is not already in the
                            # modified # reaction dictionary, add collapsed
                            # reaction
                            if len(eval_i) > 0 and len(eval_j) > 0 \
                                    and (
                                    set([rxn, i, j]) not in collapsed_options):

                                collapsed_options.append(set([rxn, i, j]))
                                add = i + '_' + rxn + '_' + j
                                changed_reactions[(rxn, i, j)] = add
                                removed_reaction.add(rxn)
                                removed_reaction.add(i)
                                removed_reaction.add(j)
                                updated_reactions[add] = {
                                    'collapsed': 'true',
                                    'collapsed_reactions': [rxn, i, j],
                                    'compartment': compartment,
                                    'id': i + '_' + rxn + '_' + j,
                                    'name': reaction_dictionary[i]['name'] \
                                        + ' // ' + name + ' // ' \
                                        + reaction_dictionary[j]['name'],
                                    'reversible': 'false',
                                    'notes': 'Compressed reaction between ' \
                                        + reaction_dictionary[i]['name'] \
                                        + ' and ' + name + ' and ' \
                                        + reaction_dictionary[j]['name'],
                                    'reactants': eval_i,
                                    'products': eval_j,
                                    'modifiers': modifiers \
                                        + reaction_dictionary[i]['modifiers'] \
                                        + reaction_dictionary[j]['modifiers'],
                                    'additional_components': \
                                        additional_components \
                                        + reaction_dictionary[i]['additional_components']\
                                        + reaction_dictionary[j]['additional_components']
                                }
                                graph = add_collapsed_components(
                                    graph=graph,
                                    rxn=add,
                                    ref=updated_reactions,
                                    samples=samples)

                            else:
                                changed_reactions[(key)] = key
                                updated_reactions[key] = {
                                    'collapsed': 'false',
                                    'collapsed_reactions': [],
                                    'compartment': compartment,
                                    'id': id,
                                    'name': name,
                                    'reversible': reversible,
                                    'notes': notes,
                                    'reactants': reactants,
                                    'products': products,
                                    'modifiers': modifiers,
                                    'additional_components': \
                                        additional_components
                                }
                else:
                    changed_reactions[(key)] = key
                    updated_reactions[key] = {
                        'collapsed': 'false',
                        'collapsed_reactions': [],
                        'compartment': compartment,
                        'id': id,
                        'name': name,
                        'reversible': reversible,
                        'notes': notes,
                        'reactants': reactants,
                        'products': products,
                        'modifiers': modifiers,
                        'additional_components': additional_components
                    }

    for k in removed_reaction:
        if k in updated_reactions:
            del updated_reactions[k]

    return graph, updated_reactions, changed_reactions, removed_reaction
