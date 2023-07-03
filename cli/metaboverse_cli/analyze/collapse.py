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

"""Import dependencies
"""
try:
    from utils import track_progress
    from analyze.utils import convert_rgba
except:
    import os
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath(os.path.join(".", "metaboverse_cli", "utils.py")))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    track_progress = utils.track_progress

    spec = importlib.util.spec_from_file_location(
        "convert_rgba", os.path.abspath("./metaboverse_cli/analyze/utils.py"))
    convert_rgba = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(convert_rgba)
    convert_rgba = convert_rgba.convert_rgba


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
        collapse_with_modifiers,
        degree_dictionary,
        collapse_threshold,
        degree_threshold):
    """Are there any values for either side of the reaction?
    - There can only be one side with values for this situation to be valid in
    this context
    """

    eval_items = []  # What will be returned
    values = []

    # Check if values exist for the new side
    # Allow modifier values to count
    pass_match = get_instructions(
        side1=current_inputs,
        side2=reaction_dictionary[neighbor]['reactants'],
        degree_dictionary=degree_dictionary,
        collapse_with_modifiers=collapse_with_modifiers,
        collapse_threshold=collapse_threshold,
        degree_threshold=degree_threshold)
    if pass_match == True:
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
    for y in (reaction_dictionary[neighbor][side_key]
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

    reaction_color = (0.75, 0.75, 0.75, 1)  # Grey
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


def get_reaction_values(
        collapse_with_modifiers,
        reaction):
    """Get real and effect inputs and outputs lists
    """

    # Get effective components, where inputs and outputs include modifiers
    # Get real components for reaction end matching
    real_reactants = reaction['reactants']
    real_products = reaction['products']
    real_modifiers = reaction['modifiers']

    if collapse_with_modifiers == "True" \
    or collapse_with_modifiers == True:
        effective_reactants = (
            reaction['reactants']
            + (
                [x[0]
                 for x in reaction['modifiers']
                    if x[1] == 'inhibitor']
            )
        )
        effective_products = (
            reaction['products']
            + (
                [x[0]
                 for x in reaction['modifiers']
                    if x[1] == 'catalyst']
            )
        )
    else:
        effective_reactants = real_reactants
        effective_products = real_products

    return (
        real_reactants,
        effective_reactants,
        real_products,
        effective_products,
        real_modifiers)


def get_unity_lists(
        reaction_list,
        neighbor_list):

    unity_len = len(list(set(reaction_list).intersection(neighbor_list)))
    _max = max([len(reaction_list), len(neighbor_list)])
    if _max <= 0:
        _max = 1
    _proportion = unity_len / _max

    return _proportion


def get_instructions(
        side1,
        side2,
        degree_dictionary,
        collapse_with_modifiers,
        collapse_threshold,
        degree_threshold):
    """See which sides were matching for collapse
    """

    # Find side that doesn't match the current reaction outputs
    short_side1 = []
    short_side2 = []
    for r in side1:
        if degree_dictionary[r] <= degree_threshold:
            short_side1.append(r)
    for p in side2:
        if degree_dictionary[p] <= degree_threshold:
            short_side2.append(p)
    unity_proportion = get_unity_lists(
        reaction_list=short_side1,
        neighbor_list=short_side2)

    if unity_proportion >= collapse_threshold:
        pass_match = True
    else:
        pass_match = False

    return pass_match


def check_neighbors(
        key,
        real_reactants,
        real_products,
        real_modifiers,
        neighbor_key,
        neighbor,
        degree_dictionary,
        input_neighbors,
        output_neighbors,
        blocklist,
        collapse_threshold,
        degree_threshold):

    # Parse potential bridge inputs and outputs
    neighbor_reactants = neighbor['reactants']
    neighbor_products = neighbor['products']
    neighbor_modifiers = neighbor['modifiers']

    # Account for no modifiers in matching
    if len(real_modifiers) == 0:
        real_modifiers = 'impossible1'
    if len(neighbor_modifiers) == 0:
        neighbor_modifiers = 'impossible2'

    # If the current reaction's neighbors inputs or outputs form a
    # complete match, append that reaction to the current
    # reaction's neighbors list
    if (real_reactants == neighbor_reactants
        or real_reactants == neighbor_products) \
            and real_reactants + real_products != \
            neighbor_reactants + neighbor_products \
            and real_modifiers != neighbor_modifiers \
            and neighbor_key != key:
        input_neighbors.append(neighbor_key)
    if (real_products == neighbor_reactants
        or real_products == neighbor_products) \
            and real_reactants + real_products != \
            neighbor_reactants + neighbor_products \
            and real_modifiers != neighbor_modifiers \
            and neighbor_key != key:
        output_neighbors.append(neighbor_key)

    # Check for partial matches with hubs removed from consideration
    short_reactants = []
    short_products = []
    short_neighbor_reactants = []
    short_neighbor_products = []

    for r in real_reactants:
        if degree_dictionary[r] <= degree_threshold and r not in blocklist:
            short_reactants.append(r)
    for p in real_products:
        if degree_dictionary[p] <= degree_threshold and p not in blocklist:
            short_products.append(p)
    for rr in neighbor_reactants:
        if degree_dictionary[rr] <= degree_threshold and rr not in blocklist:
            short_neighbor_reactants.append(rr)
    for pp in neighbor_products:
        if degree_dictionary[pp] <= degree_threshold and pp not in blocklist:
            short_neighbor_products.append(pp)

    unity_reactants_nnReactants = get_unity_lists(
        reaction_list=short_reactants,
        neighbor_list=short_neighbor_reactants)
    unity_reactants_nnProducts = get_unity_lists(
        reaction_list=short_reactants,
        neighbor_list=short_neighbor_products)
    unity_products_nnReactants = get_unity_lists(
        reaction_list=short_products,
        neighbor_list=short_neighbor_reactants)
    unity_products_nnProducts = get_unity_lists(
        reaction_list=short_products,
        neighbor_list=short_neighbor_products)

    # Check if short lists match and the length of the unity / shortest short
    # list meets threshold
    if short_reactants + short_products != \
            short_neighbor_reactants + short_neighbor_products \
    and real_modifiers != neighbor_modifiers \
    and neighbor_key != key \
    and unity_reactants_nnReactants >= collapse_threshold:
        input_neighbors.append(neighbor_key)
    if short_reactants + short_products != \
            short_neighbor_reactants + short_neighbor_products \
    and real_modifiers != neighbor_modifiers \
    and neighbor_key != key \
    and unity_reactants_nnProducts >= collapse_threshold:
        input_neighbors.append(neighbor_key)

    if short_reactants + short_products != \
            short_neighbor_reactants + short_neighbor_products \
    and real_modifiers != neighbor_modifiers \
    and neighbor_key != key \
    and unity_products_nnReactants >= collapse_threshold:
        output_neighbors.append(neighbor_key)
    if short_reactants + short_products != \
            short_neighbor_reactants + short_neighbor_products \
    and real_modifiers != neighbor_modifiers \
    and neighbor_key != key \
    and unity_products_nnProducts >= collapse_threshold:
        output_neighbors.append(neighbor_key)

    # Remove duplicates
    input_neighbors_unique = list(set(input_neighbors))
    output_neighbors_unique = list(set(output_neighbors))

    return input_neighbors, output_neighbors


def collapse_nodes(
        args_dict,
        graph,
        reaction_dictionary,
        neighbors_dictionary,
        degree_dictionary,
        samples,
        collapse_with_modifiers,
        blocklist,
        degree_threshold=50,
        collapse_threshold=0.3):
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

    print('Collapse match threshold: ' + str(collapse_threshold))
    print('Degree threshold for collapsing: ' + str(degree_threshold))

    updated_reactions = {}  # Collapsed reaction dictionary for plotting
    changed_reactions = {}  # For mapping collapsed reactions post-processing
    # This makes sure that collapsed reactions are not added twice
    collapsed_options = []
    removed_reaction = set()  # for reactions that are collapsed, make sure the
    # original reactions are removed from the final reaction dictionary

    counter = 0
    reaction_number = len(list(reaction_dictionary.keys()))

    print('Analyzing ' + str(reaction_number) + ' reactions for collapsing...')

    for rxn in list(reaction_dictionary.keys()):
        counter = track_progress(args_dict, counter, reaction_number, 10)

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

        # Collect values
        real_reactants, effective_reactants, real_products, \
            effective_products, real_modifiers = get_reaction_values(
                reaction=reaction_dictionary[rxn],
                collapse_with_modifiers=collapse_with_modifiers)

        # Collect values for effective inputs and outputs
        real_inputs = []
        for r in real_reactants:
            if degree_dictionary[r] <= degree_threshold and r not in blocklist:
                for x in graph.nodes()[r]['values']:
                    real_inputs.append(x)
        real_reactants_exist = any(
            [False if x is None else True for x in real_inputs])

        effective_inputs = []
        for rr in effective_reactants:
            if degree_dictionary[rr] <= degree_threshold and rr not in blocklist:
                for xx in graph.nodes()[rr]['values']:
                    effective_inputs.append(xx)
        effective_reactants_exist = any(
            [False if xx is None else True for xx in effective_inputs])

        real_outputs = []
        for p in real_products:
            if degree_dictionary[p] <= degree_threshold and p not in blocklist:
                for y in graph.nodes()[p]['values']:
                    real_outputs.append(y)
        real_products_exist = any(
            [False if y is None else True for y in real_outputs])

        effective_outputs = []
        for pp in effective_products:
            if degree_dictionary[pp] <= degree_threshold and pp not in blocklist:
                for yy in graph.nodes()[pp]['values']:
                    effective_outputs.append(yy)
        effective_products_exist = any(
            [False if yy is None else True for yy in effective_outputs])

        # If inputs and outputs both have at least one value, push to new dict
        # as is
        if real_reactants_exist and real_products_exist:
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

            # Check for reactions with complete and partial matching sides
            if key in neighbors_dictionary.keys():
                # all components connected to that reaction
                for neighbor_key in neighbors_dictionary[key]:
                    if key != neighbor_key \
                    and neighbor_key in reaction_dictionary.keys():
                        input_neighbors, output_neighbors = check_neighbors(
                            key=key,
                            real_reactants=real_reactants,
                            real_products=real_products,
                            real_modifiers=real_modifiers,
                            neighbor_key=neighbor_key,
                            neighbor=reaction_dictionary[neighbor_key],
                            degree_dictionary=degree_dictionary,
                            input_neighbors=input_neighbors,
                            output_neighbors=output_neighbors,
                            blocklist=blocklist,
                            collapse_threshold=collapse_threshold,
                            degree_threshold=degree_threshold)
            else:
                for neighbor_key in reaction_dictionary.keys():
                    if key != neighbor_key:
                        input_neighbors, output_neighbors = check_neighbors(
                            key=key,
                            real_reactants=real_reactants,
                            real_products=real_products,
                            real_modifiers=real_modifiers,
                            neighbor_key=neighbor_key,
                            neighbor=reaction_dictionary[neighbor_key],
                            degree_dictionary=degree_dictionary,
                            input_neighbors=input_neighbors,
                            output_neighbors=output_neighbors,
                            blocklist=blocklist,
                            collapse_threshold=collapse_threshold,
                            degree_threshold=degree_threshold)

            # Run one-sided bridging for reactions where inputs exist and
            # outputs have neighbors (could be with the neighbor's reactants
            # or products)
            if effective_reactants_exist and len(output_neighbors) != 0:

                # For all output reaction neighbor options...
                for o in output_neighbors:
                    pass_match = get_instructions(
                        side1=real_products,
                        side2=reaction_dictionary[o]['reactants'],
                        degree_dictionary=degree_dictionary,
                        collapse_with_modifiers=collapse_with_modifiers,
                        collapse_threshold=collapse_threshold,
                        degree_threshold=degree_threshold)
                    if pass_match == True:
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
                    for oo in (reaction_dictionary[o][side_key]
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
                                'name': name
                                + ' // ' + reaction_dictionary[o]['name'],
                                'reversible': 'false',
                                'notes': 'Compressed reaction between '
                                + name
                                + ' // ' + reaction_dictionary[o]['name'],
                                'reactants': reactants,
                                'products': reaction_dictionary[o]['products'],
                                'modifiers': modifiers
                                + reaction_dictionary[o]['modifiers'],
                                'additional_components': additional_components
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
            elif effective_products_exist and len(input_neighbors) != 0:

                # For all input reaction neighbor options...
                for i in input_neighbors:
                    # Find side that doesn't match the current reaction outputs
                    pass_match = get_instructions(
                        side1=reactants,
                        side2=reaction_dictionary[i]['reactants'],
                        degree_dictionary=degree_dictionary,
                        collapse_with_modifiers=collapse_with_modifiers,
                        collapse_threshold=collapse_threshold,
                        degree_threshold=degree_threshold)
                    if pass_match == True:
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
                    for ii in (reaction_dictionary[i][side_key]
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
                                'id': id + '_'
                                + reaction_dictionary[i]['id'],
                                'name': name
                                + ' // ' + reaction_dictionary[i]['name'],
                                'reversible': 'false',
                                'notes': 'Compressed reaction between '
                                + name
                                + ' and ' + reaction_dictionary[i]['name'],
                                'reactants':
                                reaction_dictionary[i]['reactants'],
                                'products': products,
                                'modifiers': modifiers
                                + reaction_dictionary[i]['modifiers'],
                                'additional_components': additional_components
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
                                collapse_with_modifiers=collapse_with_modifiers,
                                degree_dictionary=degree_dictionary,
                                collapse_threshold=collapse_threshold,
                                degree_threshold=degree_threshold)
                            eval_j = find_values(
                                graph=graph,
                                reaction_dictionary=reaction_dictionary,
                                neighbor=j,
                                current_inputs=products,
                                collapse_with_modifiers=collapse_with_modifiers,
                                degree_dictionary=degree_dictionary,
                                collapse_threshold=collapse_threshold,
                                degree_threshold=degree_threshold)

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
                                    'name': reaction_dictionary[i]['name']
                                    + ' // ' + name + ' // '
                                    + reaction_dictionary[j]['name'],
                                    'reversible': 'false',
                                    'notes': 'Compressed reaction between '
                                    + reaction_dictionary[i]['name']
                                    + ' and ' + name + ' and '
                                    + reaction_dictionary[j]['name'],
                                    'reactants': eval_i,
                                    'products': eval_j,
                                    'modifiers': modifiers
                                    + reaction_dictionary[i]['modifiers']
                                    + reaction_dictionary[j]['modifiers'],
                                    'additional_components':
                                    additional_components
                                        + reaction_dictionary[i]['additional_components']
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
                                    'additional_components':
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
