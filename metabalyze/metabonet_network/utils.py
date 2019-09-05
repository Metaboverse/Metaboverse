"""License Information
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
import csv
import copy
import re
from textwrap import dedent
import string

"""Import internal dependencies
"""
from metabalyze.utils import progress_bar

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

    if (os.path.exists(path)) \
    and (len(os.listdir(path)) < 1):
        os.rmdir(path)

"""Confirms that a path to a directory exists.
Creates a directory if it does not already exist.
arguments:
    path (str): path to directory
"""
def confirm_path_directory(
        path):

    if not os.path.exists(path):
        os.makedirs(path)

"""Confirms that a path to a directory exists.
Creates a directory if it does not already exist.
arguments:
    path (str): path to directory
"""
def confirm_file(
        file,
        type=''):

    if not os.path.isfile(file):
        raise Exception(str(type), 'can not be found')

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
    report = dedent("""\
        Curation report:
        --------------------------------------------------
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

"""Converts string of characters to lower case with only alphabetical or
numerical characters.
arguments:
    characters (str): characters in a string
returns:
    (str): characters in a string
"""
def convert_string_low_alpha_num(
        characters):

    # Convert all characters to lower case.
    characters_lower = characters.lower()

    # Remove all characters other than alphabetical or numerical characters.
    characters_novel = characters_lower

    for character in characters_lower:
        if ((character not in string.ascii_letters) \
        and (character not in string.digits)):
            characters_novel = characters_novel.replace(character, '')

    return characters_novel

"""Removes a directory if it is empty.
arguments:
    path (str): path to directory
"""
def remove_empty_directory(
        path):

    if (os.path.exists(path)) \
    and (len(os.listdir(path)) < 1):
        os.rmdir(path)

"""Reads and organizes source information from file
This function reads and organizes relevant information from file.
arguments:
    path_file (str): path to directory and file
    names (list<str>): names for values in each row of table
    delimiter (str): delimiter between values in the table
returns:
    (list<dict>): tabular information from file
"""
def read_file_table(
        path_file,
        names,
        delimiter):

    # Read information from file
    with open(path_file, 'r') as file_source:

        reader = csv.DictReader(
            file_source,
            fieldnames=names,
            delimiter=delimiter)
        table = list(map(lambda row: dict(row), list(reader)))

    return table

"""Writes information to file
arguments:
    information (list<str>): information
    path_file (str): path to directory and file
    names (list<str>): names for values in each row of table
    delimiter (str): delimiter between values in the table
"""
def write_file_table(
        information,
        path_file,
        names,
        delimiter):

    # Write information to file
    #with open(out_file_path_model, "w") as out_file:
    #    out_file.write(content_identifier)
    with open(path_file, "w") as file_product:

        writer = csv.DictWriter(
            file_product,
            fieldnames=names,
            delimiter=delimiter)
        writer.writeheader()
        writer.writerows(information)

"""Finds the first element in a sequence to match a condition, otherwise none
arguments:
    match (function): condition for elements to match
    sequence (list): sequence of elements
returns:
    (object | NoneType): first element from sequence to match condition or
        none
"""
def find_match(
        match=None,
        sequence=None):

    for element in sequence:

        if match(element):
            return element

    return None

"""Finds index of first element in sequence to match a condition, otherwise -1
arguments:
    match (function): condition for elements to match
    sequence (list): sequence of elements
returns:
    (int): index of element if it exists or -1 if it does not exist
"""
def find_index(
        match=None,
        sequence=None):

    for index, element in enumerate(sequence):

        if match(element):
            # Element matches condition
            # Return element's index
            return index

    # Not any elements match condition
    # Return -1
    return -1

"""Finds all elements in a sequence to match a condition, otherwise none.
arguments:
    match (function): condition for elements to match
    sequence (list): sequence of elements
returns:
    (list<dict> | NoneType): elements from sequence to match condition or
        none
"""
def find_all(
        match=None,
        sequence=None):

    matches = []

    for element in sequence:

        if match(element):
            matches.append(element)

    if len(matches) > 0:
        return matches

    else:

        return None

"""Collects unique elements
arguments:
    elements_original (list): sequence of elements
returns:
    (list): unique elements
"""
def collect_unique_elements(
        elements_original=None):

    elements_novel = []

    for element in elements_original:

        if element not in elements_novel:
            elements_novel.append(element)

    return elements_novel

"""Collects a single value from multiple records
arguments:
    key (str): key of value in each record
    records (list<dict>): sequence of records
returns:
    (list): values from records
"""
def collect_value_from_records(
        key=None,
        records=None):

    def access(record):

        return record[key]

    return list(map(access, records))

"""Collects values from multiple records.
arguments:
    key (str): key of value in each record
    records (list<dict>): sequence of records
returns:
    (list): values from records
"""
def collect_values_from_records(
        key=None,
        records=None):

    collection = []

    for record in records:

        collection.extend(record[key])

    return collection

"""Compares lists by inclusion
arguments:
    list_one (list): list of elements
    list_two (list): list of elements
returns:
    (bool): whether first list includes all elements from second
"""
def compare_lists_by_inclusion(
        list_one=None,
        list_two=None):

    def match(
            element_two=None):

        return element_two in list_one

    matches = list(map(match, list_two))

    return all(matches)

"""Compares lists by mutual inclusion
arguments:
    list_one (list): list of elements
    list_two (list): list of elements
returns:
    (bool): whether each list includes all elements from the other
"""
def compare_lists_by_mutual_inclusion(
        list_one=None,
        list_two=None):

    forward = compare_lists_by_inclusion(
        list_one=list_one,
        list_two=list_two)
    reverse = compare_lists_by_inclusion(
        list_one=list_two,
        list_two=list_one)

    return forward and reverse

"""
Filters elements by whether both of two lists include them
arguments:
    list_one (list): list of elements
    list_two (list): list of elements
returns:
    (list): elements that both of two lists include
"""
def filter_common_elements(
        list_one=None,
        list_two=None):

    def match(
            element_two=None):

        return element_two in list_one

    return list(filter(match, list_two))

"""Collects values of a target attribute that occur together in records with
each value of another category attribute.
Each record has a single value of the target attribute.
Each record can have either a single value or multiple values of the
category attribute.
These collections do not necessarily include only unique values of the
target attribute.
arguments:
    target (str): name of attribute in records to collect for each category
    category (str): name of attribute in records to define categories
    records (list<dict>): records with target and category attributes
returns:
    (dict<list<str>>): values of the target attribute that occur together
        in records with each value of the category attribute
"""
def collect_records_targets_by_categories(
        target=None,
        category=None,
        records=None):

    def collect_record_target_by_category(
            target_value=None,
            category_value=None,
            collection_original=None):

        collection_novel = copy.deepcopy(collection_original)

        # Determine whether collection includes the category's value.
        if category_value in collection_novel.keys():
            # Collection includes the category's value.
            target_values = collection_novel[category_value]
            target_values.append(target_value)

            # Include target's value in collection.
            collection_novel[category_value] = target_values

        else:
            # Collection does not include the category's value.
            # Include category's value and target's value in collection.
            collection_novel[category_value] = [target_value]

        return collection_novel

    collection = {}

    counter = 1
    total = len(records)

    for record in records:

        target_value = record[target]
        category_values = record[category]

        if isinstance(category_values, list):

            for category_value in category_values:

                collection = collect_record_target_by_category(
                    target_value=target_value,
                    category_value=category_value,
                    collection_original=collection)

        else:
            category_value = category_values
            collection = collect_record_target_by_category(
                target_value=target_value,
                category_value=category_value,
                collection_original=collection)

        progress_bar(
            counter,
            total,
            status='Collecting records')
        counter += 1

    return collection

"""Collects a single value from a specific record in a reference.
arguments:
    key (str): key of value in record
    identifiers (list<str>): identifiers of records in reference
    reference (dict<dict<str>>): reference of records
returns:
    (list<str>): values from records
"""
def collect_values_from_records_in_reference(
        key=None,
        identifiers=None,
        reference=None):

    values = []

    for identifier in identifiers:

        record = reference[identifier]
        value = record[key]
        values.append(value)

    return values

"""Filters nonempty elements.
arguments:
    elements_original (list<str>): sequence of elements
returns:
    (list<str>): non-empty elements
"""
def filter_nonempty_elements(
        elements_original=None):

    elements_novel = []

    for element in elements_original:

        if len(str(element)) > 0:
            elements_novel.append(element)

    return elements_novel

"""Filters nodes and links by identifiers.
arguments:
    identifiers (list<str>): identifiers of elements to keep
    entries_original (dict<dict>): entries
returns:
    (dict<dict>): entries
"""
def filter_entries_identifiers(
        identifiers=None,
        entries_original=None):

    entries_novel = {}

    for entry in entries_original.values():

        if entry['identifier'] in identifiers:
            entries_novel[entry['identifier']] = entry

    return entries_novel


# Human Metabolome Database (HMDB).
"""Matches entries from Human Metabolome Database by identifiers or names.
arguments:
    identifiers (list<str>): identifiers by which to find entries in HMDB
    names (list<str>): names by which to find entries in HMDB
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (list<str>): keys of entries in HMDB
"""
def match_hmdb_entries_by_identifiers_names(
        identifiers=None,
        names=None,
        summary_hmdb=None):

    # Ensure identifiers and names are not empty.
    identifiers_valid = filter_nonempty_elements(identifiers)
    names_valid = filter_nonempty_elements(names)

    # Determine whether measurement's record include reference to HMDB.
    if (len(identifiers_valid) > 0):
        # Measurement's record includes references to HMDB.
        # Match measurement's record to a entries in HMDB.
        # Match by identifier.
        hmdb_keys = filter_hmdb_entries_by_identifiers(
            identifiers=identifiers_valid,
            summary_hmdb=summary_hmdb)

    elif (len(names_valid) > 0):
        # Measurement's record does not include reference to HMDB.
        # Match measurement's record to an entry in HMDB.
        # Attempt to match by name.
        hmdb_keys = filter_hmdb_entries_by_synonyms(
            names=names_valid,
            summary_hmdb=summary_hmdb)

    else:
        hmdb_keys = []

    return hmdb_keys

"""Filters entries from HMDB by their identifiers.
arguments:
    identifiers (list<str>): identifiers by which to find entries in HMDB
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (list<str>): keys of entries in HMDB
"""
def filter_hmdb_entries_by_identifiers(
        identifiers=None,
        summary_hmdb=None):

    keys = []

    for key, record in summary_hmdb.items():

        hmdb_entry_identifiers = record['references_hmdb']

        # Determine whether any of entry's identifiers match the metabolite's
        # references
        checks = []

        for identifier in identifiers:

            check = identifier in hmdb_entry_identifiers
            checks.append(check)

        if any(checks):
            # The entry matches the metabolite's references.
            keys.append(key)

    return keys

"""Filters entries from HMDB by their synonyms.
arguments:
    names (list<str>): names by which to find entries in HMDB
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (list<str>): keys of entries in HMDB
"""
def filter_hmdb_entries_by_synonyms(
        names=None,
        summary_hmdb=None):

    keys = []

    for key, record in summary_hmdb.items():

        synonyms = record['synonyms']
        synonyms_comparison = []

        for synonym in synonyms:

            synonym_comparison = convert_string_low_alpha_num(synonym)
            synonyms_comparison.append(synonym_comparison)

        # Determine whether any of entry's identifiers match the metabolite's
        # references
        checks = []

        for name in names:

            name_comparison = convert_string_low_alpha_num(name)
            check = name_comparison in synonyms_comparison
            checks.append(check)

        if any(checks):
            # The entry matches the metabolite's references
            keys.append(key)

    return keys

"""Filters entries from HMDB by their identifiers.
arguments:
    reference (str): name of reference
    identifiers (list<str>): identifiers by which to find entries in HMDB
    summary_hmdb (dict<dict>): information about metabolites from Human
        Metabolome Database (HMDB)
returns:
    (list<str>): keys of entries in HMDB
"""
def filter_hmdb_entries_by_references(
        reference=None,
        identifiers=None,
        summary_hmdb=None):

    keys = []

    for key, record in summary_hmdb.items():

        reference_entry = record[reference]

        # Determine whether any of entry's references match the query.
        match = reference_entry in identifiers

        if match:
            # The entry matches the metabolite's references.
            keys.append(key)

    return keys

# Metabolic information.
"""Counts entities with any of specific references.
arguments:
    references (list<str>): identifiers of references
    entities (dict<dict>): information about entities
returns:
    (int): count of entities with specific reference
"""
def count_entities_with_references(
        references=None,
        entities=None):

    count = 0

    for entity in entities.values():

        matches = []

        for reference in references:

            if reference in entity['references'].keys():

                if len(entity['references'][reference]) > 0:
                    matches.append(True)

        if any(matches):
            count += 1

    return count

"""Collects a value from a reaction's specific participants
arguments:
    key (str): key of value to collect from each participant
    criteria (dict<list>): criteria by which to select participants
    participants (list<dict>): information about a reaction's participants
returns:
    (list<str>): values from a reaction's participants
"""
def collect_reaction_participants_value(
        key,
        criteria,
        participants):

    participants_match = filter_reaction_participants(
        criteria=criteria,
        participants=participants)

    return collect_value_from_records(
        key=key,
        records=participants_match)

"""Filters a reaction's participants by multiple criteria
arguments:
    criteria (dict<list>): criteria by which to select participants
    participants (list<dict>): information about a reaction's participants
returns:
    (list<dict>): information about a reaction's participants
"""
def filter_reaction_participants(
        criteria,
        participants):

    def match(
            participant):

        if "metabolites" in criteria:
            match_metabolite = (participant["metabolite"] in criteria["metabolites"])

        else:
            match_metabolite = True

        if "compartments" in criteria:
            match_compartment = (participant["compartment"] in criteria["compartments"])

        else:
            match_compartment = True

        if "roles" in criteria:
            match_role = participant["role"] in criteria["roles"]

        else:
            match_role = True

        return match_metabolite and match_compartment and match_role

    return list(filter(match, participants))

"""Copies and interprets content from Recon
This function copies and interprets content from a metabolic model in
Systems Biology Markup Language (SBML), a form of Extensible Markup
Language (XML).
returns:
    (object): references to definition of name space and sections within
        content
Previously called: copy_interpret_content_recon
"""
def get_recon_references(
        reference,
        range_limit=10):

    # Copy content
    reference_copy = copy.deepcopy(reference)

    # Init content definitions dictionary
    recon_references = {}
    recon_references['content'] = reference

    # Set references to overall model
    recon_references['model'] = reference_copy.getroot()[0]

    # Set references to data types
    for x in recon_references['model']:

        if 'listOfCompartments' in x.tag:
            recon_references['compartments'] = x

        elif 'listOfSpecies' in x.tag:
            recon_references['metabolites'] = x

        elif 'listOfReactions' in x.tag:
            recon_references['reactions'] = x

        else:
            pass

    # Check that references were filled
    if not 'compartments' in recon_references.keys():
        raise Exception('Compartments information not found in Recon XML file')

    if not 'metabolites' in recon_references.keys():
        raise Exception('Metabolites information not found in Recon XML file')

    if not 'reactions' in recon_references.keys():
        raise Exception('Reactions information not found in Recon XML file')

    # Get smbl version tag
    version = ''
    version = re.search('{(.*)}', recon_references['model'].tag).group(1)

    # Get syntax tag
    syntax = ''
    test_group = recon_references['compartments']

    for x in range(0,range_limit):

        for y in range(0,range_limit):

            for z in range(0,range_limit):

                try:
                    t = test_group[x][y][z]
                    t_find = re.search('{(.*)}', t.tag).group(1)

                    if 'syntax' in t_find:
                        syntax = t_find
                        break

                    else:
                        pass

                except:
                    pass
            else:
                continue # Continue if the inner loop wasn't broken
        else:
            continue # Continue if the inner loop wasn't broken

        break # Inner loop was broken, break the outer.

    if version != '' and syntax != '':
        recon_references['space'] = {
            'version': version,
            'syntax': syntax}

    else:
        raise Exception('Could not find version or syntax information from Recon XML file')

    return recon_references
