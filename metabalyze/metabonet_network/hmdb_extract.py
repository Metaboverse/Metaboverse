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
import copy
import pickle
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from metabalyze.metabonet_network.utils import write_file_table
from metabalyze.metabonet_network.utils import collect_unique_elements
from metabalyze.metabonet_network.utils import confirm_path_directory
from metabalyze.utils import progress_counter

"""Set globals
"""
event_start_ns = 'start-ns'
event_end = 'end'
event_types = (
    'start',
    event_end,
    event_start_ns,
    'end-ns')

metabolite_tag = 'metabolite'
record_id = 'identifier'
accession_tag = 'accession'
secondary_accessions_tag = 'secondary_accessions'
name_tag = 'name'
synonyms_tag = 'synonyms'
synonym_tag = 'synonym'
pubchem_tag = 'pubchem_compound_id'
chebi_tag = 'chebi_id'
kegg_tag = 'kegg_id'

"""Make HMDB tag
"""
def make_tag(
        tag,
        name_space):

    return '{' + name_space + '}' + tag


"""Extracts information about metabolites from Human Metabolome Database
(HMDB).
arguments:
    hmdb (str): path to file of HMDB in Extensible Markup Language (XML)
returns:
    (dict<dict>): information from HMDB
"""
# Can't find event_start_ns or namespace ever, not grabbing elements and parsing
# Try using Cameron's method that has hanging variables
def extract_hmdb_summary(
        hmdb_file):

    # Collect references to name space
    spaces = {}

    # Collect information about metabolites
    summary_hmdb = {}

    # Count records
    count = 0

    for record, hmdb_element in et.iterparse(
            hmdb_file,
            events=event_types):
        print(record)
        # Populate name space dictionary
        if record == event_start_ns:
            name_space = hmdb_element[1]
        else:
            name_space = 'NameSpaceNotFound'

        # Populate summary dictionary
        if record == event_end:
            element_tag = make_tag(
                tag=metabolite_tag,
                name_space=name_space)

            if hmdb_element.tag == element_tag:

                # Extract information from record.
                extracted_record = extract_hmdb_record_summary(
                    hmdb_element=hmdb_element,
                    name_space=name_space)
                summary_hmdb[extracted_record[record_id]] = extracted_record

                # Clear memory and count
                element.clear()
                count += 1
                progress_counter(
                    count,
                    status='metabolites extracted')


    # Report.
    print('Extraction complete for ' + str(count) + ' metabolites.')

    return summary_hmdb

"""Extracts information about a metabolite from Human Metabolome Database
(HMDB).
arguments:
    element (object): element within XML tree
    space (str): name of specific name space within XML document
    spaces (dict<str>): name spaces within XML document
returns:
    (dict<str>): information about a metabolite from HMDB
"""
def extract_hmdb_record_summary(
        hmdb_element,
        name_space):

    # HMDB identifiers.
    hmdb_primary = extract_subelement_value(
        element=hmdb_element,
        tag=accession_tag,
        name_space=name_space)
    hmdb_secondary = extract_subelement(
        element=hmdb_element,
        tag=secondary_accessions_tag,
        name_space=name_space)
    references_hmdb_values = extract_subelement_values(
        element=hmdb_secondary,
        tag=accession_tag,
        name_space=name_space)
    references_hmdb_values.append(hmdb_primary)
    references_hmdb = collect_unique_elements(references_hmdb_values)

    # Name.
    name = extract_subelement_value(
        element=hmdb_element,
        tag=name_tag,
        name_space=name_space)

    # Synonyms.
    synonyms_element = extract_subelement(
        element=hmdb_element,
        tag=synonyms_tag,
        name_space=name_space)
    synonyms_values = extract_subelement_values(
        element=synonyms_element,
        tag=synonym_tag,
        name_space=name_space)
    synonyms_values.append(name)
    synonyms = collect_unique_elements(synonyms_values)

    # References.
    pubchem_tentative = extract_subelement_value(
        element=hmdb_element,
        tag=pubchem_tag,
        name_space=name_space)

    # Multiple entries have references to identifier '0' for PubChem.
    # This identifier is nonsense and erroneous.
    if pubchem_tentative is not None \
    and pubchem_tentative == '0':
        reference_pubchem = None

    else:
        reference_pubchem = pubchem_tentative

    reference_chebi = extract_subelement_value(
        element=hmdb_element,
        tag=chebi_tag,
        name_space=name_space)
    reference_kegg = extract_subelement_value(
        element=hmdb_element,
        tag=kegg_tag,
        name_space=name_space)

    # Compile and return information.
    record = {
        'identifier': hmdb_primary,
        'name': name,
        'synonyms': synonyms,
        'references_hmdb': references_hmdb,
        'reference_pubchem': reference_pubchem,
        'reference_chebi': reference_chebi,
        'reference_kegg': reference_kegg}

    return record

"""Extracts a reference to an element within another element in an XML tree.
arguments:
    element (object): element within XML tree
    tag (str): name of element's tag
    space (str): name of specific name space within XML document
    spaces (dict<str>): name spaces within XML document
returns:
    (object): element within XML tree
"""
def extract_subelement(
        element,
        tag,
        name_space):

    tag_complete = make_tag(
        tag=tag,
        name_space=name_space)

    return element.find(tag_complete)

"""Extracts the text content of an element within another element in an XML
tree.
arguments:
    element (object): element within XML tree
    tag (str): name of element's tag
    space (str): name of specific name space within XML document
    spaces (dict<str>): name spaces within XML document
returns:
    (str): text content of element
"""
def extract_subelement_value(
        element,
        tag,
        name_space):

    subelement = extract_subelement(
        element,
        tag,
        name_space)

    if subelement is not None:
        return subelement.text

    else:
        return None

"""Extracts the text contents of multiple elements within another element in
an XML tree.
arguments:
    element (object): element within XML tree
    tag (str): name of element's tag
    space (str): name of specific name space within XML document
    spaces (dict<str>): name spaces within XML document
returns:
    (list<str>): text contents of elements
"""
def extract_subelement_values(
        element,
        tag,
        name_space):

    tag_complete = make_tag(
        tag=tag,
        name_space=name_space)

    values = []

    for subelement in element.findall(tag_complete):

        values.append(subelement.text)

    return values

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
    path_pickle = directory + 'hmdb_summary.pickle'
    path_text = directory + 'hmdb_summary.tsv'

    # Write information to file
    with open(path_pickle, 'wb') as file_product:
        pickle.dump(information['summary_object'], file_product)

    write_file_table(
        information=information['summary_list'],
        path_file=path_text,
        names=information['summary_list'][0].keys(),
        delimiter='\t')

"""Function to execute module's main behavior.
The purpose of this procedure is to extract relevant information from the
Human Metabolome Database.
arguments:
    directory (str): path to directory for source and product files
"""
def __main__(
        args_dict):

    # Extract information from Human Metabolome Database.
    summary_hmdb = extract_hmdb_summary(
        hmdb_file=args_dict['hmdb'])

    # Compile information.
    information = {
        'summary_object': summary_hmdb,
        'summary_list': list(summary_hmdb.values())}

    #Write product information to file
    write_product(
        directory=args_dict['extract'],
        information=information)
