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
import re
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from metabalyze.metabonet_network.utils import confirm_file
from metabalyze.metabonet_network.utils import read_file_table
from metabalyze.metabonet_network.utils import get_recon_references

"""Set globals
"""
__customization__  =  os.path.dirname(os.path.realpath(__file__)) + '/customization'

id_attribute = 'id'
species_attribute = 'species'
compartment_attribute = 'compartment'
name_attribute = 'name'
species_reference_id = 'speciesReference'
novel_metabolite_id = r'_[eciglmnrx]_boundary'
prefix_searcher = r'^M_'
prefix_hash_searcher = r'^#M_'
boundary_id = 'boundary'
species_finder = 'version:species'
reaction_finder = 'version:reaction'
compartment_finder = 'version:compartment'
description_original = 'description_original'
description_novel = 'description_novel'
identifier_original = 'identifier_original'
identifier_novel = 'identifier_novel'
annotation_searcher = 'description'
about_searcher = 'about'

reconcile_compartment = 'reconciliation_compartments.tsv'
reconcile_metabolites = 'reconciliation_metabolites.tsv'

"""Reads and organizes source information from file
arguments:
    directory (str): directory of source files
returns:
    (object): source information
"""
def read_source(
        recon_xml):

    # Specify customization directories and files
    path_compartments = __customization__ + reconcile_compartment
    path_metabolites = __customization__ + reconcile_metabolites

    # Read in recon database
    confirm_file(recon_xml)
    recon_reference = et.parse(recon_xml)

    curation_compartments = read_file_table(
        path_file=path_compartments,
        names=None,
        delimiter='\t')
    curation_metabolites = read_file_table(
        path_file=path_metabolites,
        names=None,
        delimiter='\t')

    return {
        'reference': recon_reference,
        'curation_compartments': curation_compartments,
        'curation_metabolites': curation_metabolites}

"""Changes annotations for a model's boundary
This function changes annotations of a model's boundary in compartments,
metabolites, and reactions.
arguments:
    content (object): content from file in Systems Biology Markup Language
        (XML)
returns:
    (object): content with changes
"""
def change_model_boundary(
        reference):

    # Correct designation of model's boundary in metabolites
    for metabolite in reference['metabolites'].findall(
            species_finder,
            reference['space']):

        # Determine whether metabolite's compartment is model's boundary
        if boundary_id in metabolite.attrib[id_attribute]:

            # Give novel identifier
            metabolite.attrib[id_attribute] = re.sub(
                novel_metabolite_id,
                '_b',
                metabolite.attrib[id_attribute])

            # Give novel compartment info
            metabolite.attrib[compartment_attribute] = 'b'

    # Correct designation of model's boundary in reactions.
    for reaction in reference['reactions'].findall(
            reaction_finder,
            reference['space']):

        # Search reaction's metabolites.
        iter_seacher = '{' + reference['space']['version'] + '}' + species_reference_id

        for metabolite in reaction.iter(iter_seacher):

            # Determine whether metabolite's compartment is model's boundary.
            if boundary_id in metabolite.attrib[species_attribute]:

                novel_identifier = re.sub(
                    novel_metabolite_id,
                    '_b',
                    metabolite.attrib[species_attribute])
                metabolite.attrib[species_attribute] = novel_identifier

    return reference

"""Changes annotations for a model's compartments
This function changes annotations of a model's compartments.
arguments:
    curation_compartments (list<dict<str>>): changes to information about
        compartments
    content (object): content from file in Systems Biology Markup Language
        (XML)
returns:
    (object): content with changes
"""
def change_model_compartments(
        curation_compartments=None,
        reference=None):

    # Change content for each combination of original and novel information
    for row in curation_compartments:

        # Detmerine whether to change compartment's name
        if row[description_original] != row[description_novel]:

            # Change information in compartments
            for compartment in reference['compartments'].findall(
                    compartment_finder,
                    reference['space']):

                if compartment.attrib[id_attribute] == row[identifier_original]:

                    compartment.attrib[name_attribute] = row[description_novel]

        # Determine whether to change compartment's identifier
        if row[identifier_original] != row[identifier_novel]:

            # Construct targets to recognize original and novel identifiers
            # Use underscore prefix to match complete identifiers
            original_target = ''.join(['_', row[identifier_original]])
            novel_target = ''.join(['_', row[identifier_novel]])

            # Change information in metabolites
            # Change identifiers of metabolites
            for metabolite in reference['metabolites'].findall(
                    species_finder,
                    reference['space']):

                # Determine whether to change metabolite's identifier
                if original_target in metabolite.attrib[id_attribute]:

                    metabolite.attrib[id_attribute] = metabolite.attrib[id_attribute].replace(
                        original_target,
                        novel_target)

            # Change information in reactions' metabolites
            # Change identifiers of reactions' metabolites
            for reaction in reference['reactions'].findall(
                    reaction_finder,
                    reference['space']):

                # Search reaction's metabolites
                for metabolite in reaction.iter(
                    '{' + reference['space']['version'] + '}',
                    species_reference_id):

                    # Determine whether to change metabolite's identifier.
                    if original_target in metabolite.attrib[species_attribute]:

                        metabolite.attrib[species_attribute] = metabolite.attrib[species_attribute].replace(
                                original_target,
                                novel_target)

    return reference

"""Removes unnecessary prefixes from identifiers for model's entities
This function removes unnecessary prefixes from identifiers for
metabolites.
arguments:
    content (object): content from file in Systems Biology Markup Language
        (XML)
returns:
    (object): content with changes
"""
def remove_model_metabolite_prefix(
        reference=None):

    # Remove prefixes from identifiers for metabolites.
    for metabolite in reference['metabolites'].findall(
            species_finder,
            reference['space']):

        # Remove prefix from metabolite's identifier.
        novel_identifier = re.sub(
            prefix_searcher,
            '',
            metabolite.attrib[id_attribute])
        metabolite.attrib[id_attribute] = novel_identifier

        # Search metabolite's annotation.
        for description in metabolite.iter(
                '{' + reference['space']['syntax'] + '}' + annotation_searcher):

            # Remove prefix from metabolite's identifier.
            novel_identifier = re.sub(
                prefix_hash_searcher,
                '#',
                description.attrib['{' + reference['space']['syntax'] + '}' + about_searcher])
            description.attrib['{' + reference['space']['syntax'] + '}' + about_searcher] = novel_identifier

    # Remove prefixes from identifiers for reactions' metabolites.
    for reaction in reference['reactions'].findall(
            reaction_finder,
            reference['space']):

        # Search reaction's metabolites.
        for metabolite in reaction.iter(
                '{' + reference['space']['version'] + '}' + species_reference_id):

            # Remove prefix from metabolite's identifier.
            novel_identifier = re.sub(
                prefix_searcher,
                '',
                metabolite.attrib[species_attribute])
            metabolite.attrib[species_attribute] = novel_identifier

    return reference

"""Changes metabolites' identifiers
This function changes metabolites' identifiers according to information
about translation.
arguments:
    curation_metabolites (list<dict<str>>): changes to information about
        metabolites
    content (object): content from file in Systems Biology Markup Language
        (XML)
returns:
    (object): content with changes
"""
def change_model_metabolites(
        curation_metabolites=None,
        reference=None):

    # Change content for each combination of original and novel identifiers.
    for row in curation_metabolites:

        # Construct targets to recognize original and novel identifiers.
        # Use trailing underscore to match complete identifiers.
        original_target = ''.join([row[identifier_original], '_'])
        novel_target = ''.join([row[identifier_novel], '_'])

        # Change identifiers of metabolites.
        for metabolite in reference['metabolites'].findall(
                species_finder,
                reference['space']):

            # Determine whether to change metabolite's identifier.
            if original_target in metabolite.attrib[id_attribute]:
                metabolite.attrib[id_attribute] = metabolite.attrib[id_attribute].replace(
                    original_target,
                    novel_target)

        # Change identifiers of reactions' metabolites.
        for reaction in reference['reactions'].findall(
                reaction_finder,
                reference['space']):

            # Search reaction's metabolites.
            for metabolite in reaction.iter(
                    '{' + reference['space']['version'] + '}' + species_reference_id):

                # Determine whether to change metabolite's identifier.
                if original_target in metabolite.attrib[species_attribute]:

                    metabolite.attrib[species_attribute] = (
                        metabolite.attrib[species_attribute].replace(
                            original_target,
                            novel_target))

    return reference

"""Counts compartments, reactions, and metabolites in model.
arguments:
    content (object): content from file in Systems Biology Markup Language
        (XML)
returns:
    (dict<int>): summary
"""
def count_model_sets_entities(
        reference=None):

    # Count compartments.
    compartments = 0

    for compartment in reference['compartments'].findall(
            compartment_finder,
            reference['space']):

        compartments = compartments + 1

    # Count reactions.
    reactions = 0

    for reaction in reference['reactions'].findall(
            reaction_finder,
            reference['space']):

        reactions = reactions + 1

    # Count metabolites.
    metabolites = 0

    for metabolite in reference['metabolites'].findall(
            species_finder,
            reference['space']):

        metabolites = metabolites + 1

    return {
        'compartments': compartments,
        'reactions': reactions,
        'metabolites': metabolites}

"""Writes product information to file
arguments:
    directory (str): directory for product files
    information (object): information to write to file
"""
def write_product(
        output=None,
        reference=None):

    # Specify directories and files
    path_file = output + 'recon_reconciled.xml'

    # Write information to file.
    reference['content'].write(
        path_file,
        xml_declaration=False)

"""Function to execute module's main behavior.
The purpose of this procedure is to reconcile the metabolic model to be
compatible with MetaNetX.
arguments:
    directory (str): path to directory for source and product files
"""
def __main__(
        args_dict):

    # Read source information from file
    recon = read_source(
        args_dict['recon'])

    # Copy and interpret content
    recon_reference = get_recon_references(
        recon['reference'])

    # Change model's content and correct content where necessary
    reference_updated_boundary = change_model_boundary(
        reference=recon_reference)
    reference_updated_compartments = change_model_compartments(
        curation_compartments=recon['curation_compartments'],
        reference=reference_updated_boundary)
    reference_prefix_removed = remove_model_metabolite_prefix(
        reference=reference_updated_compartments)
    reference_update_metabolites = change_model_metabolites(
        curation_metabolites=recon['curation_metabolites'],
        reference=reference_prefix_removed)

    #Write product information to file.
    write_product(
        output=args_dict['reconcile'],
        reference=reference_update_metabolites)

    # Summary.
    summary = count_model_sets_entities(
        reference=reference_update_metabolites)

    # Report.
    print('Metabolic model summary:')
    print('--------------------------------------------------')
    print('compartments: ' + str(summary['compartments']))
    print('reactions: ' + str(summary['reactions']))
    print('metabolites: ' + str(summary['metabolites']))
    print('--------------------------------------------------')
