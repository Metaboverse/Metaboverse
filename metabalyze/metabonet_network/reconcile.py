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
import shutil
import re
import csv
import xml.etree.ElementTree as et
import copy

"""Import internal dependencies
"""
from metabalyze.metabonet_network.utils import confirm_file
from metabalyze.metabonet_network.utils import read_file_table
from metabalyze.metabonet_network.utils import confirm_path_directory
from metabalyze.metabonet_network.collect import copy_interpret_content_recon2m2

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'

"""Reads and organizes source information from file
arguments:
    directory (str): directory of source files
returns:
    (object): source information
"""
def read_source(
        recon_xml):

    # Specify customization directories and files
    path_customization = __path__ + 'customization'
    path_compartments = os.path.join(path_customization, 'reconciliation_compartments.tsv')
    path_metabolites = os.path.join(path_customization, 'reconciliation_metabolites.tsv')

    # Read in recon database
    confirm_file(recon_xml)
    recon_content = et.parse(recon_xml)

    curation_compartments = read_file_table(
        path_file=path_compartments,
        names=None,
        delimiter='\t')
    curation_metabolites = read_file_table(
        path_file=path_metabolites,
        names=None,
        delimiter='\t')

    return {
        'content': recon_content,
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
        content):

    # Copy and interpret content
    reference = copy_interpret_content_recon(
        content=content)

    # Correct designation of model's boundary in metabolites
    for metabolite in reference['metabolites'].findall(
            'version:species',
            reference['space']):

        # Determine whether metabolite's compartment is model's boundary
        if 'boundary' in metabolite.attrib['id']:

            novel_identifier = re.sub(r'_[eciglmnrx]_boundary', '_b', metabolite.attrib['id'])
            metabolite.attrib['id'] = novel_identifier
            novel_compartment = 'b'
            metabolite.attrib['compartment'] = novel_compartment

    # Correct designation of model's boundary in reactions.
    for reaction in reference['reactions'].findall(
            'version:reaction', reference['space']):

        # Search reaction's metabolites.
        """
        Remove this hard-coding
        """
        for metabolite in reaction.iter(
                '{http://www.sbml.org/sbml/level2/version4}speciesReference'):

            # Determine whether metabolite's compartment is model's boundary.
            if 'boundary' in metabolite.attrib['species']:

                novel_identifier = re.sub(
                    r'_[eciglmnrx]_boundary',
                    '_b',
                    metabolite.attrib['species'])
                metabolite.attrib['species'] = novel_identifier

    return reference['content']

"""Counts compartments, reactions, and metabolites in model.
arguments:
    content (object): content from file in Systems Biology Markup Language
        (XML)
returns:
    (dict<int>): summary
"""
def count_model_sets_entities(
        content=None):

    # Copy and interpret content.
    reference = copy_interpret_content_recon2m2(
        content=content)

    # Count compartments.
    compartments = 0

    for compartment in reference['compartments'].findall(
            'version:compartment', reference['space']):

        compartments = compartments + 1

    # Count reactions.
    reactions = 0

    for reaction in reference['reactions'].findall(
            'version:reaction', reference['space']):

        reactions = reactions + 1

    # Count metabolites.
    metabolites = 0

    for metabolite in reference['metabolites'].findall(
            'version:species', reference['space']):

        metabolites = metabolites + 1

    return {
        'compartments': compartments,
        'reactions': reactions,
        'metabolites': metabolites}

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
        content=None):

    # Copy and interpret content.
    reference = copy_interpret_content_recon2m2(
        content=content)

    # Change content for each combination of original and novel information.
    for row in curation_compartments:

        # Detmerine whether to change compartment's name.
        if row['description_original'] != row['description_novel']:

            # Change information in compartments.
            for compartment in reference['compartments'].findall(
                    'version:compartment', reference['space']):

                if compartment.attrib['id'] == row['identifier_original']:

                    compartment.attrib['name'] = row['description_novel']

        # Determine whether to change compartment's identifier.
        if row['identifier_original'] != row['identifier_novel']:

            # Construct targets to recognize original and novel identifiers.
            # Use underscore prefix to match complete identifiers.
            original_elements = ['_', row['identifier_original']]
            original_target = ''.join(original_elements)
            novel_elements = ['_', row['identifier_novel']]
            novel_target = ''.join(novel_elements)

            # Change information in metabolites.
            # Change identifiers of metabolites.
            for metabolite in reference['metabolites'].findall(
                    'version:species', reference['space']):

                # Determine whether to change metabolite's identifier.
                if original_target in metabolite.attrib['id']:

                    metabolite.attrib['id'] = metabolite.attrib['id'].replace(
                        original_target,
                        novel_target)

            # Change information in reactions' metabolites.
            # Change identifiers of reactions' metabolites.
            for reaction in reference['reactions'].findall(
                    'version:reaction', reference['space']):

                # Search reaction's metabolites.
                """
                Remove hard-coding
                """
                for metabolite in reaction.iter(
                    '{http://www.sbml.org/sbml/level2/version4}'
                    'speciesReference'):

                    # Determine whether to change metabolite's identifier.
                    if original_target in metabolite.attrib['species']:

                        metabolite.attrib['species'] = (
                            metabolite.attrib['species'].replace(
                                original_target,
                                novel_target))

    return reference['content']

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
        content=None):

    # Copy and interpret content.
    reference = copy_interpret_content_recon2m2(
        content=content)

    # Remove prefixes from identifiers for metabolites.
    for metabolite in reference['metabolites'].findall(
            'version:species', reference['space']):

        # Remove prefix from metabolite's identifier.
        novel_identifier = re.sub(r'^M_', '', metabolite.attrib['id'])
        metabolite.attrib['id'] = novel_identifier

        # Search metabolite's annotation.
        for description in metabolite.iter(
                '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}description'):

            # Remove prefix from metabolite's identifier.
            novel_identifier = re.sub(
                r'^#M_',
                '#',
                description.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about'])
            description.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about'] = novel_identifier

    # Remove prefixes from identifiers for reactions' metabolites.
    for reaction in reference['reactions'].findall(
            'version:reaction', reference['space']):

        # Search reaction's metabolites.
        for metabolite in reaction.iter(
                '{http://www.sbml.org/sbml/level2/version4}speciesReference'):

            # Remove prefix from metabolite's identifier.
            novel_identifier = re.sub(r'^M_', '', metabolite.attrib['species'])
            metabolite.attrib['species'] = novel_identifier

    return reference['content']

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
        content=None):

    # Copy and interpret content.
    reference = copy_interpret_content_recon2m2(
        content=content)

    # Change content for each combination of original and novel identifiers.
    for row in curation_metabolites:

        # Construct targets to recognize original and novel identifiers.
        # Use trailing underscore to match complete identifiers.
        original_elements = [row['identifier_original'], '_']
        original_target = ''.join(original_elements)
        novel_elements = [row['identifier_novel'], '_']
        novel_target = ''.join(novel_elements)

        # Change identifiers of metabolites.
        for metabolite in reference['metabolites'].findall(
                'version:species',
                reference['space']):

            # Determine whether to change metabolite's identifier.
            if original_target in metabolite.attrib['id']:
                metabolite.attrib['id'] = metabolite.attrib['id'].replace(
                    original_target,
                    novel_target)

        # Change identifiers of reactions' metabolites.
        for reaction in reference['reactions'].findall(
                'version:reaction',
                reference['space']):

            # Search reaction's metabolites.
            for metabolite in reaction.iter(
                    '{http://www.sbml.org/sbml/level2/version4}speciesReference'):

                # Determine whether to change metabolite's identifier.
                if original_target in metabolite.attrib['species']:

                    metabolite.attrib['species'] = (
                        metabolite.attrib['species'].replace(
                            original_target,
                            novel_target))

    return reference['content']

"""Writes product information to file
arguments:
    directory (str): directory for product files
    information (object): information to write to file
"""
def write_product(
        directory,
        information=None):

    # Specify directories and files.
    confirm_path_directory(directory)
    path_file = os.path.join(
        directory,
        'recon2m2_reconciliation.xml')

    # Write information to file.
    information.write(
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

    # Change model's content and correct content where necessary
    content_boundary = change_model_boundary(
        content=recon['content'])
    content_compartments = change_model_compartments(
        curation_compartments=recon['curation_compartments'],
        content=content_boundary)
    content_prefix = remove_model_metabolite_prefix(
        content=content_compartments)
    content_metabolites = change_model_metabolites(
        curation_metabolites=recon['curation_metabolites'],
        content=content_prefix)

    #Write product information to file.
    write_product(
        args_dict['reconcile'],
        information=content_metabolites)

    # Summary.
    summary = count_model_sets_entities(
        content=content_metabolites)

    # Report.
    print('compartments: ' + str(summary['compartments']))
    print('reactions: ' + str(summary['reactions']))
    print('metabolites: ' + str(summary['metabolites']))
