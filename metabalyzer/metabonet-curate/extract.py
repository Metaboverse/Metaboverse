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
    alias: metabalyzer
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

###############################################################################
# Notes

# The purpose of this procedure is to enhance information about metabolic sets
# and entities.

###############################################################################
# Installation and importation of packages and modules


# Packages and modules from the python standard library

import os
#import sys
import shutil
#import importlib
import csv
import copy
import pickle
import xml.etree.ElementTree as et

# Packages and modules from third parties

#import numpy
#import pandas
#import scipy

# Packages and modules from local source

import metabonet.utility as utility

###############################################################################
# Functionality


def read_source(directory=None):
    """
    Reads and organizes source information from file

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_source = os.path.join(directory, "source")
    path_hmdb = os.path.join(path_source, "hmdb_metabolites.xml")
    # Read information from file.
    #hmdb = et.parse(path_hmdb)
    # Compile and return information.
    return {
        "hmdb": path_hmdb
    }


def extract_hmdb_summary(hmdb=None):
    """
    Extracts information about metabolites from Human Metabolome Database
    (HMDB).

    arguments:
        hmdb (str): path to file of HMDB in Extensible Markup Language (XML)

    raises:

    returns:
        (dict<dict>): information from HMDB

    """

    # Collect references to name space.
    spaces = {}
    # Collect information about metabolites.
    summary_hmdb = {}
    # Count records.
    count = 0
    for event, element in et.iterparse(
        hmdb, events=('start', 'end', 'start-ns', 'end-ns')
    ):
        if event == "start-ns":
            space = element[0]
            reference = element[1]
            spaces[space] = reference
        if event == "end":
            if element.tag == construct_tag(
                tag="metabolite", space=space, spaces=spaces
            ):
                # Parse complete for a new metabolite.
                # Count.
                count = count + 1
                # Extract information from record.
                record = extract_hmdb_record_summary(
                    element=element,
                    space=space,
                    spaces=spaces
                )
                summary_hmdb[record["identifier"]] = record
                # Clear memory.
                element.clear()
    # Report.
    print("Extraction complete for " + str(count) + " metabolites.")
    # Return information.
    return summary_hmdb


def construct_tag(tag=None, space=None, spaces=None):
    """
    Constructs complete tag for name space in Extensible Markup Language (XML).

    arguments:
        tag (str): name of element's tag
        space (str): name of specific name space within XML document
        spaces (dict<str>): name spaces within XML document

    raises:

    returns:
        (str): complete name of tag

    """

    return "{" + spaces[space] + "}" + tag


def extract_hmdb_record_summary(element=None, space=None, spaces=None):
    """
    Extracts information about a metabolite from Human Metabolome Database
    (HMDB).

    arguments:
        element (object): element within XML tree
        space (str): name of specific name space within XML document
        spaces (dict<str>): name spaces within XML document

    raises:

    returns:
        (dict<str>): information about a metabolite from HMDB

    """

    # HMDB identifiers.
    hmdb_primary = extract_subelement_value(
        element=element,
        tag="accession",
        space=space,
        spaces=spaces
    )
    hmdb_secondary = extract_subelement(
        element=element,
        tag="secondary_accessions",
        space=space,
        spaces=spaces
    )
    references_hmdb_values = extract_subelement_values(
        element=hmdb_secondary,
        tag="accession",
        space=space,
        spaces=spaces
    )
    references_hmdb_values.append(hmdb_primary)
    references_hmdb = utility.collect_unique_elements(references_hmdb_values)
    # Name.
    name = extract_subelement_value(
        element=element,
        tag="name",
        space=space,
        spaces=spaces
    )
    # Synonyms.
    synonyms_element = extract_subelement(
        element=element,
        tag="synonyms",
        space=space,
        spaces=spaces
    )
    synonyms_values = extract_subelement_values(
        element=synonyms_element,
        tag="synonym",
        space=space,
        spaces=spaces
    )
    synonyms_values.append(name)
    synonyms = utility.collect_unique_elements(synonyms_values)
    # References.
    pubchem_tentative = extract_subelement_value(
        element=element,
        tag="pubchem_compound_id",
        space=space,
        spaces=spaces
    )
    # Multiple entries have references to identifier "0" for PubChem.
    # This identifier is nonsense and erroneous.
    if (
        (pubchem_tentative is not None) and
        (pubchem_tentative == "0")
    ):
        reference_pubchem = None
    else:
        reference_pubchem = pubchem_tentative
    reference_chebi = extract_subelement_value(
        element=element,
        tag="chebi_id",
        space=space,
        spaces=spaces
    )
    reference_kegg = extract_subelement_value(
        element=element,
        tag="kegg_id",
        space=space,
        spaces=spaces
    )
    # Compile and return information.
    record = {
        "identifier": hmdb_primary,
        "name": name,
        "synonyms": synonyms,
        "references_hmdb": references_hmdb,
        "reference_pubchem": reference_pubchem,
        "reference_chebi": reference_chebi,
        "reference_kegg": reference_kegg
    }
    return record


def extract_subelement(element=None, tag=None, space=None, spaces=None):
    """
    Extracts a reference to an element within another element in an XML tree.

    arguments:
        element (object): element within XML tree
        tag (str): name of element's tag
        space (str): name of specific name space within XML document
        spaces (dict<str>): name spaces within XML document

    raises:

    returns:
        (object): element within XML tree

    """

    tag_complete = construct_tag(
        tag=tag, space=space, spaces=spaces
    )
    return element.find(tag_complete)


def extract_subelement_value(element=None, tag=None, space=None, spaces=None):
    """
    Extracts the text content of an element within another element in an XML
    tree.

    arguments:
        element (object): element within XML tree
        tag (str): name of element's tag
        space (str): name of specific name space within XML document
        spaces (dict<str>): name spaces within XML document

    raises:

    returns:
        (str): text content of element

    """

    #name = element.find("{http://www.hmdb.ca}name").text

    subelement = extract_subelement(
        element=element,
        tag=tag,
        space=space,
        spaces=spaces
    )
    if subelement is not None:
        return subelement.text
    else:
        return None


def extract_subelement_values(element=None, tag=None, space=None, spaces=None):
    """
    Extracts the text contents of multiple elements within another element in
    an XML tree.

    arguments:
        element (object): element within XML tree
        tag (str): name of element's tag
        space (str): name of specific name space within XML document
        spaces (dict<str>): name spaces within XML document

    raises:

    returns:
        (list<str>): text contents of elements

    """

    tag_complete = construct_tag(
        tag=tag, space=space, spaces=spaces
    )
    values = []
    for subelement in element.findall(tag_complete):
        values.append(subelement.text)
    return values


def write_product(directory=None, information=None):
    """
    Writes product information to file

    arguments:
        directory (str): directory for product files
        information (object): information to write to file

    raises:

    returns:

    """

    # Specify directories and files.
    path = os.path.join(directory, "extraction")
    utility.confirm_path_directory(path)
    path_pickle = os.path.join(path, "hmdb_summary.pickle")
    path_text = os.path.join(path, "hmdb_summary.tsv")
    # Write information to file.
    with open(path_pickle, "wb") as file_product:
        pickle.dump(information["summary_object"], file_product)
    utility.write_file_table(
        information=information["summary_list"],
        path_file=path_text,
        names=information["summary_list"][0].keys(),
        delimiter="\t"
    )


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to extract relevant information from the
    Human Metabolome Database.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)
    # Extract information from Human Metabolome Database.
    summary_hmdb = extract_hmdb_summary(hmdb=source["hmdb"])
    # Compile information.
    information = {
        "summary_object": summary_hmdb,
        "summary_list": list(summary_hmdb.values())
    }
    #Write product information to file
    write_product(directory=directory, information=information)
