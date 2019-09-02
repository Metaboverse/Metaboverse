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
import csv
import copy
import pickle
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from metabalyze.metabonet_network.utils import confirm_path_directory
from metabalyze.metabonet_network.utils import prepare_curation_report
from metabalyze.metabonet_network.utils import read_file_table
from metabalyze.metabonet_network.utils import collect_value_from_records
from metabalyze.metabonet_network.utils import search_source # originally called "find"

"""Read data
Reads and organizes source information from file
returns:
    (dict): source information
"""
def read_source(
        directory):

    # Specify directories and files.
    path_model = directory + 'recon2m2_reconciliation.xml'
    path_genes = directory + 'recon2m2_metanetx_genes.tsv'
    path_compartments = directory + 'recon2m2_metanetx_compartments.tsv'
    path_metabolites = directory + 'recon2m2_metanetx_metabolites.tsv'
    path_reactions = directory + 'recon2m2_metanetx_reactions.tsv'

    # Read information from file.
    model = et.parse(path_model)

    compartments = read_file_table(
        path_file=path_compartments,
        names=[
            'identifier',
            'name',
            'source'],
        delimiter='\t')

    genes = read_file_table(
        path_file=path_genes,
        names=[
            'reaction',
            'genes',
            'low_bound',
            'up_bound',
            'direction'],
        delimiter='\t')

    reactions = read_file_table(
        path_file=path_reactions,
        names=[
            'identifier',
            'equation',
            'recon2m2',
            'metanetx',
            'enzyme_commission',
            'processes',
            'references'],
        delimiter='\t')

    metabolites = read_file_table(
        path_file=path_metabolites,
        names=[
            'identifier',
            'name',
            'source',
            'formula',
            'mass',
            'charge',
            'references'],
        delimiter='\t')

    # Compile and return information.
    return {
        'model': model,
        'compartments': compartments,
        'genes': genes,
        'reactions': reactions,
        'metabolites': metabolites}

"""Copies and interprets content from Recon
This function copies and interprets content from a metabolic model in
Systems Biology Markup Language (SBML), a form of Extensible Markup
Language (XML).
returns:
    (object): references to definition of name space and sections within
        content
"""
def copy_interpret_content_recon(
        content=None):

    # Copy content
    content_copy = copy.deepcopy(content)

    # Get xml file smdb and xml encoding
    space = {
        "version": "http://www.sbml.org/sbml/level2/version4", # Get this automatically from file
        "syntax": "http://www.w3.org/1999/02/22-rdf-syntax-ns#" # Get this automatically from file
    }

    # Set references to sections within content
    model = content_copy.getroot()[0]
    compartments = model[1]
    metabolites = model[2]
    reactions = model[3]

    return {
        'space': space,
        'content': content,
        'model': model,
        'compartments': compartments,
        'metabolites': metabolites,
        'reactions': reactions}

"""Extracts information from source about compartments
arguments:
    compartments_source (list<dict>): source information about compartments
returns:
    (dict<dict>): information about compartments
"""
def extract_compartments(
        compartments_source=None):

    compartments = {}

    for compartment in compartments_source:

        record = {
            'identifier': compartment['identifier'],
            'name': compartment['name']}

        compartments[compartment['identifier']] = record

    return compartments

"""Extracts information from source about processes
arguments:
    reactions_source (list<dict>): source information about reactions
returns:
    (dict<dict>): information about processes
"""
def extract_processes(
        reactions_source=None):

    processes = {}

    for reaction in reactions_source:

        reaction_processes_names = extract_reaction_processes_names(
            reaction_source=reaction)

        for name in reaction_processes_names:

            # Determine whether a record exists for the process
            novelty = determine_process_name_novelty(
                name=name,
                processes=processes)

            if novelty:

                # Create and include a record for the process
                index = len(processes.keys())
                identifier = 'P' + str(index + 1)
                record = {
                    'identifier': identifier,
                    'name': name}
                processes[identifier] = record

    return processes

"""Extracts names of a reaction's metabolic processes
arguments:
    reaction_source (dict): source information about a reaction
returns:
    (list<str>): names of a reaction's process
"""
def extract_reaction_processes_names(
        reaction_source=None):

    # Separate references
    processes_source = reaction_source['processes']
    processes = extract_reference_information(
        key='model:',
        references_source=processes_source)

    return processes

"""Extracts reference information
arguments:
    key (str): identifier of reference information for a specific type
    references_source (str): source information about references
returns:
    (list<str>): reference information for a specific type
"""
def extract_reference_information(
        key=None,
        references_source=None):

    # Separate information for references
    references = references_source.split(';')

    # Filter identifiers for the reference
    pairs = list(filter(lambda pair: key in pair, references))

    # Remove key from identifiers
    identifiers = list(map(lambda pair: pair.replace(key, ''), pairs))

    return identifiers

"""Determines whether a process' name is novel
arguments:
    name (str): name of a process
    processes (dict<dict>): information about processes
returns:
    (bool): whether process' name is novel
"""
def determine_process_name_novelty(
        name=None,
        processes=None):

    for record in processes.values():

        if name == record['name']:

            return False

    return True

"""Extracts reactions' names from Recon
arguments:
    model (object): content from Recon 2M.2 in SBML
returns:
    (dict<str>): names of reactions
"""
def extract_reactions_names(
        model=None):

    # Copy and interpret content
    reference = copy_interpret_content_recon2m2(
        content=model)
    reactions_names = {}

    reactions = reference['reactions'].findall('version:reaction', reference['space'])

    for reaction in reactions:

        identifier = reaction.attrib['id']
        name = reaction.attrib['name']
        reactions_names[identifier] = name

    return reactions_names

"""Extracts information from source about reactions
arguments:
    reactions_source (list<dict>): source information about reactions
    reactions_names (dict<str>): names of reactions
    genes_source (list<dict>): source information about genes
    processes (dict<dict>): information about processes
returns:
    (dict<dict>): information about reactions
"""
def extract_reactions(
        reactions_source=None,
        reactions_names=None,
        genes_source=None,
        processes=None):

    reactions = {}

    for reaction_source in reactions_source:

        record = extract_reaction(
            reaction_source=reaction_source,
            reactions_names=reactions_names,
            genes_source=genes_source,
            processes=processes)

        reactions[record['identifier']] = record

    return reactions

"""Extracts information from source about a reaction
arguments:
    reaction_source (dict): source information about a reaction
    reactions_names (dict<str>): names of reactions
    genes_source (list<dict>): source information about genes
    processes (dict<dict>): information about processes
returns:
    (dict): information about a reaction
"""
def extract_reaction(
        reaction_source=None,
        reactions_names=None,
        genes_source=None,
        processes=None):

    # Determine information
    identifier = reaction_source['identifier']

    name = match_reaction_name(
        reaction_source=reaction_source,
        reactions_names=reactions_names)
    equation = reaction_source['equation']
    reversibility = extract_reaction_reversibility(
        equation=equation)
    participants = extract_reaction_participants(
        equation=equation)
    processes = extract_reaction_processes(
        reaction_source=reaction_source,
        processes=processes)
    references = extract_reaction_references(
        identifier=identifier,
        recon2m2=reaction_source['recon2m2'],
        metanetx=reaction_source['metanetx'],
        enzyme_commission=reaction_source['enzyme_commission'],
        genes=genes_source,
        references_source=reaction_source['references'])

    return {
        'identifier': identifier,
        'name': name,
        'equation': equation,
        'reversibility': reversibility,
        'participants': participants,
        'processes': processes,
        'references': references}

"""Enhances reactions by including names from Recon
arguments:
    reaction_source (dict): source information about a reaction
    reactions_names (dict<str>): names of reactions
returns:
    (dict<dict>): information about reactions
"""
def match_reaction_name(
        reaction_source=None,
        reactions_names=None):

    # Reconciliation merges some multiple reactions from Recon to a single
    # reaction in MetaNetX.
    identifiers = reaction_source["recon2m2"]
    identifiers_split = identifiers.split(";")
    count = len(identifiers_split)

    if count > 1:

        # Use first non-empty name.
        for index in range(count):

            identifier = identifiers_split[index]

            if identifier in reactions_names.keys():

                name_temporary = reactions_names[identifier]

            else:

                name_temporary = ''

            if len(name_temporary) > 1:

                name = name_temporary
                break

            if index == (count -1):

                name = name_temporary

    else:

        identifier = identifiers_split[0]

        if identifier in reactions_names.keys():

            name = reactions_names[identifier]

        else:

            name = ''

    return name

"""Extracts information about a reaction's reversibility
arguments:
    equation (str): a reaction's equation from MetaNetX
returns:
    (bool): whether reaction is reversible
"""
def extract_reaction_reversibility(
        equation=None):

    if '<==>' in equation:

        return True

    else:

        return False

"""Extracts information about a reaction's participants
arguments:
    equation (str): a reaction's equation from MetaNetX
returns:
    (list<dict>): information about a reaction's participants
"""
def extract_reaction_participants(
        equation=None):

    # Extract raw information about reaction's participants
    participants_raw = extract_reaction_equation_raw_participants_by_role(
        equation=equation)

    # Extract information about participants' role, coefficient, metabolite,
    # and compartment
    reactants = extract_reaction_participants_by_role(
        participants_raw=participants_raw['reactants'],
        role='reactant')
    products = extract_reaction_participants_by_role(
        participants_raw=participants_raw['products'],
        role='product')

    return reactants + products

"""Extracts raw information about a reaction's participants from its equation
arguments:
    equation (str): a reaction's equation from MetaNetX
returns:
    (dict<list<str>>): raw information about a reaction's participants by
        role
"""
def extract_reaction_equation_raw_participants_by_role(
        equation=None):

    # Separate information about participants' reactants from products
    # Determine reaction's directionality
    if '<==>' in equation:

        equation_sides = equation.split(' <==> ')
        reactants_side = equation_sides[0]
        products_side = equation_sides[1]

    elif '-->' in equation:

        equation_sides = equation.split(' --> ')
        reactants_side = equation_sides[0]
        products_side = equation_sides[1]

    elif '<--' in equation:

        equation_sides = equation.split(' <-- ')
        reactants_side = equation_sides[1]
        products_side = equation_sides[0]
    else:
        pass

    # Separate information about individual participants
    reactants = reactants_side.split(' + ')
    products = products_side.split(' + ')

    return {
        'reactants': reactants,
        'products': products}

"""Extracts information about a reaction's participants by role
arguments:
    participants_raw (list<str>): raw information about a reaction's
        participants
    role (str): participants' role in reaction
returns:
    (list<dict>): information about a reaction's participants
"""
def extract_reaction_participants_by_role(
        participants_raw=None,
        role=None):

    # Extract information about participants
    participants = []

    for participant_raw in participants_raw:

        record = extract_reaction_participant_by_role(
            participant_raw=participant_raw,
            role=role)
        participants.append(record)

    return participants

"""Extracts information about a reaction's participant by role
arguments:
    participant_raw (str): raw information about a reaction's participant
    role (str): participant's role in reaction
returns:
    (dict): information about a reaction's participant
"""
def extract_reaction_participant_by_role(
        participant_raw=None,
        role=None):

    # Separate information
    participant_split = participant_raw.split(' ')
    coefficient = float(participant_split[0])
    metabolite_compartment = participant_split[1].split('@')
    metabolite = metabolite_compartment[0]
    compartment = metabolite_compartment[1]

    return {
        'metabolite': metabolite,
        'compartment': compartment,
        'coefficient': coefficient,
        'role': role}

"""Extracts information about a reaction's genes
arguments:
    identifier (str): identifier of a reaction
    genes_source (list<dict>): source information about genes
returns:
    (list<str>): identifiers of a reaction's genes
"""
def extract_reaction_genes(
        identifier=None,
        genes_source=None):

    def match_reaction_gene(
            gene_record):

        return gene_record['reaction'] == identifier

    gene_source = search_source(
        match_reaction_gene,
        genes_source)
    genes_references = gene_source['genes']
    genes_split_one = genes_references.split(';')
    genes_split_two = genes_references.split('+')
    genes = []

    for gene_pair in genes_split_two:

        identifier = gene_pair.replace('gene:', '')

        if len(identifier) > 0:

            genes.append(identifier)

    return genes

"""Extracts identifiers of a reaction's metabolic processes
arguments:
    reaction_source (dict): source information about a reaction
    processes (dict<dict>): information about processes
returns:
    (list<str>): identifiers of a reaction's process
"""
def extract_reaction_processes(
        reaction_source=None,
        processes=None):

    # Find process
    def match_reaction_process(
            record):

        return record['name'] == name

    reaction_processes_names = extract_reaction_processes_names(
        reaction_source=reaction_source)

    reaction_processes = []

    for name in reaction_processes_names:

        process_record = search_source(
            match_reaction_process,
            processes.values())
        identifier = process_record['identifier']
        reaction_processes.append(identifier)

    return reaction_processes

"""Extracts references from source about a reaction
arguments:
    identifier (str): identifier of a reaction
    recon2m2 (str): identifier of reaction in original version of model,
        Recon 2M.2
    metanetx (str): identifier of reaction in MetaNetX
    enzyme_commission (str): identifier of reaction in Enzyme Commission
    genes (list<dict>): source information about genes
    references_source (str): source information about a reaction's
        references
returns:
    (dict<str>): references about a reaction
"""
def extract_reaction_references(
        identifier=None,
        recon2m2=None,
        metanetx=None,
        enzyme_commission=None,
        genes=None,
        references_source=None):

    # Collect identifiers for each reference.
    gene = extract_reaction_genes(
        identifier=identifier,
        genes_source=genes)

    if ';' in enzyme_commission:

        enzyme = enzyme_commission.split(';')

    else:

        enzyme = []

    # Extract information
    rhea = extract_reference_information(
        references_source=references_source,
        key='rhea:')
    bigg = extract_reference_information(
        references_source=references_source,
        key='bigg:')
    metanetx_prior = extract_reference_information(
        references_source=references_source,
        key='deprecated:')

    if len(metanetx) > 0:

        metanetx_current = [metanetx] + metanetx_prior

    else:

        metanetx_current = metanetx_prior

    # Extract KEGG info
    kegg = extract_reference_information(
        references_source=references_source,
        key='kegg:')
    metacyc = extract_reference_information(
        references_source=references_source,
        key='metacyc:')
    reactome = extract_reference_information(
        references_source=references_source,
        key='reactome:')
    sabiork = extract_reference_information(
        references_source=references_source,
        key='sabiork:')
    seed = extract_reference_information(
        references_source=references_source,
        key='seed:')

    return {
        'recon2m2': [recon2m2],
        'metanetx': metanetx_current,
        'gene': gene,
        'enzyme': enzyme,
        'rhea': rhea,
        'bigg': bigg,
        'kegg': kegg,
        'metacyc': metacyc,
        'reactome': reactome,
        'sabiork': sabiork,
        'seed': seed}

"""Extracts information from source about metabolites
arguments:
    metabolites_source (list<dict>): source information about metabolites
returns:
    (dict<dict>): information about metabolites\
"""
def extract_metabolites(
        metabolites_source=None):

    metabolites = {}

    for metabolite_source in metabolites_source:

        record = extract_metabolite(
            metabolite_source=metabolite_source)
        metabolites[record['identifier']] = record

    return metabolites

"""Extracts information from source about a metabolite
arguments:
    metabolite_source (dict): source information about a metabolite
returns:
    (dict): information about a metabolite
"""
def extract_metabolite(
        metabolite_source=None):

    # Determine information
    identifier = metabolite_source['identifier']
    name = metabolite_source['name']
    formula = metabolite_source['formula']
    mass = metabolite_source['mass']
    charge = metabolite_source['charge']
    references = extract_metabolite_references(
        identifier=identifier,
        references_source=metabolite_source['references'])

    return {
        'identifier': identifier,
        'name': name,
        'formula': formula,
        'mass': mass,
        'charge': charge,
        'references': references}

"""Extracts references from source about a metabolite
arguments:
    identifier (str): identifier of metabolite in MetaNetX
    references_source (str): source information about a metabolite's
        references
returns:
    (dict<str>): references about a metabolite
"""
def extract_metabolite_references(
        identifier=None,
        references_source=None):

    # Collect identifiers for each reference.
    metanetx_prior = extract_reference_information(
        references_source=references_source, key='deprecated:')

    if not ('MNXMK' in identifier):

        metanetx = [identifier] + metanetx_prior

    else:

        metanetx = metanetx_prior

    hmdb = extract_reference_information(
        references_source=references_source,
        key='hmdb:')
    pubchem = extract_reference_information(
        references_source=references_source,
        key='pubchem:')
    kegg = extract_reference_information(
        references_source=references_source,
        key='kegg:')
    chebi = extract_reference_information(
        references_source=references_source,
        key='chebi:')
    bigg = extract_reference_information(
        references_source=references_source,
        key='bigg:')
    envipath = extract_reference_information(
        references_source=references_source,
        key='envipath:')
    lipidmaps = extract_reference_information(
        references_source=references_source,
        key='lipidmaps:')
    metacyc = extract_reference_information(
        references_source=references_source,
        key='metacyc:')
    reactome = extract_reference_information(
        references_source=references_source,
        key='reactome:')
    sabiork = extract_reference_information(
        references_source=references_source,
        key='sabiork:')
    seed = extract_reference_information(
        references_source=references_source,
        key='seed:')
    slm = extract_reference_information(
        references_source=references_source,
        key='slm:')

    return {
        'metanetx': metanetx,
        'hmdb': hmdb,
        'pubchem': pubchem,
        'kegg': kegg,
        'chebi': chebi,
        'bigg': bigg,
        'metanetx': metanetx,
        'envipath': envipath,
        'lipidmaps': lipidmaps,
        'metacyc': metacyc,
        'reactome': reactome,
        'sabiork': sabiork,
        'seed': seed,
        'slm': slm}

"""Prepares report of information about metabolites for review.
arguments:
    metabolites (dict<dict>): information about metabolites
returns:
    (list<dict>): information about metabolites
"""
def prepare_report_metabolites(
        metabolites=None):

    records = []

    for metabolite in metabolites.values():

        record = {
            'identifier': metabolite['identifier'],
            'name': metabolite['name'],
            'formula': metabolite['formula'],
            'mass': metabolite['mass'],
            'charge': metabolite['charge'],
            'reference_metanetx': metabolite['references']['metanetx'],
            'reference_hmdb': metabolite['references']['hmdb'],
            'reference_chebi': metabolite['references']['chebi'],
            'reference_bigg': metabolite['references']['bigg'],
            'reference_kegg': metabolite['references']['kegg'],
            'reference_metacyc': metabolite['references']['metacyc'],
            'reference_reactome': metabolite['references']['reactome'],
            'reference_lipidmaps': metabolite['references']['lipidmaps'],
            'reference_sabiork': metabolite['references']['sabiork'],
            'reference_seed': metabolite['references']['seed'],
            'reference_slm': metabolite['references']['slm'],
            'reference_envipath': metabolite['references']['envipath']}
        records.append(record)

    return records

"""Prepares report of information about reactions for review.
arguments:
    reactions (dict<dict>): information about reactions
returns:
    (list<dict>): information about reactions
"""
def prepare_report_reactions(
        reactions=None):

    records = []

    for reaction in reactions.values():

        compartments = collect_value_from_records(
            key='compartment',
            records=reaction['participants'])
        metabolites = collect_value_from_records(
            key='metabolite',
            records=reaction['participants'])
        record = {
            'identifier': reaction['identifier'],
            'name': reaction['name'],
            'reversibility': reaction['reversibility'],
            'metabolites': metabolites,
            'compartments': compartments,
            'processes': reaction['processes'],
            'reference_gene': reaction['references']['gene'],
            'reference_metanetx': reaction['references']['metanetx']}
        records.append(record)

    return records

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
    path_compartments = directory + 'compartments.pickle'
    path_processes = directory + 'processes.pickle'
    path_reactions = directory + 'reactions.pickle'
    path_metabolites = directory + 'metabolites.pickle'
    path_metabolites_report = directory + 'metabolites.tsv'
    path_reactions_report = directory + 'reactions.tsv'

    # Write information to file.
    with open(path_compartments, 'wb') as file_product:
        pickle.dump(information['compartments'], file_product)
    with open(path_processes, 'wb') as file_product:
        pickle.dump(information['processes'], file_product)
    with open(path_reactions, 'wb') as file_product:
        pickle.dump(information['reactions'], file_product)
    with open(path_metabolites, 'wb') as file_product:
        pickle.dump(information['metabolites'], file_product)

    write_file_table(
        information=information['metabolites_report'],
        path_file=path_metabolites_report,
        names=information['metabolites_report'][0].keys(),
        delimiter='\t')
    write_file_table(
        information=information['reactions_report'],
        path_file=path_reactions_report,
        names=information['reactions_report'][0].keys(),
        delimiter='\t')

"""Run database collection
The purpose of this procedure is to extract information about metabolic
entities and sets from MetaNetX.
arguments:
    directory (str): path to directory for source and product files
"""
def execute_procedure(
        args_dict):

    # Read source information from file.
    source = read_source(
        directory=args_dict['source'])

    # Extract information about compartments.
    compartments = extract_compartments(
        compartments_source=source['compartments'])

    # Extract information about processes.
    processes = extract_processes(
        reactions_source=source['reactions'])

    # Extract reactions' names from metabolic model.
    reactions_names = extract_reactions_names(
        model=source['model'])

    # Extract information about reactions.
    reactions = extract_reactions(
        reactions_source=source['reactions'],
        reactions_names=reactions_names,
        genes_source=source['genes'],
        processes=processes)

    # Extract information about metabolites.
    metabolites = extract_metabolites(
        metabolites_source=source['metabolites'])

    # Prepare reports of information for review.
    metabolites_report = prepare_report_metabolites(
        metabolites=metabolites)
    reactions_report = prepare_report_reactions(
        reactions=reactions)

    # Compile information.
    information = {
        'compartments': compartments,
        'processes': processes,
        'reactions': reactions,
        'metabolites': metabolites,
        'metabolites_report': metabolites_report,
        'reactions_report': reactions_report,
        'reactions_names': list(reactions_names.values())}

    #Write product information to file
    write_product(
        args_dict['collect'],
        information=information)

    # Report.
    report = prepare_curation_report(
        compartments=compartments,
        processes=processes,
        reactions=reactions,
        metabolites=metabolites)

    print(report)
