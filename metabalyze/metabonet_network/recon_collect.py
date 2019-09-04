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
import pickle
import xml.etree.ElementTree as et

"""Import internal dependencies
"""
from metabalyze.metabonet_network.utils import confirm_path_directory
from metabalyze.metabonet_network.utils import prepare_curation_report
from metabalyze.metabonet_network.utils import read_file_table
from metabalyze.metabonet_network.utils import write_file_table
from metabalyze.metabonet_network.utils import collect_value_from_records
from metabalyze.metabonet_network.utils import find_match # originally called "find" in metabonet
from metabalyze.metabonet_network.utils import get_recon_references

"""Set globals
"""
# Source files
source_recon = 'recon'
source_metanetx = 'metanetx'
recon_reconciled_file = 'recon_reconciled.xml'

# Gene files
gene_file = 'enzymes.tsv'
gene_reaction = 'reaction'
gene_genes = 'genes'
genes_list = [
    'reaction',
    'genes',
    'low_bound',
    'up_bound',
    'direction']
gene_key = 'gene:'

# Compartment files
compartment_file = 'compartments.tsv'
compartment_id = 'identifier'
compartment_name = 'name'
compartment_list = [
    compartment_id,
    compartment_name,
    'source']

# Metabolite files
metabolite_file = 'chemicals.tsv'
metabolite_id = 'identifier'
metabolite_name = 'name'
metabolite_formula = 'formula'
metabolite_mass = 'mass'
metabolite_charge = 'charge'
metabolite_references = 'references'
metabolite_list = [
    metabolite_id,
    metabolite_name,
    'source',
    metabolite_formula,
    metabolite_mass,
    metabolite_charge,
    metabolite_references]
metabolite_metanetx_key = 'deprecated:'
metabolite_metanetx_id = 'MNXMK'
metabolite_hmdb_key = 'hmdb:'
metabolite_pubchem_key = 'pubchem:'
metabolite_kegg_key = 'kegg:'
metabolite_chebi_key = 'chebi:'
metabolite_bigg_key = 'bigg:'
metabolite_envipath_key = 'envipath:'
metabolite_lipidmaps_key = 'lipidmaps:'
metabolite_metacyc_key = 'metacyc:'
metabolite_reactome_key = 'reactome:'
metabolite_sabiork_key = 'sabiork:'
metabolite_seed_key = 'seed:'
metabolite_slm_key = 'slm:'

# Reaction files
reaction_file = 'reactions.tsv'
reaction_id = 'identifier'
reaction_equation = 'equation'
reaction_recon = 'recon'
reaction_metanetx = 'metanetx'
reaction_enzyme = 'enzyme_commission'
reaction_processes = 'processes'
reaction_references = 'references'
reaction_list = [
    reaction_id,
    reaction_equation,
    reaction_recon,
    reaction_metanetx,
    reaction_enzyme,
    reaction_processes,
    reaction_references]
reaction_finder = 'version:reaction'
reaction_key = 'model:'
reaction_rhea_key = 'rhea:'
reaction_bigg_key = 'bigg:'
reaction_metanetx_key = 'deprecated:'
reaction_kegg_key = 'kegg:'
reaction_metacyc_key = 'metacyc:'
reaction_reactome_key = 'reactome:'
reaction_sabiork_key = 'sabiork:'
reaction_seed_key = 'seed:'

# Other variables
id_attribute = 'id'
name_attribute = 'name'
equation_symbol = '<==>'
forward_symbol = '-->'
backward_symbol = '<--'
reactants_label = 'reactants'
reactant_label = 'reactant'
products_label = 'products'
product_label = 'product'
metabolite_compartment_raw_split_symbol = ' '
metabolite_compartment_split_symbol = '@'

"""Read data
Reads and organizes source information from file
returns:
    (dict): source information
"""
def read_source(
        directory):

    # Specify directories and files
    path_model = directory + recon_reconciled_file
    path_genes = directory + gene_file
    path_compartments = directory + compartment_file
    path_metabolites = directory + metabolite_file
    path_reactions = directory + reaction_file

    # Read information from file
    model = et.parse(path_model)

    genes = read_file_table(
        path_file=path_genes,
        names=genes_list,
        delimiter='\t')

    compartments = read_file_table(
        path_file=path_compartments,
        names=compartment_list,
        delimiter='\t')

    metabolites = read_file_table(
        path_file=path_metabolites,
        names=metabolite_list,
        delimiter='\t')

    reactions = read_file_table(
        path_file=path_reactions,
        names=reaction_list,
        delimiter='\t')

    return {
        'model': model,
        'compartments': compartments,
        'genes': genes,
        'reactions': reactions,
        'metabolites': metabolites}

"""Extracts information from source about compartments
arguments:
    compartments_source (list<dict>): source information about compartments
returns:
    (dict<dict>): information about compartments
"""
def extract_compartments(
        compartments_source):

    compartments = {}

    for compartment in compartments_source:

        record = {
            'identifier': compartment[compartment_id],
            'name': compartment[compartment_name]}

        compartments[compartment[compartment_id]] = record

    return compartments

"""Extracts information from source about processes
arguments:
    reactions_source (list<dict>): source information about reactions
returns:
    (dict<dict>): information about processes
"""
def extract_processes(
        reactions_source):

    processes = {}

    for reaction in reactions_source:

        reaction_processes_names = extract_reaction_processes_names(
            reaction_source=reaction)

        for name in reaction_processes_names:

            # Determine whether a record exists for the process
            novel = determine_process_name_novelty(
                name=name,
                processes=processes)

            if novel == True:

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
        reaction_source):

    # Separate references
    processes_source = reaction_source['processes']
    processes = extract_reference_information(
        key=reaction_key,
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
        key,
        references_source):

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
        name,
        processes):

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
        reference):

    reactions_names = {}

    reactions = reference['reactions'].findall(
        reaction_finder,
        reference['space'])

    for reaction in reactions:

        identifier = reaction.attrib[id_attribute]
        name = reaction.attrib[name_attribute]
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
        reactions_source,
        reactions_names,
        genes_source,
        processes):

    reactions = {}

    for reaction_source in reactions_source:

        record = extract_reaction(
            reaction_source=reaction_source,
            reactions_names=reactions_names,
            genes_source=genes_source,
            processes=processes)

        reactions[record[reaction_id]] = record

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
        reaction_source,
        reactions_names,
        genes_source,
        processes):

    # Determine information
    identifier = reaction_source[reaction_id]

    name = match_reaction_name(
        reaction_source=reaction_source,
        reactions_names=reactions_names)
    equation = reaction_source[reaction_equation]
    reversibility = extract_reaction_reversibility(
        equation=equation)
    participants = extract_reaction_participants(
        equation=equation)
    processes = extract_reaction_processes(
        reaction_source=reaction_source,
        processes=processes)
    references = extract_reaction_references(
        identifier=identifier,
        recon=reaction_source[reaction_recon],
        metanetx=reaction_source[reaction_metanetx],
        enzyme_commission=reaction_source[reaction_enzyme],
        genes=genes_source,
        references_source=reaction_source[reaction_references])

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
        reaction_source,
        reactions_names):

    # Reconciliation merges some multiple reactions from Recon to a single
    # reaction in MetaNetX.
    identifiers = reaction_source[reaction_recon]
    identifiers_split = identifiers.split(';')
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
        equation):

    if equation_symbol in equation:
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
        equation):

    # Extract raw information about reaction's participants
    participants_raw = extract_reaction_equation_raw_participants_by_role(
        equation=equation)

    # Extract information about participants' role, coefficient, metabolite,
    # and compartment
    reactants = extract_reaction_participants_by_role(
        participants_raw=participants_raw[reactants_label],
        role=reactant_label)
    products = extract_reaction_participants_by_role(
        participants_raw=participants_raw[products_label],
        role=product_label)

    return reactants + products

"""Extracts raw information about a reaction's participants from its equation
arguments:
    equation (str): a reaction's equation from MetaNetX
returns:
    (dict<list<str>>): raw information about a reaction's participants by
        role
"""
def extract_reaction_equation_raw_participants_by_role(
        equation):

    # Separate information about participants' reactants from products
    # Determine reaction's directionality
    if equation_symbol in equation:
        equation_sides = equation.split(' ' + equation_symbol + ' ')
        reactants_side = equation_sides[0]
        products_side = equation_sides[1]

    elif forward_symbol in equation:
        equation_sides = equation.split(' ' + forward_symbol + ' ')
        reactants_side = equation_sides[0]
        products_side = equation_sides[1]

    elif backward_symbol in equation:
        equation_sides = equation.split(' ' + backward_symbol + ' ')
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
        participants_raw,
        role):

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
        participant_raw,
        role):

    # Separate information
    participant_split = participant_raw.split(metabolite_compartment_raw_split_symbol)
    coefficient = float(participant_split[0])
    metabolite_compartment = participant_split[1].split(metabolite_compartment_split_symbol)
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
        identifier,
        genes_source):

    def match_reaction_gene(
            gene_record):

        return gene_record[gene_reaction] == identifier

    gene_source = find_match(
        match_reaction_gene,
        genes_source)
    genes_references = gene_source[gene_genes]
    genes_split_one = genes_references.split(';')
    genes_split_two = genes_references.split('+')

    genes = []

    for gene_pair in genes_split_two:

        identifier = gene_pair.replace(gene_key, '')

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
        reaction_source,
        processes):

    # Find process
    def match_reaction_process(
            record):

        return record['name'] == name

    reaction_processes_names = extract_reaction_processes_names(
        reaction_source=reaction_source)

    reaction_processes = []

    for name in reaction_processes_names:

        process_record = find_match(
            match_reaction_process,
            processes.values())
        identifier = process_record[reaction_id]
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
        identifier,
        recon,
        metanetx,
        enzyme_commission,
        genes,
        references_source):

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
        key=reaction_rhea_key)
    bigg = extract_reference_information(
        references_source=references_source,
        key=reaction_bigg_key)
    metanetx_prior = extract_reference_information(
        references_source=references_source,
        key=reaction_metanetx_key)

    if len(metanetx) > 0:
        metanetx_current = [metanetx] + metanetx_prior

    else:
        metanetx_current = metanetx_prior

    # Extract KEGG info
    kegg = extract_reference_information(
        references_source=references_source,
        key=reaction_kegg_key)
    metacyc = extract_reference_information(
        references_source=references_source,
        key=reaction_metacyc_key)
    reactome = extract_reference_information(
        references_source=references_source,
        key=reaction_reactome_key)
    sabiork = extract_reference_information(
        references_source=references_source,
        key=reaction_sabiork_key)
    seed = extract_reference_information(
        references_source=references_source,
        key=reaction_seed_key)

    return {
        'recon': [recon],
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
        metabolites_source):

    metabolites = {}

    for metabolite_source in metabolites_source:

        record = extract_metabolite(
            metabolite_source=metabolite_source)
        metabolites[record[metabolite_id]] = record

    return metabolites

"""Extracts information from source about a metabolite
arguments:
    metabolite_source (dict): source information about a metabolite
returns:
    (dict): information about a metabolite
"""
def extract_metabolite(
        metabolite_source):

    # Determine information
    identifier = metabolite_source[metabolite_id]
    name = metabolite_source[metabolite_name]
    formula = metabolite_source[metabolite_formula]
    mass = metabolite_source[metabolite_mass]
    charge = metabolite_source[metabolite_charge]
    references = extract_metabolite_references(
        identifier=identifier,
        references_source=metabolite_source[metabolite_references])

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
        identifier,
        references_source):

    # Collect identifiers for each reference.
    metanetx_prior = extract_reference_information(
        references_source=references_source,
        key=metabolite_metanetx_key)

    if not (metabolite_metanetx_id in identifier):
        metanetx = [identifier] + metanetx_prior

    else:
        metanetx = metanetx_prior

    hmdb = extract_reference_information(
        references_source=references_source,
        key=metabolite_hmdb_key)
    pubchem = extract_reference_information(
        references_source=references_source,
        key=metabolite_pubchem_key)
    kegg = extract_reference_information(
        references_source=references_source,
        key=metabolite_kegg_key)
    chebi = extract_reference_information(
        references_source=references_source,
        key=metabolite_chebi_key)
    bigg = extract_reference_information(
        references_source=references_source,
        key=metabolite_bigg_key)
    envipath = extract_reference_information(
        references_source=references_source,
        key=metabolite_envipath_key)
    lipidmaps = extract_reference_information(
        references_source=references_source,
        key=metabolite_lipidmaps_key)
    metacyc = extract_reference_information(
        references_source=references_source,
        key=metabolite_metacyc_key)
    reactome = extract_reference_information(
        references_source=references_source,
        key=metabolite_reactome_key)
    sabiork = extract_reference_information(
        references_source=references_source,
        key=metabolite_sabiork_key)
    seed = extract_reference_information(
        references_source=references_source,
        key=metabolite_seed_key)
    slm = extract_reference_information(
        references_source=references_source,
        key=metabolite_slm_key)

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
        metabolites):

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
        reactions):

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
        information):

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
def __main__(
        args_dict):

    # Read source information from file
    model = read_source(
        directory=args_dict['reconcile'])

    # Read in recon reference
    reference = get_recon_references(
        reference=model['model'])

    # Extract information about compartments
    compartments = extract_compartments(
        compartments_source=model['compartments'])

    # Extract information about processes.
    processes = extract_processes(
        reactions_source=model['reactions'])

    # Extract reactions' names from metabolic model.
    reactions_names = extract_reactions_names(
        reference=reference)

    # Extract information about reactions.
    reactions = extract_reactions(
        reactions_source=model['reactions'],
        reactions_names=reactions_names,
        genes_source=model['genes'],
        processes=processes)

    # Extract information about metabolites.
    metabolites = extract_metabolites(
        metabolites_source=model['metabolites'])

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
