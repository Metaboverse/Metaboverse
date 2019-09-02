"""License Information
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

"""Import dependencies
"""

"""Curate network model
"""
def curate(
        args_dict):

    args_dict = set_paths(args_dict)

    # Move start files to args_dict['output'] sub-directory 'source' for downstream use

    return args_dict

"""Specify output sub-directory paths
"""
def set_paths(
        args_dict):

    args_dict['reconcile'] = args_dict['output'] + 'reconcile/'
    args_dict['collect'] = args_dict['output'] + 'collect/'
    args_dict['extract'] = args_dict['output'] + 'extract/'
    args_dict['enhance'] = args_dict['output'] + 'enhance/'
    args_dict['curate'] = args_dict['output'] + 'curate/'
    args_dict['convert'] = args_dict['output'] + 'convert/'
    args_dict['provision'] = args_dict['output'] + 'provision/'

    # Step 1

    # Step 2

    # Step 3

    # ...

    return args_dict













"""Notes
"""

"""

        The root directory must contain several sources of information before
        execution of the procedure.
        /root/source/
            recon2m2.xml
            hmdb_metabolites.xml
        /root/source/customization/
            curation_compartments.tsv
            curation_metabolites.tsv
            curation_processes.tsv
            curation_reactions.tsv
            filtration_compartments.tsv
            filtration_processes.tsv
            reconciliation_compartments.tsv
            reconciliation_metabolites.tsv
            simplification_metabolites.tsv
            simplification_reactions.tsv


        --------------------------------------------------
        reconciliation
        Reconcile information in human metabolic model to MetaNetX.
        --------------------------------------------------
        collection
        Collect relevant information from metabolic model and from MetaNetX
        about compartments, processes, reactions, and metabolites.
        --------------------------------------------------
        extraction
        Extract relevant information about metabolites from entries in Human
        Metabolome Database (HMDB).
        --------------------------------------------------
        enhancement
        Enhance information about metabolites and reactions.
        --------------------------------------------------
        curation
        Curate information about compartments, processes, reactions, and
        metabolites.
        --------------------------------------------------
        conversion
        Convert information about compartments, processes, reactions, and
        metabolites to versatile formats.
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

"""
