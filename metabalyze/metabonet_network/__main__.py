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

"""Import internal dependencies
"""
from metabalyze.metabonet_network.reconcile import __main__ as run_reconcile

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

    return args_dict

"""Step 0 -- reconcile
"""
def _reconcile(
        args_dict):

    print('\nReconciling information in human metabolic model...')

    run_reconcile(
        args_dict)

    # Print instructions
    print('')
    print('--> Recon reconiliation for MetaNetX found at: ' + args_dict['reconcile'] + 'recon_reconciled.xml')
    print('--> Please upload this file to https://www.metanetx.org/cgi-bin/mnxweb/import_mnet and run Import')
    print('--> Output files need to be loaded into ' + args_dict['reconcile'] + ' with the prefix \"metanetx\"')

    print('=== Reconiliation complete ===')

"""Step 1 -- collect
"""
def _collect(args_dict):

    print('Collecting relevant information from metabolic model about compartments, processes, reactions, and metabolites...')

"""Step 2 -- extract
"""
def _extract():

    print('Extracting relevant information about metabolites...')

"""Step 3 -- enhance
"""
def _enhance():

    print('Enhancing information about metabolites and reactions...')

"""Step 4 -- curate
"""
def _curate():

    print('Curating information about compartments, processes, reactions, and metabolites...')

"""Step 5 -- convert
"""
def _convert():

    print('Converting information about compartments, processes, reactions, and metabolites to versatile formats...')

"""Curate network model
Required inputs:
    - recon2m2.xml
    - hmdb_metabolites.xml
"""
def __main__(
        args_dict):

    args_dict = set_paths(args_dict)

    # Move start files to args_dict['output'] sub-directory 'source' for downstream use
    _reconcile(args_dict)
    _collect(args_dict)

    return args_dict
