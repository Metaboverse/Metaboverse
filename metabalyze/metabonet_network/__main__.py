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
import sys
import shutil

"""Import internal dependencies
"""
from metabalyze.metabonet_network.recon_reconcile import __main__ as run_reconcile
from metabalyze.metabonet_network.recon_collect import __main__ as run_collect
from metabalyze.metabonet_network.hmdb_extract import __main__ as run_extract

"""Set globals
"""
reconciliation_files = [
    'compartments.tsv',
    'enzymes.tsv', # genes
    'chemicals.tsv', # metabolites
    'reactions.tsv']

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

"""Remove intermediate network model files and directories
"""
def remove_intermediates(
        args_dict):

    # Move final files to args_dict['output']
    input = ''
    shutil.move(input, args_dict['output'] + 'curated_model.xml')

    # Delete all intermediate folders
    shutil.rmtree(args_dict['reconcile'])
    shutil.rmtree(args_dict['collect'])
    shutil.rmtree(args_dict['extract'])
    shutil.rmtree(args_dict['enhance'])
    shutil.rmtree(args_dict['curate'])
    shutil.rmtree(args_dict['convert'])
    shutil.rmtree(args_dict['provision'])

"""Check all members of file list exist
"""
def check_file_list(
        path,
        files):

    missing = False

    for f in files:
        if not os.path.isfile(path + f):
            print('! ' + f + ' not found')
            missing = True

    return missing



"""Check reconciliation files from metanetx have been added
"""
def check_reconciliation_success(
        path,
        files):

    # Print instructions
    print('')
    print('--> Recon reconiliation file for MetaNetX curation at:\n')
    print(path)
    print('└── recon_reconciled.xml\n')
    print('--> Please upload this file to https://www.metanetx.org/cgi-bin/mnxweb/import_mnet and run Import')
    print('--> Output files can be found by clicking the link under \"Mnet\" and downloading the relevant files')
    print('--> MetaNetX output files must then be loaded as follows:')
    print('')
    print(path)

    # Print reconciliation file tree structure
    for x in reconciliation_files:
        print('└── ' + x)

    print('')
    print('--> Access the mapping summary by clicking on the Import -> Mapping summary links and save as follows:')
    print('')
    print(path)
    print('└── mapping_summary.txt')

    # Get user input to continue
    input('\n--> Press enter to continue once these files have been processed and added to the specified location')
    print('')
    print('Confirming directory contents...')
    os.system('tree ' + path)
    input('\n--> Press enter to continue if the directory contents are correct')

    # Check files exist
    files_missing = check_file_list(
        path,
        files)

    if files_missing == True:
        print('--> Add missing files to continue')
        sys.exit(1)

"""Step 0 -- reconcile
!- Maybe consider trying to automatically load this to metanetx through some kind of API and saving output files
"""
def _reconcile(
        args_dict):

    print('\nPreparing Recon metabolic model for MetaNetX reconciliation...\n')

    run_reconcile(
        args_dict)

    # Check metanetx reconciliation files added
    check_reconciliation_success(
        args_dict['reconcile'],
        reconciliation_files)

    print('\n=== Reconciliation complete ===\n')

"""Step 1 -- collect
"""
def _collect(
        args_dict):

    print('Collecting relevant information from reconciled Recon metabolic model about compartments, processes, reactions, and metabolites...\n')

    run_collect(
        args_dict)

    print('=== Collection complete ===\n')

"""Step 2 -- extract
"""
def _extract(
        args_dict):

    print('Extracting relevant information about metabolites from HMDB...\n')

    run_extract(
        args_dict)

    print('=== Extraction complete ===\n')

"""Step 3 -- enhance
"""
def _enhance(
        args_dict):

    print('Enhancing information about metabolites and reactions from Recon and HMDB...\n')

"""Step 4 -- curate
"""
def _curate(
        args_dict):

    print('Curating information about compartments, processes, reactions, and metabolites...\n')

"""Step 5 -- convert
"""
def _convert(
        args_dict):

    print('Converting information about compartments, processes, reactions, and metabolites to versatile formats...\n')

"""Curate network model
Required inputs:
    - recon.xml
    - hmdb_metabolites.xml
"""
def __main__(
        args_dict):

    args_dict = set_paths(args_dict)

    # Reconcile
    files_missing = check_file_list(
        args_dict['reconcile'],
        reconciliation_files)

    run_reconcile = False

    if files_missing == False:
        print('\nIt appears you already have files that were reconciled using MetaNetX\n')
        os.system('tree ' + args_dict['reconcile'])
        user_input = input('\n--> Do you want to use these reconciliation files? (Y/n) ')

        if user_input.lower() == 'n':
            _reconcile(args_dict)

        else:
            print('')

    elif 'recon' in args_dict \
    and args_dict['recon'] != None:
        _reconcile(args_dict)

    else:
        print('Skipping reconciliation, no Recon files found...')

    # Collect
    if 'recon' in args_dict \
    and args_dict['recon'] != None:
        _collect(args_dict)

    # Extract
    if 'hmdb' in args_dict \
    and args_dict['hmdb'] != None:
        _extract(args_dict)

    # Enhance
    _enhance(args_dict)

    # Curate
    _curate(args_dict)

    # Convert
    _convert(args_dict)

    # Provision?

    # Clean up intermediate files
    #remove_intermediates(args_dict)

    sys.exit(1)
