"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
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
import argparse
from textwrap import dedent

"""Import internal dependencies
"""
from __init__ import __version__
from __init__ import __dependencies__
from utils import check_directories
from utils import check_curate
from utils import check_analyze
from utils import argument_checks

__path__  =  os.path.dirname(os.path.realpath(__file__))
url = 'https://raw.githubusercontent.com/j-berg/Metaboverse/master/metaboverse/__init__.py'

description_table  =  """\
    The Metaboverse sub-modules can be accessed by executing:
        'metaboverse sub-module_name arg1 arg2 ...'
    Sub-module help can be displayed by executing:
    'metaboverse sub-module_name --help'
    Sub-module descriptions:
        +-----------------------+--------------------------------------------+
        |    preprocess         |   Preprocess a dataframe                   |
        +-----------------------+--------------------------------------------+
        |    curate             |   Curate network with optional user data   |
        +-----------------------+--------------------------------------------+
"""

"""Check arguments
"""
def check_arguments(
        args_dict):

    # Run general checks
    args_dict = argument_checks(args_dict)

    # Run sub-module specific checks
    if args_dict['cmd'] == 'curate':
        check_curate(args_dict)

    elif args_dict['cmd'] == 'preprocess':
        check_preprocess(args_dict)

    elif args_dict['cmd'] == 'analyze':
        check_analyze(args_dict)

    else:
        raise Exception('Invalid sub-module selected')

    if 'output' not in args_dict \
    or args_dict['output'] == None \
    or args_dict['output'].lower() == 'none':
        args_dict['output'] = args_dict['output_file'].rsplit('/', 1)[0] + '/'



    # Print out user commands to log file
    print('Metaboverse version: ' + str(__version__))
    print('======================\nUser commands summary:\n======================')

    for key, value in args_dict.items():
        print(str(key) + ': ' + str(value))
    print('=====================\nEnd commands summary\n=====================\n')

    return args_dict

"""Parse arguments
Will print help menu if no arguments are provided
"""
def parse_arguments(
        args,
        __version__):

    # Require user input
    if args is None:
        args = sys.argv[1:]

    # Initialize main parser
    parser = argparse.ArgumentParser(
        prog = 'metaboverse',
        description = dedent(description_table),
        formatter_class = argparse.RawTextHelpFormatter)

    # Optional main arguments
    parser.add_argument(
        '-v', '--version',
        help = 'Print installed version to stout',
        action = 'version',
        version = '%(prog)s ' + str(__version__))

    # Add sub-parsers
    subparser = parser.add_subparsers(dest = 'cmd')

    # Preprocess parser
    preprocess_parser = subparser.add_parser(
        'preprocess',
        description = 'Preprocess a dataframe',
        add_help = False)

    # Curate required arguments
    preprocess_reqs = preprocess_parser.add_argument_group('required arguments')
    preprocess_reqs.add_argument(
        '-d', '--data',
        help = 'Path and filename of dataset -- refer to documentation for details on formatting',
        metavar = '<path/filename>',
        type = str,
        required = True)
    preprocess_reqs.add_argument(
        '-m', '--metadata',
        help = 'Path and filename of metadata -- refer to documentation for details on formatting',
        metavar = '<path/filename>',
        type = str,
        required = True)
    preprocess_reqs.add_argument(
        '-t', '--type',
        help = 'Omics type (options: transcriptomics, proteomics, metabolomics)',
        metavar = '<transcriptomics/proteomics/metabolomics>',
        type = str,
        required = True)

    # Curate parser
    curate_parser = subparser.add_parser(
        'curate',
        description = 'Curate biological network',
        add_help = False)

    # Curate required arguments
    curate_reqs = curate_parser.add_argument_group('required arguments')
    curate_reqs.add_argument(
        '-o', '--output',
        help = 'Path to output directory (default: current working directory)',
        metavar = '<path>',
        type = str,
        required = True)
    curate_reqs.add_argument(
        '-s', '--species_id',
        help = 'Reactome species ID',
        metavar = '<species_id>',
        type = str,
        default = 'HSA',
        required = True)

    # Curate optional arguments
    curate_opts = curate_parser.add_argument_group('optional arguments')
    curate_opts.add_argument(
        '-h', '--help',
        action = 'help',
        help = 'Show help message and exit')
    curate_opts.add_argument(
        '-c', '--organism_curation',
        help = 'Path and name for organism curation file',
        metavar = '<path/filename.pickle>',
        type = str,
        required = False)
    curate_opts.add_argument(
        '-i', '--model_file',
        help = 'Path and name for organism curation file output. If --organism_curation is used, this argument will be ignored.',
        metavar = '<path/filename.pickle>',
        type = str,
        required = False)
    curate_opts.add_argument(
        '-f', '--output_file',
        help = 'Path and name for output database file',
        metavar = '<path/filename.json>',
        type = str,
        required = False)
    curate_opts.add_argument(
        '-r', '--transcriptomics',
        help = 'Path and filename of RNA-Seq data -- refer to documentation for details on formatting and normalization',
        metavar = '<path/filename>',
        type = str,
        default = 'None',
        required = False)
    curate_opts.add_argument(
        '-p', '--proteomics',
        help = 'Path and filename of proteomics data -- refer to documentation for details on formatting and normalization',
        metavar = '<path/filename>',
        type = str,
        default = 'None',
        required = False)
    curate_opts.add_argument(
        '-m', '--metabolomics',
        help = 'Path and filename of metabolomics data -- refer to documentation for details on formatting and normalization',
        metavar = '<path/filename>',
        type = str,
        default = 'None',
        required = False)
    #curate_opts.add_argument(
    #    '-a', '--additional_reactions',
    #    help = 'Path and filename of additional reaction table. See #documentation for more details on appropriate file formatting.',
    #    metavar = '<file.txt, file.tsv>',
    #    type = str,
    #    default = 'None',
    #    required = False)
    curate_opts.add_argument(
        '--collapse_with_modifiers',
        help = 'Include modifiers when considering a potential reaction collapse.',
        action = 'store_true',
        required = False)
    curate_opts.add_argument(
        '--broadcast_genes',
        help = 'Broadcast gene values to proteins where protein values not available.',
        action = 'store_true',
        required = False)
    curate_opts.add_argument(
        '-e', '--experiment_type',
        help = 'Specify experiment type',
        metavar = '<default/timecourse/flux/multi-condition>',
        type = str,
        required = False)
    curate_opts.add_argument(
        '--experiment_name',
        help = 'Specify experiment name',
        metavar = '<default/timecourse/flux/multi-condition>',
        type = str,
        required = False)
    curate_opts.add_argument(
        '--labels',
        help = 'Comma separated list of labels to use for multi-condition or timecourse experiments',
        metavar = 'cond1, cond2, cond3, ...',
        type = str,
        required = False)
    curate_opts.add_argument(
        '--blacklist',
        help = 'Comma separated list of names for metabolites, etc., to ignore in the analysis and visualization.',
        metavar = 'H+, H2O, etc...',
        type = str,
        required = False)
    curate_opts.add_argument(
        '-n', '--session_data',
        help = 'Path and filename to session data file',
        metavar = '<path/filename>',
        type = str,
        required = False)
    curate_opts.add_argument(
        '-l', '--progress_log',
        help = 'Path and filename to progress log file',
        metavar = '<path/filename>',
        type = str,
        required = False)

    # Get arguments are print help if no arguments provided
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse arguments into NameSpace
    args = parser.parse_args(args)

    # Collect subargs and package, add metaboverse script path to parameter dictionary
    args_dict = vars(args)
    args_dict['path'] = str(__path__) + '/'

    # Check inputs validity
    args_dict = check_arguments(
        args_dict)

    return args, args_dict
