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
from metaboverse.__init__ import __version__
from metaboverse.__init__ import __dependencies__
from metaboverse.utils import check_directories
from metaboverse.utils import check_curate
from metaboverse.utils import check_analyze
from metaboverse.utils import generate_log
from metaboverse.utils import argument_checks

"""Set global variables
"""
DEFAULT_MAX_PROCESSORS = None
DEFAULT_BLACKLIST = [
    'H+', # proton
    'H2O', # water
    'O2', # dioxygen
    # phosphate
    # diphosphate
    'CO2', # carbon dioxide
    # sulfate
    # hydrogen peroxide
    # ammonium
    # sulfite
    # sodium
    # hydrogen carbonate
    # hydroxide
    ]

__path__  =  os.path.dirname(os.path.realpath(__file__))
url = 'https://raw.githubusercontent.com/j-berg/Metaboverse/master/metaboverse/__init__.py'

description_table  =  """\
    The Metaboverse sub-modules can be accessed by executing:
        'metaboverse sub-module_name arg1 arg2 ...'
    Sub-module help can be displayed by executing:
    'metaboverse sub-module_name --help'
    Sub-module descriptions:
        +-----------------------+--------------------------------------------------------------------------------------+
        |    curate             |   Curate network model                                                               |
        |-----------------------|--------------------------------------------------------------------------------------|
        |    preprocess         |   Preprocess -omics data                                                             |
        |-----------------------|--------------------------------------------------------------------------------------|
        |    analyze            |   Analyze -omics data using network model                                            |
        +-----------------------+--------------------------------------------------------------------------------------+
"""

"""Check dependencies
"""
def get_dependencies(
        args_dict,
        dependencies):

    # Make dependencies log
    os.system(
        'echo \"Conda dependencies:\n===================\" >> '
        + str(args_dict['log_loc']) + 'dependencies.log')

    for d in dependencies:
        os.system(
            'pip freeze | grep -i \"' + str(d) + '\" >> '
            + str(args_dict['log_loc']) + 'dependencies.log')

    os.system(
        'echo \"R dependencies:\n===================\" >> '
        + str(args_dict['log_loc']) + 'dependencies.log')

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

    args_dict = generate_log(args_dict)

    # Print out user commands to log file
    os.system(
        'echo \"======================\nUser commands summary:\n======================\"'
        + str(args_dict['log']))
    print('======================\nUser commands summary:\n======================')

    os.system(
        'echo \"Metaboverse version: ' + str(__version__) + '\"'
        + str(args_dict['log']))
    print('Metaboverse version: ' + str(__version__))

    for key, value in args_dict.items():

        os.system(
            'echo \"' + str(key) + ': ' + str(value) + '\"'
            + str(args_dict['log']))
        print(str(key) + ': ' + str(value))

    os.system(
        'echo \"=====================\nEnd commands summary\n=====================\n\"'
        + str(args_dict['log']))
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
        prog = 'metabalyze',
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

    # Curate parser
    curate_parser = subparser.add_parser(
        'curate',
        description = 'Curate biological network',
        add_help = False)

    # Curate required arguments
    curate_reqs = curate_parser.add_argument_group('required arguments')

    # Curate optional arguments
    curate_opts = curate_parser.add_argument_group('optional arguments')
    curate_opts.add_argument(
        '-h', '--help',
        action = 'help',
        help = 'Show help message and exit')
    curate_opts.add_argument(
        '-o', '--output',
        help = 'Path to output directory (default: current working directory)',
        metavar = '<path>',
        type = str,
        required = False)
    curate_opts.add_argument(
        '-s', '--species',
        help = 'Reactome species ID',
        metavar = '<species_id>',
        type = str,
        default = 'HSA',
        required = False)
    curate_opts.add_argument(
        '-m', '--max_processors',
        help = 'Number of max processors to use for tasks (default: No limit)',
        metavar = '<processors>',
        type = int,
        default = DEFAULT_MAX_PROCESSORS,
        required = False)

    # Preprocess parser
    preprocess_parser = subparser.add_parser(
        'preprocess',
        description = 'Preprocess data for network analysis',
        add_help = False)

    # Preprocess required arguments
    preprocess_reqs = preprocess_parser.add_argument_group('required arguments')

    # Preprocess optional arguments
    preprocess_opts = preprocess_parser.add_argument_group('optional arguments')
    preprocess_opts.add_argument(
        '-h', '--help',
        action = 'help',
        help = 'Show help message and exit')

    # Analyze parser
    analyze_parser = subparser.add_parser(
        'analyze',
        description = 'Analyze data on biological network',
        add_help = False)

    # Analyze required arguments
    analyze_reqs = analyze_parser.add_argument_group('required arguments')
    analyze_reqs.add_argument(
        '-d', '--model',
        help = 'Path to metaboverse network model for NetworkX (should be at output_dir/_network/network.pickle)',
        metavar = '<path/filename>',
        type = str,
        required = True)
    analyze_reqs.add_argument(
        '-t', '--metadata',
        help = 'Path and filename of metadata -- refer to documentation for details on formatting',
        metavar = '<path/filename>',
        type = str,
        required = True)

    # Analyze optional arguments
    analyze_opts = analyze_parser.add_argument_group('optional arguments')
    analyze_opts.add_argument(
        '-h', '--help',
        action = 'help',
        help = 'Show help message and exit')
    analyze_opts.add_argument(
        '-s', '--species',
        help = 'Provide Reactome species ID (default: HSA for human)',
        metavar = '<id>',
        type = str,
        default = 'HSA',
        required = False)
    analyze_opts.add_argument(
        '-r', '--rnaseq',
        help = 'Path and filename of RNA-Seq data -- refer to documentation for details on formatting and normalization',
        metavar = '<path/filename>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '-p', '--proteomics',
        help = 'Path and filename of proteomics data -- refer to documentation for details on formatting and normalization',
        metavar = '<path/filename>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '-b', '--metabolomics',
        help = 'Path and filename of metabolomics data -- refer to documentation for details on formatting and normalization',
        metavar = '<path/filename>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '--time',
        help = 'Analyze timecourse data',
        action = 'store_true',
        required = False)
    analyze_opts.add_argument(
        '--static',
        help = 'Analyze static data (two conditions allowed)',
        action = 'store_true',
        required = False)
    analyze_opts.add_argument(
        '--normalize',
        help = 'Normalize expression values on standard scale (z-score), otherwise will display log$_2$(Fold change) values on plotting',
        action = 'store_true',
        required = False)
    analyze_opts.add_argument(
        '--pathway',
        help = 'Name of pathway to analyze. If not provided, will analyze all available pathways. Use the parameter PROVIDE to print a list of available pathways to analyze',
        metavar = '<pathway name>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '--blacklist',
        help = 'Provide space-seperated list of analytes to blacklist',
        default = DEFAULT_BLACKLIST,
        type = str,
        nargs = '+',
        required = False)
    analyze_opts.add_argument(
        '-m', '--max_processors',
        help = 'Number of max processors to use for tasks (default: No limit)',
        metavar = '<processors>',
        type = int,
        default = DEFAULT_MAX_PROCESSORS,
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
