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
import argparse
from textwrap import dedent

"""Import internal dependencies
"""
from metabalyze.__init__ import __version__
from metabalyze.__init__ import __dependencies__
from metabalyze.utils import check_directories
from metabalyze.utils import check_curate
from metabalyze.utils import check_analyze
from metabalyze.utils import generate_log
from metabalyze.utils import argument_checks

"""Set global variables
"""
DEFAULT_MAX_PROCESSORS = None
DEFAULT_HUB_STRINGENCY = 50

__path__  =  os.path.dirname(os.path.realpath(__file__))
url = 'https://raw.githubusercontent.com/j-berg/MetaboNet-Analyzer/master/metabalyzer/__init__.py'

description_table  =  """\
    The metabonet-analyzer sub-modules can be accessed by executing:
        'metabalyzer sub-module_name arg1 arg2 ...'
    Sub-module help can be displayed by executing:
    'metabalyzer sub-module_name --help'
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
        'echo \"MetaboNet-Analyzer version: ' + str(__version__) + '\"'
        + str(args_dict['log']))
    print('MetaboNet-Analyzer version: ' + str(__version__))

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
    curate_reqs.add_argument(
        '--recon',
        help = 'Path and filename of Recon database',
        metavar = '<path/filename>',
        type = str,
        required = True)
    curate_reqs.add_argument(
        '--hmdb',
        help = 'Path and filename of HMDB database',
        metavar = '<path/filename>',
        type = str,
        required = True)

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
        '-c', '--component',
        help = 'Provide argument if you wish to only select the main network component for output',
        action = 'store_true',
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
        help = 'Path to metabonet network model for NetworkX (should be at output_dir/_network/network.pickle)',
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
        '-r', '--rnaseq',
        help = 'Path and filename of RNA-Seq data -- refer to documentation for details on formatting',
        metavar = '<path/filename>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '-p', '--proteomics',
        help = 'Path and filename of proteomics data -- refer to documentation for details on formatting',
        metavar = '<path/filename>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '-b', '--metabolomics',
        help = 'Path and filename of metabolomics data -- refer to documentation for details on formatting',
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
        '--pathway',
        help = 'Name of pathway to analyze',
        metavar = '<pathway name>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '-n', '--normalization',
        help = 'Normalization to perform on data',
        metavar = '<type>',
        type = str,
        required = False)
    analyze_opts.add_argument(
        '--compartments',
        help = 'Model network with compartments',
        action = 'store_true',
        required = False)
    analyze_opts.add_argument(
        '--hubs',
        help = 'Model network with hubs',
        action = 'store_true',
        required = False)
    analyze_opts.add_argument(
        '--hub_stringency',
        help = 'Change connectedness degree for hub cut-off (default: %s)' % DEFAULT_HUB_STRINGENCY,
        default = DEFAULT_HUB_STRINGENCY,
        type = int,
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

    # Collect subargs and package, add metabonet-analyzer script path to parameter dictionary
    args_dict = vars(args)
    args_dict['path'] = str(__path__) + '/'

    # Check inputs validity
    args_dict = check_arguments(
        args_dict)

    return args, args_dict
